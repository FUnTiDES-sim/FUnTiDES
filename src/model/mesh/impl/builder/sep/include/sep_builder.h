#pragma once

#include <model.h>
#include <model_unstruct.h>
#include <sep_reader.h>

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace model
{

/**
 * @class GllQuadrature
 * @brief Gauss-Lobatto-Legendre quadrature points and weights
 * @details Precomputed GLL points for polynomial orders 1-5
 */
class GllQuadrature
{
 public:
  /// Get GLL points in [-1, 1] for given order
  static const std::vector<double>& getGllPoints(int order)
  {
    static const std::vector<std::vector<double>> gll_points = {
        {},  // order 0 (unused)
        {-1.0, 1.0},
        {-1.0, 0.0, 1.0},
        {-1.0, -0.447213595499958, 0.447213595499958, 1.0},
        {-1.0, -0.654653670707977, 0.0, 0.654653670707977, 1.0},
        {-1.0, -0.764884927562167, -0.285231768779026, 0.285231768779026,
         0.764884927562167, 1.0}};
    if (order < 1 || order >= (int)gll_points.size())
      throw std::invalid_argument("GLL points only defined for orders 1-5");
    return gll_points[order];
  }
};

/**
 * @class SepUnstructBuilder
 * @brief Builder for constructing unstructured spectral element models from
 * SEP files
 *
 * @tparam FloatType Numeric type for coordinates and properties (float or
 * double)
 * @tparam ScalarType Integer type for indices (int, int64_t, etc.)
 *
 * @details
 * Implements proper SEM workflow:
 * 1. Read SEP metadata and velocity data
 * 2. Create spectral elements from SEP grid
 * 3. Generate GLL nodes within each element using reference coordinates
 * 4. Map nodes to physical space via isoparametric mapping
 * 5. Interpolate velocity to GLL nodes
 *
 * @see ModelBuilderBase, SepReader, ModelUnstruct, GllQuadrature
 */
template <typename FloatType, typename ScalarType>
class SepUnstructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  using ModelBuilderBase<FloatType, ScalarType>::MAX_ORDER;

  /// @brief Default constructor
  SepUnstructBuilder() : sepfile_(""), order_(-1) {}

  /**
   * @brief Constructor with SEP file path and polynomial order
   * @param sepfile Path to the SEP file (.H or .H@)
   * @param order Polynomial order (1-5, or -1 to use header value)
   * @throws std::invalid_argument if order exceeds MAX_ORDER
   */
  SepUnstructBuilder(const std::string& sepfile, const int order)
      : sepfile_(sepfile), order_(order)
  {
    if (order > 0 && order > MAX_ORDER)
      throw std::invalid_argument("Polynomial order must be between 1 and " +
                                  std::to_string(MAX_ORDER));
  }

  /**
   * @brief Set the SEP file path
   * @param sepfile Path to the SEP file to read
   */
  void setSepFile(const std::string& sepfile) { sepfile_ = sepfile; }

  /**
   * @brief Set the polynomial order
   * @param order Polynomial order (1-5, -1 uses header value)
   * @throws std::invalid_argument if order exceeds MAX_ORDER
   */
  void setOrder(int order)
  {
    if (order > 0 && order > MAX_ORDER)
      throw std::invalid_argument("Polynomial order must be between 1 and " +
                                  std::to_string(MAX_ORDER));
    order_ = order;
  }

  /**
   * @brief Build and return the ModelUnstruct
   * @return std::shared_ptr to ready-to-use ModelUnstruct
   * @throws std::runtime_error if no SEP file configured
   * @throws std::runtime_error if SEP file cannot be read
   * @throws std::runtime_error if mesh generation fails
   */
  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const override
  {
    if (sepfile_.empty())
      throw std::runtime_error(
          "No SEP file configured. Call setSepFile() first.");

    loadSepData();
    return std::make_shared<model::ModelUnstruct<FloatType, ScalarType>>(
        model_data_);
  }

 private:
  mutable std::string sepfile_;  ///< Path to SEP file
  mutable int order_;            ///< Polynomial order
  mutable model::ModelUnstructData<FloatType, ScalarType>
      model_data_;                           ///< Model data
  mutable SepHeader header_;                 ///< SEP metadata
  mutable std::vector<float> raw_velocity_;  ///< Raw velocity grid

  /// Load SEP file and generate mesh
  void loadSepData() const
  {
    header_ = SepReader::readHeader(sepfile_);
    header_.print();

    int effective_order = order_;
    if (effective_order <= 0) effective_order = header_.order;
    if (effective_order <= 0 || effective_order > MAX_ORDER)
      effective_order = 4;

    raw_velocity_ = SepReader::readData(header_);
    buildSemMesh(effective_order);
  }

  /// Build SEM mesh: elements first, then GLL nodes
  void buildSemMesh(int order) const
  {
    // Compute number of spectral elements and nodes
    int n_elem_i = header_.n1;
    int n_node_i = (n_elem_i * order) + 1;
    int n_elem_j = header_.n2;
    int n_node_j = (n_elem_j * order) + 1;
    int n_elem_k = header_.n3;
    int n_node_k = (n_elem_k * order) + 1;
    const int n_element = n_elem_i * n_elem_j * n_elem_k;
    const int n_nodes_per_elem = (order + 1) * (order + 1) * (order + 1);
    const int n_nodes = n_node = (n_elem_i * order) + 1 * (n_elem_j * order) +
                                 1 * (n_elem_k * order) + 1;

    // GLL points in reference element [-1, 1]^3
    const auto& gll_pts = GllQuadrature::getGllPoints(order);

    // Storage for global node data
    // std::vector<FloatType> global_coords_x, global_coords_y, global_coords_z;
    auto global_coords_x =
        allocateVector<vectorReal>(n_nodes) std::vector<FloatType>
            global_velocity;
    auto global_node_index =
        allocateArray2D<arrayInt>(n_element, n_nodes_per_elem);

    // Iterate over spectral elements
    for (int elem_k = 0; elem_k < n_elem_k; ++elem_k)
    {
      for (int elem_j = 0; elem_j < n_elem_j; ++elem_j)
      {
        for (int elem_i = 0; elem_i < n_elem_i; ++elem_i)
        {
          int elem_idx =
              elem_k * n_elem_i * n_elem_j + elem_j * n_elem_i + elem_i;

          // Element corners in SEP grid
          int corner_i_start = elem_i * order;
          int corner_j_start = elem_j * order;
          int corner_k_start = elem_k * order;
          int corner_i_end = std::min(corner_i_start + order, n_node_i - 1);
          int corner_j_end = std::min(corner_j_start + order, n_node_j - 1);
          int corner_k_end = std::min(corner_k_start + order, n_node_k - 1);

          // Physical coordinates of element corners
          FloatType x_min = header_.o1 + corner_i_start * header_.d1;
          FloatType x_max = x_min + header_.d1;
          FloatType y_min = header_.o2 + corner_j_start * header_.d2;
          FloatType y_max = y_min + header_.d2;
          FloatType z_min = header_.o3 + corner_k_start * header_.d3;
          FloatType z_max = z_max + header_.d3;

          // Generate GLL nodes for this element
          int local_idx = 0;
          for (int k_idx = 0; k_idx <= order; ++k_idx)
          {
            for (int j_idx = 0; j_idx <= order; ++j_idx)
            {
              for (int i_idx = 0; i_idx <= order; ++i_idx)
              {
                auto node_key = std::make_tuple(elem_i, elem_j, elem_k, i_idx,
                                                j_idx, k_idx);

                ScalarType node_id;
                // generate GLL coordinates and interpolate velocity

                // Map reference element [-1,1]^3 to physical coordinates
                double xi = gll_pts[i_idx];
                double eta = gll_pts[j_idx];
                double zeta = gll_pts[k_idx];

                // Isoparametric mapping: bilinear in each direction
                FloatType x_phys = x_min + 0.5 * (xi + 1.0) * (x_max - x_min);
                FloatType y_phys = y_min + 0.5 * (eta + 1.0) * (y_max - y_min);
                FloatType z_phys = z_min + 0.5 * (zeta + 1.0) * (z_max - z_min);

                // Store the coordinate
                global_coords_x.push_back(x_phys);
                global_coords_y.push_back(y_phys);
                global_coords_z.push_back(z_phys);

                // Interpolate velocity at GLL node
                FloatType velocity =
                    interpolateVelocity(x_phys, y_phys, z_phys);
                global_velocity.push_back(velocity);

                node_id = global_node_id;
                node_map[node_key] = global_node_id++;

                global_node_index(elem_idx, local_idx++) = node_id;
              }
            }
          }
        }
      }
    }

    // Store in model data
    ScalarType n_node_scalar = static_cast<ScalarType>(global_node_id);
    FloatType lx = (n_node_i - 1) * header_.d1;
    FloatType ly = (n_node_j - 1) * header_.d2;
    FloatType lz = (n_node_k - 1) * header_.d3;

    auto coords_x = allocateVector<vectorReal>(n_node_scalar);
    auto coords_y = allocateVector<vectorReal>(n_node_scalar);
    auto coords_z = allocateVector<vectorReal>(n_node_scalar);
    auto velocity = allocateVector<vectorReal>(n_node_scalar);

    auto empty = allocateVector<vectorReal>(0);

    model_data_ = model::ModelUnstructData<FloatType, ScalarType>(
        static_cast<ScalarType>(order), static_cast<ScalarType>(n_element),
        n_node_scalar, lx, ly, lz, true, false, global_node_index, coords_x,
        coords_y, coords_z,
        velocity,  // model_vp_node
        empty,     // model_vp_element
        empty,     // model_rho_node
        empty,     // model_rho_element
        empty,     // model_vs_node
        empty,     // model_vs_element
        empty,     // model_delta_node
        empty,     // model_delta_element
        empty,     // model_epsilon_node
        empty,     // model_epsilon_element
        empty,     // model_gamma_node
        empty,     // model_gamma_element
        empty,     // model_theta_node
        empty,     // model_theta_element
        empty,     // model_phi_node
        empty,     // model_phi_element
        empty);    // boundaries_t
  }

  /// Trilinear interpolation of velocity on the SEP grid
  FloatType interpolateVelocity(FloatType x, FloatType y, FloatType z) const
  {
    // Map physical coordinates to SEP grid indices
    double i_real = (x - header_.o1) / header_.d1;
    double j_real = (y - header_.o2) / header_.d2;
    double k_real = (z - header_.o3) / header_.d3;

    // Clamp to grid bounds
    i_real = std::max(0.0, std::min(i_real, (double)(header_.n1 - 1)));
    j_real = std::max(0.0, std::min(j_real, (double)(header_.n2 - 1)));
    k_real = std::max(0.0, std::min(k_real, (double)(header_.n3 - 1)));

    int i0 = (int)i_real, i1 = std::min(i0 + 1, header_.n1 - 1);
    int j0 = (int)j_real, j1 = std::min(j0 + 1, header_.n2 - 1);
    int k0 = (int)k_real, k1 = std::min(k0 + 1, header_.n3 - 1);

    double di = i_real - i0;
    double dj = j_real - j0;
    double dk = k_real - k0;

    // Trilinear interpolation
    auto getVel = [this](int i, int j, int k) {
      int idx = k * header_.n1 * header_.n2 + j * header_.n1 + i;
      if (idx >= 0 && idx < (int)raw_velocity_.size())
        return (FloatType)raw_velocity_[idx];
      return FloatType(0);
    };

    FloatType v000 = getVel(i0, j0, k0);
    FloatType v100 = getVel(i1, j0, k0);
    FloatType v010 = getVel(i0, j1, k0);
    FloatType v110 = getVel(i1, j1, k0);
    FloatType v001 = getVel(i0, j0, k1);
    FloatType v101 = getVel(i1, j0, k1);
    FloatType v011 = getVel(i0, j1, k1);
    FloatType v111 = getVel(i1, j1, k1);

    FloatType v00 = v000 * (1.0 - di) + v100 * di;
    FloatType v10 = v010 * (1.0 - di) + v110 * di;
    FloatType v01 = v001 * (1.0 - di) + v101 * di;
    FloatType v11 = v011 * (1.0 - di) + v111 * di;

    FloatType v0 = v00 * (1.0 - dj) + v10 * dj;
    FloatType v1 = v01 * (1.0 - dj) + v11 * dj;

    return v0 * (1.0 - dk) + v1 * dk;
  }
};

}  // namespace model
