#pragma once

#include <model.h>
#include <model_unstruct.h>
#include <sep_reader.h>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace model
{
/**
 * @class SepUnstructBuilder
 * @brief Builder for constructing unstructured spectral element models from SEP
 * files
 *
 * @tparam FloatType Numeric type for coordinates and properties (float or
 * double)
 * @tparam ScalarType Integer type for indices (int, int64_t, etc.)
 *
 * @details
 * Inherits from ModelBuilderBase and provides implementation for building
 * ModelUnstruct instances from Stanford Exploration Project (SEP) files. The
 * builder automates:
 * - Reading SEP file metadata and binary data
 * - Validating grid parameters
 * - Creating spectral element mesh
 * - Mapping element connectivity
 * - Loading velocity properties
 *
 * @see ModelBuilderBase, SepReader, ModelUnstruct
 */
template <typename FloatType, typename ScalarType>
class SepUnstructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  using ModelBuilderBase<FloatType, ScalarType>::MAX_ORDER;

  /// @brief Default constructor - requires setSepFile() before getModel()
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
      model_data_;            ///< Model data
  mutable SepHeader header_;  ///< SEP metadata

  /// Load and process SEP file
  void loadSepData() const
  {
    header_ = SepReader::readHeader(sepfile_);
    header_.print();

    int effective_order = order_;
    if (effective_order <= 0) effective_order = header_.order;
    if (effective_order <= 0 || effective_order > MAX_ORDER)
      effective_order = 4;

    std::vector<float> raw_data = SepReader::readData(header_);
    fillMeshData(raw_data, effective_order);
  }

  /// Build mesh from SEP data
  void fillMeshData(const std::vector<float>& raw_data, int order) const
  {
    const int n1 = header_.n1, n2 = header_.n2, n3 = header_.n3;
    int n_elem_i = (n1 - 1) / order;
    int n_elem_j = (n2 - 1) / order;
    int n_elem_k = (n3 - 1) / order;

    const int n_element = n_elem_i * n_elem_j * n_elem_k;
    const int n_node = n1 * n2 * n3;

    auto global_node_index = allocateArray2D<arrayInt>(
        n_element, (order + 1) * (order + 1) * (order + 1));
    auto nodes_coords_x = allocateVector<vectorReal>(n_node);
    auto nodes_coords_y = allocateVector<vectorReal>(n_node);
    auto nodes_coords_z = allocateVector<vectorReal>(n_node);
    auto model_vp_node = allocateVector<vectorReal>(n_node);

    generateCoordinates(nodes_coords_x, nodes_coords_y, nodes_coords_z, n1, n2,
                        n3);
    generateConnectivity(global_node_index, model_vp_node, raw_data, n_element,
                         n_node, n1, n2, n3, order, n_elem_i, n_elem_j,
                         n_elem_k);

    FloatType lx = (n1 - 1) * header_.d1;
    FloatType ly = (n2 - 1) * header_.d2;
    FloatType lz = (n3 - 1) * header_.d3;

    auto empty = allocateVector<vectorReal>(0);

    model_data_ = model::ModelUnstructData<FloatType, ScalarType>(
        order, n_element, n_node, lx, ly, lz, true, false, global_node_index,
        nodes_coords_x, nodes_coords_y, nodes_coords_z, model_vp_node, empty,
        empty, empty, empty, empty, empty, empty, empty, empty, empty, empty,
        empty, empty, empty, empty, empty);
  }

  /// Generate Cartesian node coordinates
  void generateCoordinates(vectorReal& coords_x, vectorReal& coords_y,
                           vectorReal& coords_z, int n1, int n2, int n3) const
  {
    int idx = 0;
    for (int k = 0; k < n3; ++k)
      for (int j = 0; j < n2; ++j)
        for (int i = 0; i < n1; ++i)
        {
          coords_x[idx] = header_.o1 + i * header_.d1;
          coords_y[idx] = header_.o2 + j * header_.d2;
          coords_z[idx] = header_.o3 + k * header_.d3;
          idx++;
        }
  }

  /// Build element-to-node connectivity
  void generateConnectivity(arrayInt& global_node_index, vectorReal& model_vp,
                            const std::vector<float>& raw_data, int n_element,
                            int n_node, int n1, int n2, int n3, int order,
                            int n_elem_i, int n_elem_j, int n_elem_k) const
  {
    for (int i = 0; i < n_node && i < (int)raw_data.size(); ++i)
      model_vp[i] = static_cast<FloatType>(raw_data[i]);

    int elem_idx = 0;
    for (int elem_k = 0; elem_k < n_elem_k; ++elem_k)
      for (int elem_j = 0; elem_j < n_elem_j; ++elem_j)
        for (int elem_i = 0; elem_i < n_elem_i; ++elem_i)
        {
          int start_i = elem_i * order, start_j = elem_j * order,
              start_k = elem_k * order;
          int local_idx = 0;
          for (int k = 0; k <= order; ++k)
            for (int j = 0; j <= order; ++j)
              for (int i = 0; i <= order; ++i)
              {
                int gi = std::min(start_i + i, n1 - 1);
                int gj = std::min(start_j + j, n2 - 1);
                int gk = std::min(start_k + k, n3 - 1);
                global_node_index(elem_idx, local_idx++) =
                    gk * n1 * n2 + gj * n1 + gi;
              }
          if (++elem_idx >= n_element) break;
        }
  }
};

}  // namespace model
