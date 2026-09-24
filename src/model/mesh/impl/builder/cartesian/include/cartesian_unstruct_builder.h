#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BUILDER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BUILDER_H_

#include <stdexcept>

#include "builder.h"
#include "cartesian_model_file_reader.h"
#include "cartesian_params.h"
#include "cartesian_unstruct_boundary_classifier.h"
#include "gllpoints.h"
#include "model_unstruct.h"

namespace model {
/**
 * @brief Builds an unstructured-storage model on a Cartesian grid of hexahedral elements.
 *
 * The constructor generates the element-to-global-node table, the GLL node coordinates and the
 * material arrays of the local subdomain. getModel() then assembles the model with boundary flags
 * classified against the global domain.
 *
 * @tparam FloatType Floating-point type of the model data.
 * @tparam ScalarType Integer type of the model indices.
 */
template <typename FloatType, typename ScalarType>
class CartesianUnstructBuilder : public ModelBuilderBase<FloatType, ScalarType> {
 public:
  using ModelBuilderBase<FloatType, ScalarType>::MAX_ORDER;

  /// @brief Default constructor: builds nothing, members other than the origins are left uninitialized.
  CartesianUnstructBuilder() {}

  /**
   * @brief Builds the node table, node coordinates and material arrays of the subdomain.
   *
   * If a global size is not positive, the global size and origin of that axis fall back to the
   * local ones.
   *
   * @param[in] p Local subdomain description, global domain description and material options.
   * @throws std::runtime_error If the order is outside [1, MAX_GLL_ORDER], or if a model file is
   *         given with node-based models or does not match the element count.
   */
  CartesianUnstructBuilder(const CartesianParams<FloatType, ScalarType>& p)
      : ex_(p.ex),
        ey_(p.ey),
        ez_(p.ez),
        lx_(p.lx),
        ly_(p.ly),
        lz_(p.lz),
        order_(p.order),
        isModelOnNodes_(p.isModelOnNodes),
        isElastic_(p.isElastic),
        isAcoustoElastic_(p.isAcoustoElastic),
        acoustoElasticBoundaryZ_(p.acoustoElasticBoundaryZ),
        DgSemBoundaryZ_(p.DgSemBoundaryZ),
        ox_(p.origin_x),
        oy_(p.origin_y),
        oz_(p.origin_z),
        global_lx_(p.global_lx > 0 ? p.global_lx : p.lx),
        global_ly_(p.global_ly > 0 ? p.global_ly : p.ly),
        global_lz_(p.global_lz > 0 ? p.global_lz : p.lz),
        global_ox_(p.global_lx > 0 ? p.global_origin_x : p.origin_x),
        global_oy_(p.global_ly > 0 ? p.global_origin_y : p.origin_y),
        global_oz_(p.global_lz > 0 ? p.global_origin_z : p.origin_z),
        model_file_(p.model_file) {
    initGlobalNodeList();
    initNodesCoords();
    initModels();
  }

  /**
   * @brief Assembles the model of the local subdomain.
   *
   * @param[in] free_surface_on_top Forwarded to the boundary classifier, which decides how the
   *            top face of the global domain is flagged.
   * @return Model holding the local mesh, the material arrays, the boundary flags and the face
   *         connectivity.
   */
  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel(bool free_surface_on_top) const override {
    const int n_node = (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);

    // Tolerance for the boundary classifier: 1e-4 of the smallest element size.
    const FloatType tol = std::min({lx_ / ex_, ly_ / ey_, lz_ / ez_}) * static_cast<FloatType>(1e-4);

    // Boundaries are classified against the global domain, not the local subdomain.
    auto boundaries_t = CartesianUnstructBoundaryClassifier<FloatType, ScalarType>(
                            global_ox_, global_ox_ + global_lx_, global_oy_, global_oy_ + global_ly_, global_oz_,
                            global_oz_ + global_lz_, tol, free_surface_on_top)
                            .classify(n_node, nodes_coords_x_, nodes_coords_y_, nodes_coords_z_);

    model::ModelUnstructData<FloatType, ScalarType> modelData(
        order_, ex_ * ey_ * ez_, n_node, lx_, ly_, lz_, isModelOnNodes_, isElastic_, global_node_index_,
        nodes_coords_x_, nodes_coords_y_, nodes_coords_z_, model_vp_node_, model_vp_element_, model_rho_node_,
        model_rho_element_, model_vs_node_, model_vs_element_, model_delta_node_, model_delta_element_,
        model_epsilon_node_, model_epsilon_element_, model_gamma_node_, model_gamma_element_, model_theta_node_,
        model_theta_element_, model_phi_node_, model_phi_element_, model_C_tensor_element_, boundaries_t);

    modelData.ox_ = ox_;
    modelData.oy_ = oy_;
    modelData.oz_ = oz_;

    auto model = std::make_shared<model::ModelUnstruct<FloatType, ScalarType>>(modelData);
    model->buildFaceConnectivity();

    return model;
  }

  ~CartesianUnstructBuilder() = default;

 private:
  FloatType ox_{0}, oy_{0}, oz_{0};                       ///< Local subdomain origin.
  FloatType global_ox_{0}, global_oy_{0}, global_oz_{0};  ///< Global domain origin.
  ScalarType ex_, ey_, ez_;                               ///< Number of local elements per axis.
  FloatType lx_, ly_, lz_;                                ///< Local subdomain size per axis.
  FloatType global_lx_{0}, global_ly_{0}, global_lz_{0};  ///< Global domain size per axis.

  int order_;
  bool isModelOnNodes_;  ///< True: material stored per node; false: per element.
  bool isElastic_;
  bool isAcoustoElastic_{false};
  FloatType acoustoElasticBoundaryZ_{static_cast<FloatType>(0)};  ///< z at and above which the medium is fluid.
  FloatType DgSemBoundaryZ_{static_cast<FloatType>(0)};
  std::string model_file_;

  arrayInt global_node_index_;  ///< Size numElements x (order+1)^3, global node of each local node.
  vectorReal nodes_coords_x_;   ///< Size total node count, indexed by global node number.
  vectorReal nodes_coords_y_;
  vectorReal nodes_coords_z_;

  vectorReal model_vp_node_;
  vectorReal model_vp_element_;
  vectorReal model_rho_node_;
  vectorReal model_rho_element_;
  vectorReal model_vs_node_;
  vectorReal model_vs_element_;
  vectorReal model_delta_node_;
  vectorReal model_delta_element_;
  vectorReal model_epsilon_node_;
  vectorReal model_epsilon_element_;
  vectorReal model_gamma_node_;
  vectorReal model_gamma_element_;
  vectorReal model_theta_node_;
  vectorReal model_theta_element_;
  vectorReal model_phi_node_;
  vectorReal model_phi_element_;
  array3DReal model_C_tensor_element_;

  /// Fills global_node_index_. Elements are numbered i + j*ex + k*ex*ey, local nodes l + n*(order+1) + m*(order+1)^2
  /// (x fastest), global nodes i + j*nx + k*nx*ny with nx = ex*order+1.
  void initGlobalNodeList() {
    int nodes_x = order_ + 1;
    int nodes_y = order_ + 1;
    int nodes_z = order_ + 1;
    int total_nodes = nodes_x * nodes_y * nodes_z;
    global_node_index_ = allocateArray2D<arrayInt>(ex_ * ey_ * ez_, total_nodes, "global node index");
    int nx = ex_ * order_ + 1;
    int ny = ey_ * order_ + 1;
    int nz = ez_ * order_ + 1;

    for (int k = 0; k < ez_; k++) {
      for (int j = 0; j < ey_; j++) {
        for (int i = 0; i < ex_; i++) {
          int elementNum = i + j * ex_ + k * ex_ * ey_;
          int offset = i * order_ + j * order_ * nx + k * order_ * nx * ny;

          for (int m = 0; m < order_ + 1; m++) {
            for (int n = 0; n < order_ + 1; n++) {
              for (int l = 0; l < order_ + 1; l++) {
                int dofLocal = l + n * (order_ + 1) + m * (order_ + 1) * (order_ + 1);
                int dofGlobal = offset + l + n * nx + m * nx * ny;
                global_node_index_(elementNum, dofLocal) = dofGlobal;
              }
            }
          }
        }
      }
    }
  }

  /// Writes the order+1 GLL node coordinates of element n_element along one axis into coord.
  /// h is the element size along that axis.
  void getCoordInOneDirection(FloatType h, const int& n_element, float* coord, FloatType offset) {
    if (order_ < 1 || order_ > MAX_GLL_ORDER)
      throw std::runtime_error("Cartesian unstruct builder error: order not supported.");

    const FloatType elementStart = n_element * h;
    for (int j = 0; j < order_ + 1; j++) {
      coord[j] = elementStart + (GLLPoints::get(order_, j) + 1.0f) * h * 0.5f + offset;
    }
  }

  /// Fills nodes_coords_{x,y,z}_ with GLL node coordinates, without the subdomain origin.
  void initNodesCoords() {
    int nodes_x = ex_ * order_ + 1;
    int nodes_y = ey_ * order_ + 1;
    int nodes_z = ez_ * order_ + 1;
    int total_nodes = nodes_x * nodes_y * nodes_z;

    nodes_coords_x_ = allocateVector<vectorReal>(total_nodes, "nodes coords x");
    nodes_coords_y_ = allocateVector<vectorReal>(total_nodes, "nodes coords y");
    nodes_coords_z_ = allocateVector<vectorReal>(total_nodes, "nodes coords z");

    float coord_x[MAX_ORDER + 1];
    float coord_y[MAX_ORDER + 1];
    float coord_z[MAX_ORDER + 1];

    auto hx = lx_ / ex_;
    auto hy = ly_ / ey_;
    auto hz = lz_ / ez_;

    for (int n = 0; n < ez_; n++) {
      getCoordInOneDirection(hz, n, coord_z, 0);
      for (int m = 0; m < ey_; m++) {
        getCoordInOneDirection(hy, m, coord_y, 0);
        for (int l = 0; l < ex_; l++) {
          getCoordInOneDirection(hx, l, coord_x, 0);

          for (int k = 0; k < order_ + 1; k++) {
            for (int j = 0; j < order_ + 1; j++) {
              for (int i = 0; i < order_ + 1; i++) {
                int global_i = l * order_ + i;
                int global_j = m * order_ + j;
                int global_k = n * order_ + k;

                int global_node_index = global_i + global_j * nodes_x + global_k * nodes_x * nodes_y;

                if (global_i < nodes_x && global_j < nodes_y && global_k < nodes_z) {
                  nodes_coords_x_(global_node_index) = coord_x[i];
                  nodes_coords_y_(global_node_index) = coord_y[j];
                  nodes_coords_z_(global_node_index) = coord_z[k];
                }
              }
            }
          }
        }
      }
    }
  }

  /// Fills the material arrays (per node or per element) with uniform or two-layer values, then
  /// overrides per-element properties from model_file_ if one is given.
  void initModels() {
    // TODO: only uniform and two-layer models are generated; no user-defined material
    // variation other than the element-based model file.
    int n_element = ex_ * ey_ * ez_;
    int n_node = (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);
    if (isModelOnNodes_) {
      model_rho_node_ = allocateVector<vectorReal>(n_node, "model rho node");
      model_vp_node_ = allocateVector<vectorReal>(n_node, "model vp node");

      if (isAcoustoElastic_) {
        // Two layers: fluid for z >= boundary, solid below. The coupled solver classifies
        // elements by vs: vs = 0 means fluid.
        model_vs_node_ = allocateVector<vectorReal>(n_node, "model vs node");
        model_delta_node_ = allocateVector<vectorReal>(n_node, "model delta node");
        model_gamma_node_ = allocateVector<vectorReal>(n_node, "model gamma node");
        model_epsilon_node_ = allocateVector<vectorReal>(n_node, "model epsilon node");
        model_theta_node_ = allocateVector<vectorReal>(n_node, "model theta node");
        model_phi_node_ = allocateVector<vectorReal>(n_node, "model phi node");

        for (int n = 0; n < n_node; ++n) {
          FloatType z = nodes_coords_z_(n);
          bool const is_fluid = (z >= acoustoElasticBoundaryZ_);
          model_rho_node_[n] = is_fluid ? static_cast<FloatType>(1020) : static_cast<FloatType>(2500);
          model_vp_node_[n] = is_fluid ? static_cast<FloatType>(1500) : static_cast<FloatType>(3000);
          model_vs_node_[n] = is_fluid ? static_cast<FloatType>(0) : static_cast<FloatType>(1500);
          model_delta_node_[n] = 0;
          model_epsilon_node_[n] = 0;
          model_gamma_node_[n] = 0;
          model_theta_node_[n] = 0;
          model_phi_node_[n] = 0;
        }
      } else {
        for (int i = 0; i < n_node; i++) {
          model_rho_node_[i] = 1;
          model_vp_node_[i] = 1500;
        }
        if (isElastic_) {
          model_vs_node_ = allocateVector<vectorReal>(n_node, "model vs node");
          model_delta_node_ = allocateVector<vectorReal>(n_node, "model delta node");
          model_gamma_node_ = allocateVector<vectorReal>(n_node, "model gamma node");
          model_epsilon_node_ = allocateVector<vectorReal>(n_node, "model epsilon node");
          model_theta_node_ = allocateVector<vectorReal>(n_node, "model theta node");
          model_phi_node_ = allocateVector<vectorReal>(n_node, "model phi node");

          for (int i = 0; i < n_node; i++) {
            model_vs_node_[i] = 755;
            model_delta_node_[i] = 0.0;
            model_epsilon_node_[i] = 0.0;
            model_gamma_node_[i] = 0.0;
            model_theta_node_[i] = 0.0;
            model_phi_node_[i] = 0.0;
          }
        }
      }
    }

    else {
      model_rho_element_ = allocateVector<vectorReal>(n_element, "model rho elem");
      model_vp_element_ = allocateVector<vectorReal>(n_element, "model vp elem");

      if (isAcoustoElastic_) {
        // Two layers: fluid when the element centroid z >= boundary, solid below. The coupled
        // solver classifies elements by vs: vs = 0 means fluid.
        model_vs_element_ = allocateVector<vectorReal>(n_element, "model vs element");
        model_delta_element_ = allocateVector<vectorReal>(n_element, "model delta element");
        model_gamma_element_ = allocateVector<vectorReal>(n_element, "model gamma element");
        model_epsilon_element_ = allocateVector<vectorReal>(n_element, "model epsilon element");
        model_theta_element_ = allocateVector<vectorReal>(n_element, "model theta element");
        model_phi_element_ = allocateVector<vectorReal>(n_element, "model phi element");

        FloatType const hz = lz_ / ez_;
        for (int k = 0; k < ez_; ++k) {
          FloatType const centroid_z = oz_ + (k + static_cast<FloatType>(0.5)) * hz;
          bool const is_fluid = (centroid_z >= acoustoElasticBoundaryZ_);
          for (int j = 0; j < ey_; ++j) {
            for (int i = 0; i < ex_; ++i) {
              int const e = i + j * ex_ + k * ex_ * ey_;
              model_rho_element_[e] = is_fluid ? static_cast<FloatType>(1000) : static_cast<FloatType>(2000);
              model_vp_element_[e] = is_fluid ? static_cast<FloatType>(1500) : static_cast<FloatType>(3000);
              model_vs_element_[e] = is_fluid ? static_cast<FloatType>(0) : static_cast<FloatType>(1500);
              model_delta_element_[e] = 0;
              model_epsilon_element_[e] = 0;
              model_gamma_element_[e] = 0;
              model_theta_element_[e] = 0;
              model_phi_element_[e] = 0;
            }
          }
        }
      } else {
        for (int i = 0; i < n_element; i++) {
          model_rho_element_[i] = 1;
          model_vp_element_[i] = 1500;
        }

        if (isElastic_) {
          model_vs_element_ = allocateVector<vectorReal>(n_element, "model vs element");
          model_delta_element_ = allocateVector<vectorReal>(n_element, "model delta element");
          model_gamma_element_ = allocateVector<vectorReal>(n_element, "model gamma element");
          model_epsilon_element_ = allocateVector<vectorReal>(n_element, "model epsilon element");
          model_theta_element_ = allocateVector<vectorReal>(n_element, "model theta element");
          model_phi_element_ = allocateVector<vectorReal>(n_element, "model phi element");

          for (int i = 0; i < n_element; i++) {
            model_vs_element_[i] = 755;
            model_delta_element_[i] = 0.0;
            model_epsilon_element_[i] = 0.0;
            model_gamma_element_[i] = 0.0;
            model_theta_element_[i] = 0.0;
            model_phi_element_[i] = 0.0;
          }
        }
      }
    }
    if (!model_file_.empty()) {
      if (isModelOnNodes_) {
        throw std::runtime_error(
            "[CartesianUnstructBuilder] model_file only supported with "
            "isModelOnNodes=false for now.");
      }
      model::CartesianModelFileReader reader(model_file_);
      const size_t n = reader.count();
      auto fill_view = [&](vectorReal& view, const std::string& prop, const std::string& label) {
        if (!reader.has(prop)) return;
        view = allocateVector<vectorReal>(static_cast<int>(n), label.c_str());
        const auto& buf = reader.get(prop);
        auto host = Kokkos::create_mirror_view(view);
        for (size_t i = 0; i < n; ++i) host[i] = static_cast<FloatType>(buf[i]);
        Kokkos::deep_copy(view, host);
      };
      fill_view(model_vp_element_, "Vp", "model_vp_element");
      fill_view(model_rho_element_, "Rho", "model_rho_element");
      fill_view(model_vs_element_, "Vs", "model_vs_element");
      fill_view(model_delta_element_, "Delta", "model_delta_element");
      fill_view(model_epsilon_element_, "Epsilon", "model_epsilon_element");
      fill_view(model_gamma_element_, "Gamma", "model_gamma_element");
      fill_view(model_theta_element_, "Theta", "model_theta_element");
      fill_view(model_phi_element_, "Phi", "model_phi_element");

      std::cout << "[Model] vp[0]=" << model_vp_element_[0] << " vp[mid]=" << model_vp_element_[n_element / 2]
                << " vp[last]=" << model_vp_element_[n_element - 1] << std::endl;

      // The file must provide one value per element of the local mesh.
      if (model_vp_element_.extent(0) != static_cast<size_t>(n_element))
        throw std::runtime_error("[CartesianUnstructBuilder] model_file has " +
                                 std::to_string(model_vp_element_.extent(0)) + " elements but mesh has " +
                                 std::to_string(n_element));
    }
  }
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_BUILDER_H_
