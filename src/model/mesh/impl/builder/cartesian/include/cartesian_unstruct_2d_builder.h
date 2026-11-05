#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_2D_BUILDER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_2D_BUILDER_H_

#include <algorithm>
#include <limits>
#include <stdexcept>

#include "builder.h"
#include "cartesian_params.h"
#include "cartesian_unstruct_boundary_classifier.h"
#include "gllpoints.h"
#include "model_unstruct.h"

namespace model {

/**
 * @brief Builds a 2-D unstructured-storage model on a Cartesian grid of quadrilaterals in the
 * (x, z) plane (ModelUnstruct with dimension 2).
 *
 * Uses order, ex, ez, lx, lz, isModelOnNodes, origin_x/z and global_lx/lz (global_origin_x/z)
 * of CartesianParams; the y members are ignored. Acoustic only: isElastic, isAcoustoElastic and
 * model_file are not supported.
 *
 * Numbering: element i + k*ex; local node l + m*(order+1) (l along x, m along z); global node
 * gi + gk*nx with nx = ex*order+1. The node coordinates include the subdomain origin; the y
 * coordinates are 0.
 *
 * Material: vp = 1500 m/s and rho = 1 everywhere, per element or per node. Heterogeneous models
 * are given directly as ModelUnstructData arrays.
 *
 * @tparam FloatType Floating-point type of the model data.
 * @tparam ScalarType Integer type of the model indices.
 */
template <typename FloatType, typename ScalarType>
class CartesianUnstruct2DBuilder : public ModelBuilderBase<FloatType, ScalarType> {
 public:
  /**
   * @brief Builds the node table, node coordinates and material arrays.
   * @param[in] p Mesh description (see the class documentation for the members used).
   * @throws std::invalid_argument If the order is outside [1, MAX_GLL_ORDER], a size is not
   *         positive, or an unsupported option (elastic, acousto-elastic, model file) is set.
   */
  explicit CartesianUnstruct2DBuilder(const CartesianParams<FloatType, ScalarType>& p)
      : order_(p.order),
        ex_(p.ex),
        ez_(p.ez),
        lx_(p.lx),
        lz_(p.lz),
        ox_(p.origin_x),
        oz_(p.origin_z),
        global_lx_(p.global_lx > 0 ? p.global_lx : p.lx),
        global_lz_(p.global_lz > 0 ? p.global_lz : p.lz),
        global_ox_(p.global_lx > 0 ? p.global_origin_x : p.origin_x),
        global_oz_(p.global_lz > 0 ? p.global_origin_z : p.origin_z),
        isModelOnNodes_(p.isModelOnNodes) {
    if (order_ < 1 || order_ > MAX_GLL_ORDER)
      throw std::invalid_argument("CartesianUnstruct2DBuilder: order not supported");
    if (ex_ < 1 || ez_ < 1 || !(lx_ > 0) || !(lz_ > 0))
      throw std::invalid_argument("CartesianUnstruct2DBuilder: element counts and sizes must be positive");
    if (p.isElastic || p.isAcoustoElastic || !p.model_file.empty())
      throw std::invalid_argument("CartesianUnstruct2DBuilder: acoustic models only, without model file");
    initGlobalNodeList();
    initNodesCoords();
    initModels();
  }

  /**
   * @brief Assembles the 2-D model.
   * @param[in] free_surface_on_top true to flag the nodes of the global z-max edge as Surface,
   *            false to flag them as Damping like the other global boundary edges.
   * @return ModelUnstruct of dimension 2 (no face connectivity).
   */
  std::shared_ptr<ModelApi<FloatType, ScalarType>> getModel(bool free_surface_on_top) const override {
    FloatType const tol = std::min(lx_ / ex_, lz_ / ez_) * static_cast<FloatType>(1e-4);
    // No y boundary in 2-D: infinite y bounds are never reached.
    FloatType const inf = std::numeric_limits<FloatType>::infinity();
    auto boundaries_t = CartesianUnstructBoundaryClassifier<FloatType, ScalarType>(
                            global_ox_, global_ox_ + global_lx_, -inf, inf, global_oz_, global_oz_ + global_lz_,
                            tol, free_surface_on_top)
                            .classify(numberOfNodes(), nodes_coords_x_, nodes_coords_y_, nodes_coords_z_);

    ModelUnstructData<FloatType, ScalarType> data;
    data.dim_ = 2;
    data.order_ = order_;
    data.n_element_ = ex_ * ez_;
    data.n_node_ = numberOfNodes();
    data.lx_ = lx_;
    data.ly_ = 0;
    data.lz_ = lz_;
    data.ox_ = ox_;
    data.oy_ = 0;
    data.oz_ = oz_;
    data.isModelOnNodes_ = isModelOnNodes_;
    data.isElastic_ = false;
    data.global_node_index_ = global_node_index_;
    data.nodes_coords_x_ = nodes_coords_x_;
    data.nodes_coords_y_ = nodes_coords_y_;
    data.nodes_coords_z_ = nodes_coords_z_;
    data.model_vp_node_ = model_vp_node_;
    data.model_vp_element_ = model_vp_element_;
    data.model_rho_node_ = model_rho_node_;
    data.model_rho_element_ = model_rho_element_;
    data.boundaries_t_ = boundaries_t;
    return std::make_shared<ModelUnstruct<FloatType, ScalarType>>(data);
  }

 private:
  int order_;
  ScalarType ex_, ez_;
  FloatType lx_, lz_;
  FloatType ox_, oz_;
  FloatType global_lx_, global_lz_;
  FloatType global_ox_, global_oz_;
  bool isModelOnNodes_;

  dualArrayInt global_node_index_;
  dualVectorReal nodes_coords_x_;
  dualVectorReal nodes_coords_y_;
  dualVectorReal nodes_coords_z_;
  dualVectorReal model_vp_node_;
  dualVectorReal model_vp_element_;
  dualVectorReal model_rho_node_;
  dualVectorReal model_rho_element_;

  int nodesX() const { return ex_ * order_ + 1; }
  int nodesZ() const { return ez_ * order_ + 1; }
  int numberOfNodes() const { return nodesX() * nodesZ(); }

  /// Fills global_node_index_ (see the class documentation for the numbering).
  void initGlobalNodeList() {
    int const n1d = order_ + 1;
    int const nx = nodesX();
    global_node_index_ = allocateArray2D<dualArrayInt>(ex_ * ez_, n1d * n1d, "global node index");
    auto gni = global_node_index_.view_host();
    for (int k = 0; k < ez_; ++k)
      for (int i = 0; i < ex_; ++i) {
        int const e = i + k * ex_;
        for (int m = 0; m < n1d; ++m)
          for (int l = 0; l < n1d; ++l) gni(e, l + m * n1d) = (i * order_ + l) + (k * order_ + m) * nx;
      }
    global_node_index_.modify_host();
  }

  /// Fills the node coordinates: GLL points of each element along x and z, plus the origin.
  void initNodesCoords() {
    int const nx = nodesX();
    int const nz = nodesZ();
    nodes_coords_x_ = allocateVector<dualVectorReal>(nx * nz, "nodes coords x");
    nodes_coords_y_ = allocateVector<dualVectorReal>(nx * nz, "nodes coords y");
    nodes_coords_z_ = allocateVector<dualVectorReal>(nx * nz, "nodes coords z");
    auto x = nodes_coords_x_.view_host();
    auto z = nodes_coords_z_.view_host();
    FloatType const hx = lx_ / ex_;
    FloatType const hz = lz_ / ez_;
    // Coordinate of the global node line g: element g / order, local node g % order (the last
    // line, g = n_elements * order, is local node 0 of a virtual element: the domain end)
    auto coord1d = [&](int g, FloatType h, FloatType origin) {
      return origin + (g / order_) * h + (GLLPoints::get(order_, g % order_) + 1) * h / 2;
    };
    for (int gk = 0; gk < nz; ++gk) {
      FloatType const zk = coord1d(gk, hz, oz_);
      for (int gi = 0; gi < nx; ++gi) {
        x(gi + gk * nx) = coord1d(gi, hx, ox_);
        z(gi + gk * nx) = zk;
      }
    }
    nodes_coords_x_.modify_host();
    nodes_coords_y_.modify_host();  // zero-initialized
    nodes_coords_z_.modify_host();
  }

  /// Homogeneous acoustic material: vp = 1500 m/s, rho = 1.
  void initModels() {
    int const n = isModelOnNodes_ ? numberOfNodes() : ex_ * ez_;
    auto& vp = isModelOnNodes_ ? model_vp_node_ : model_vp_element_;
    auto& rho = isModelOnNodes_ ? model_rho_node_ : model_rho_element_;
    vp = allocateVector<dualVectorReal>(n, "model vp");
    rho = allocateVector<dualVectorReal>(n, "model rho");
    Kokkos::deep_copy(vp.view_host(), static_cast<FloatType>(1500));
    Kokkos::deep_copy(rho.view_host(), static_cast<FloatType>(1));
    vp.modify_host();
    rho.modify_host();
  }
};

}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_UNSTRUCT_2D_BUILDER_H_
