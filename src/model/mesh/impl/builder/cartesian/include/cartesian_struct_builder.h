#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_

#include <builder.h>
#include <model_struct.h>

#include <algorithm>

#include "cartesian_model_file_reader.h"
#include "cartesian_struct_boundary_classifier.h"

namespace model {
/**
 * @brief Builds a structured Cartesian hexahedral model for one subdomain.
 *
 * The builder describes the local subdomain (element counts, sizes, origin) and the
 * global domain it belongs to. The global extents are used to classify boundary nodes,
 * so that only faces on the global boundary are tagged as boundaries.
 *
 * @tparam FloatType Floating-point type of coordinates and material values.
 * @tparam ScalarType Integer type of element counts.
 * @tparam Order Polynomial order of the spectral elements.
 */
template <typename FloatType, typename ScalarType, int Order>
class CartesianStructBuilder : public ModelBuilderBase<FloatType, ScalarType> {
 public:
  /**
   * @brief Stores the description of the local subdomain and of the global domain.
   *
   * @param[in] ex,ey,ez Number of local elements along x, y, z.
   * @param[in] lx,ly,lz Local subdomain size along x, y, z (meters).
   * @param[in] isModelOnNodes True if the material model is given per node, false if per element.
   * @param[in] isElastic True for an elastic model, false for an acoustic one.
   * @param[in] ox,oy,oz Coordinates of the local subdomain origin.
   * @param[in] global_lx,global_ly,global_lz Global domain size; a negative value means
   *            "same as the local size".
   * @param[in] global_ox,global_oy,global_oz Coordinates of the global domain origin.
   * @param[in] isAcoustoElastic If true, getModel() builds a two-layer fluid/solid model.
   * @param[in] acoustoElasticBoundaryZ Z coordinate of the fluid/solid interface; the fluid
   *            is the part with z at or above it. Only used when isAcoustoElastic is true.
   * @param[in] DgSemBoundaryZ @todo VERIFY: meaning and unit; stored but read by no method of this class.
   * @param[in] model_file Path of a per-element model file; empty means no file.
   *            Requires isModelOnNodes to be false.
   */
  CartesianStructBuilder(ScalarType ex, FloatType lx, ScalarType ey, FloatType ly, ScalarType ez, FloatType lz,
                         bool isModelOnNodes, bool isElastic, FloatType ox = 0.0, FloatType oy = 0.0,
                         FloatType oz = 0.0, FloatType global_lx = -1.0, FloatType global_ly = -1.0,
                         FloatType global_lz = -1.0, FloatType global_ox = 0.0, FloatType global_oy = 0.0,
                         FloatType global_oz = 0.0, bool isAcoustoElastic = false,
                         FloatType acoustoElasticBoundaryZ = static_cast<FloatType>(0),
                         FloatType DgSemBoundaryZ = static_cast<FloatType>(0), std::string model_file = "")
      : ex_(ex),
        ey_(ey),
        ez_(ez),
        lx_(lx),
        ly_(ly),
        lz_(lz),
        isModelOnNodes_(isModelOnNodes),
        isElastic_(isElastic),
        ox_(ox),
        oy_(oy),
        oz_(oz),
        global_lx_(global_lx < 0 ? lx : global_lx),
        global_ly_(global_ly < 0 ? ly : global_ly),
        global_lz_(global_lz < 0 ? lz : global_lz),
        global_ox_(global_ox),
        global_oy_(global_oy),
        global_oz_(global_oz),
        isAcoustoElastic_(isAcoustoElastic),
        acoustoElasticBoundaryZ_(acoustoElasticBoundaryZ),
        DgSemBoundaryZ_(DgSemBoundaryZ),
        model_file_(std::move(model_file)) {}

  ~CartesianStructBuilder() = default;

  /**
   * @brief Builds the structured model of the local subdomain, with its face connectivity.
   *
   * Boundary flags are classified against the global domain. If isAcoustoElastic is set,
   * per-node or per-element vp, vs and rho are filled with a fluid layer above the
   * interface and a solid below. If a model file was given, per-element vp, rho and vs
   * are read from it and replace the values of the corresponding properties.
   *
   * @param[in] free_surface_on_top Forwarded to the boundary classifier for the top (z max) face.
   * @return The built model.
   * @throws std::runtime_error If a model file is given with isModelOnNodes set, or if the
   *         number of vp values in the file differs from the number of elements.
   */
  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel(bool free_surface_on_top) const override {
    model::ModelStructData<FloatType, ScalarType> data;
    data.ex_ = ex_;
    data.ey_ = ey_;
    data.ez_ = ez_;
    data.dx_ = lx_;
    data.dy_ = ly_;
    data.dz_ = lz_;
    data.ox_ = ox_;
    data.oy_ = oy_;
    data.oz_ = oz_;
    data.isModelOnNodes_ = isModelOnNodes_;
    data.isElastic_ = isElastic_;

    const int nx = static_cast<int>(ex_) * Order + 1;
    const int ny = static_cast<int>(ey_) * Order + 1;
    const int nz = static_cast<int>(ez_) * Order + 1;
    const int n_node = nx * ny * nz;
    const FloatType tol = std::min({lx_ / ex_, ly_ / ey_, lz_ / ez_}) * static_cast<FloatType>(1e-4);

    data.boundaries_t_ = CartesianStructBoundaryClassifier<FloatType, ScalarType>(
                             global_ox_, global_ox_ + global_lx_, global_oy_, global_oy_ + global_ly_, global_oz_,
                             global_oz_ + global_lz_, tol, free_surface_on_top)
                             .classify(n_node, nx, ny, nz, ox_, oy_, oz_, lx_, ly_, lz_);

    // Two-layer model: fluid (z >= boundary) over solid. vs = 0 in the fluid, which is
    // what TagElements uses to classify an element as acoustic.
    if (isAcoustoElastic_) {
      auto temp_model = model::ModelStruct<FloatType, ScalarType, Order>(data);

      if (isModelOnNodes_) {
        data.model_vp_node_ = allocateVector<vectorReal>(n_node, "model_vp_node");
        data.model_vs_node_ = allocateVector<vectorReal>(n_node, "model_vs_node");
        data.model_rho_node_ = allocateVector<vectorReal>(n_node, "model_rho_node");

        for (int n = 0; n < n_node; ++n) {
          bool const is_fluid = (temp_model.nodeCoord(n, 2) >= acoustoElasticBoundaryZ_);
          data.model_vp_node_[n] = is_fluid ? static_cast<FloatType>(1500) : static_cast<FloatType>(3400);
          data.model_vs_node_[n] = is_fluid ? static_cast<FloatType>(0) : static_cast<FloatType>(1963);
          data.model_rho_node_[n] = is_fluid ? static_cast<FloatType>(1020) : static_cast<FloatType>(2500);
        }
      } else {
        int const n_elem = ex_ * ey_ * ez_;
        FloatType const hz = lz_ / ez_;

        data.model_vp_element_ = allocateVector<vectorReal>(n_elem, "model_vp_elem");
        data.model_vs_element_ = allocateVector<vectorReal>(n_elem, "model_vs_elem");
        data.model_rho_element_ = allocateVector<vectorReal>(n_elem, "model_rho_elem");

        // An element is fluid if its centroid is at or above the interface.
        for (int k = 0; k < ez_; ++k) {
          FloatType const centroid_z = oz_ + (k + static_cast<FloatType>(0.5)) * hz;
          bool const is_fluid = (centroid_z >= acoustoElasticBoundaryZ_);
          for (int j = 0; j < ey_; ++j) {
            for (int i = 0; i < ex_; ++i) {
              int const e = i + j * ex_ + k * ex_ * ey_;
              data.model_vp_element_[e] = is_fluid ? static_cast<FloatType>(1500) : static_cast<FloatType>(3400);
              data.model_vs_element_[e] = is_fluid ? static_cast<FloatType>(0) : static_cast<FloatType>(1963);
              data.model_rho_element_[e] = is_fluid ? static_cast<FloatType>(1020) : static_cast<FloatType>(2500);
            }
          }
        }
      }
    }

    if (!model_file_.empty()) {
      if (isModelOnNodes_)
        throw std::runtime_error(
            "[CartesianStructBuilder] model_file only supported with isModelOnNodes=false for now.");
      int const n_elem = ex_ * ey_ * ez_;
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
      fill_view(data.model_vp_element_, "Vp", "model_vp_element");
      fill_view(data.model_rho_element_, "Rho", "model_rho_element");
      fill_view(data.model_vs_element_, "Vs", "model_vs_element");
      if (data.model_vp_element_.extent(0) != static_cast<size_t>(n_elem))
        throw std::runtime_error("[CartesianStructBuilder] model_file has " +
                                 std::to_string(data.model_vp_element_.extent(0)) + " elements but mesh has " +
                                 std::to_string(n_elem));
      std::cout << "[Model] vp[0]=" << data.model_vp_element_[0] << " vp[mid]=" << data.model_vp_element_[n_elem / 2]
                << " vp[last]=" << data.model_vp_element_[n_elem - 1] << std::endl;
    }

    // The model uses the local origin and sizes; boundary flags were classified
    // above against the global domain.
    auto model = std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(data);

    model->buildFaceConnectivity();

    return model;
  }

 private:
  FloatType ox_, oy_, oz_;                       ///< Local subdomain origin.
  FloatType global_ox_, global_oy_, global_oz_;  ///< Global domain origin.
  ScalarType ex_, ey_, ez_;                      ///< Number of local elements along each axis.
  FloatType lx_, ly_, lz_;                       ///< Local subdomain size (meters).
  FloatType global_lx_, global_ly_, global_lz_;  ///< Global domain size (meters).
  bool isModelOnNodes_;                          ///< Material model given per node (true) or per element (false).
  bool isElastic_;                               ///< Elastic (true) or acoustic (false) model.
  bool isAcoustoElastic_{false};                 ///< Build a fluid/solid two-layer model.
  FloatType acoustoElasticBoundaryZ_{static_cast<FloatType>(0)};  ///< Z of the fluid/solid interface.
  FloatType DgSemBoundaryZ_{static_cast<FloatType>(0)};           ///< Unused, see constructor.
  std::string model_file_;                       ///< Per-element model file path, empty if none.
};
}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_
