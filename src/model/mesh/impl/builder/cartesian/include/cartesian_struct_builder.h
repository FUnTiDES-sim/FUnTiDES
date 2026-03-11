#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_

#include <builder.h>
#include <model_struct.h>

#include <algorithm>

#include "cartesian_struct_boundary_classifier.h"

namespace model
{
template <typename FloatType, typename ScalarType, int Order>
class CartesianStructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  CartesianStructBuilder(ScalarType ex, FloatType lx, ScalarType ey,
                         FloatType ly, ScalarType ez, FloatType lz,
                         bool isModelOnNodes, bool isElastic,
                         FloatType ox = 0.0, FloatType oy = 0.0,
                         FloatType oz = 0.0, FloatType global_lx = -1.0,
                         FloatType global_ly = -1.0, FloatType global_lz = -1.0,
                         FloatType global_ox = 0.0, FloatType global_oy = 0.0,
                         FloatType global_oz = 0.0)
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
        global_oz_(global_oz)
  {
  }

  ~CartesianStructBuilder() = default;

  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel(
      bool free_surface_on_top) const override
  {
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
    const FloatType tol = std::min({lx_ / ex_, ly_ / ey_, lz_ / ez_}) *
                          static_cast<FloatType>(1e-4);

    data.boundaries_t_ =
        CartesianStructBoundaryClassifier<FloatType, ScalarType>(
            global_ox_, global_ox_ + global_lx_, global_oy_,
            global_oy_ + global_ly_, global_oz_, global_oz_ + global_lz_, tol,
            free_surface_on_top)
            .classify(n_node, nx, ny, nz, ox_, oy_, oz_, lx_, ly_, lz_);

    // -------------------------------------------------------------------------
    // Construct model with local coordinates and dimensions, but use global
    // boundaries for boundary classification.
    // -------------------------------------------------------------------------
    auto model =
        std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(
            data);

    model->buildFaceConnectivity();

    return model;
  }

 private:
  FloatType ox_, oy_, oz_;  // Local origin coordinate in 3D
  FloatType global_ox_, global_oy_, global_oz_;  // Global origin
  ScalarType ex_, ey_, ez_;  // Number of elements for each axis (local)
  FloatType lx_, ly_, lz_;   // Domain size (local)
  FloatType global_lx_, global_ly_, global_lz_;  // Domain size (global)
  bool isModelOnNodes_;
  bool isElastic_;
};
}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_
