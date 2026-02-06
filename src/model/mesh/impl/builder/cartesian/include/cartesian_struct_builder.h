#pragma once

#include <builder.h>
#include <model_struct.h>

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
                         FloatType oz = 0.0,
                         // AJOUTER : paramètres globaux (optionnels)
                         FloatType global_lx = -1.0, FloatType global_ly = -1.0,
                         FloatType global_lz = -1.0, FloatType global_ox = 0.0,
                         FloatType global_oy = 0.0, FloatType global_oz = 0.0)
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
        global_lx_(global_lx < 0 ? lx : global_lx),  // Si pas fourni, = local
        global_ly_(global_ly < 0 ? ly : global_ly),
        global_lz_(global_lz < 0 ? lz : global_lz),
        global_ox_(global_ox),
        global_oy_(global_oy),
        global_oz_(global_oz)
  {
  }

  ~CartesianStructBuilder() = default;

  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const override
  {
    model::ModelStructData<FloatType, ScalarType> data;
    data.ex_ = ex_;
    data.ey_ = ey_;
    data.ez_ = ez_;

    data.dx_ = lx_;
    data.dy_ = ly_;
    data.dz_ = lz_;

    // Local origin
    data.ox_ = ox_;
    data.oy_ = oy_;
    data.oz_ = oz_;

    // AJOUTER : Global bounds
    data.ox_global_ = global_ox_;
    data.oy_global_ = global_oy_;
    data.oz_global_ = global_oz_;
    data.lx_global_ = global_lx_;
    data.ly_global_ = global_ly_;
    data.lz_global_ = global_lz_;

    data.isModelOnNodes_ = isModelOnNodes_;
    data.isElastic_ = isElastic_;

    auto model =
        std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(
            data);

    model->buildFaceConnectivity();
    model->initializeBoundaryFlags(true);  // true = free surface on top

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