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
                         FloatType oz = 0.0)
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
        oz_(oz)
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

    // Pass origin to data
    data.ox_ = ox_;
    data.oy_ = oy_;
    data.oz_ = oz_;

    data.isModelOnNodes_ = isModelOnNodes_;
    data.isElastic_ = isElastic_;

    auto model =
        std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(
            data);

    if (isElastic_ && !isModelOnNodes_)
    {
      model->initElasticityTensors();
    }

    return model;
  }

 private:
  FloatType ox_, oy_, oz_;   //< origin coordinate in 3D
  ScalarType ex_, ey_, ez_;  //< number of elements for each axis
  FloatType lx_, ly_, lz_;   //< domain size
  bool isModelOnNodes_;
  bool isElastic_;
};
}  // namespace model
