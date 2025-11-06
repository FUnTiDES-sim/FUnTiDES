#pragma once

#include <builder.h>
#include <model_struct.h>

namespace model
{
template <typename FloatType, typename ScalarType, int Order>
class CartesianStructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  CartesianStructBuilder(ScalarType ex, FloatType lx, ScalarType ey, FloatType ly, ScalarType ez, FloatType lz,
                         bool isModelOnNodes, bool isElastic)
      : ex_(ex), ey_(ey), ez_(ez), lx_(lx), ly_(ly), lz_(lz), isModelOnNodes_(isModelOnNodes), isElastic_(isElastic)
  {
  }

  ~CartesianStructBuilder() = default;

  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel() const override
  {
    model::ModelStructData<FloatType, ScalarType> data;
    data.ex_ = ex_;
    data.ey_ = ey_;
    data.ez_ = ez_;

    data.dx_ = lx_;
    data.dy_ = ly_;
    data.dz_ = lz_;
    data.isModelOnNodes_ = isModelOnNodes_;
    data.isElastic_ = isElastic_;

    return std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(data);
  }

 private:
  ScalarType ex_, ey_, ez_;
  FloatType lx_, ly_, lz_;
  bool isModelOnNodes_;
  bool isElastic_;
};
}  // namespace model
