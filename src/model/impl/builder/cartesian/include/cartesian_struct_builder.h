#pragma once

#include <builder.h>
#include <model_struct.h>

namespace model
{
template <typename FloatType, typename ScalarType, int Order>
class CartesianStructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  CartesianStructBuilder(ScalarType ex, FloatType hx, ScalarType ey,
                         FloatType hy, ScalarType ez, FloatType hz,
                         bool isModelOnNodes)
      : ex_(ex),
        ey_(ey),
        ez_(ez),
        hx_(hx),
        hy_(hy),
        hz_(hz),
        isModelOnNodes_(isModelOnNodes)
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

    data.dx_ = hx_;
    data.dy_ = hy_;
    data.dz_ = hz_;
    data.isModelOnNodes_ = isModelOnNodes_;

    return std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(
        data);
  }

 private:
  ScalarType ex_, ey_, ez_;
  FloatType hx_, hy_, hz_;
  bool isModelOnNodes_;
};
}  // namespace model
