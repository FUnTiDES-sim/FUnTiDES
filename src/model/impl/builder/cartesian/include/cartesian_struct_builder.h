#pragma once

#include <builder.h>
#include <model_struct.h>

namespace model {
    template<typename FloatType, typename ScalarType, int Order>
    class CartesianStructBuilder: public ModelBuilderBase<FloatType, ScalarType> {
     public:
       CartesianStructBuilder() = default;
       ~CartesianStructBuilder() = default;

       model::ModelStruct<FloatType, ScalarType, Order> getModel(ScalarType e, FloatType h) const {
          model::ModelStructData<FloatType, ScalarType> data;
          data.ex_ = e;
          data.ey_ = e;
          data.ez_ = e;

          data.dx_ = h;
          data.dy_ = h;
          data.dz_ = h;

          return model::ModelStruct<FloatType, ScalarType, Order>(data);
       }
    };
}
