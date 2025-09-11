#pragma once

#include <builder.h>
#include <model_struct.h>

namespace model {
    template<typename FloatType, typename ScalarType, int Order>
    class CartesianStructBuilder: public ModelBuilderBase<FloatType, ScalarType> {
     public:
       CartesianStructBuilder() = default;
       ~CartesianStructBuilder() = default;

       model::ModelStruct<FloatType, ScalarType, Order> getModel(ScalarType ex, FloatType hx,
                                                                 ScalarType ey, FloatType hy,
                                                                 ScalarType ez, FloatType hz) const {
          model::ModelStructData<FloatType, ScalarType> data;
          data.ex_ = ex;
          data.ey_ = ey;
          data.ez_ = ez;

          data.dx_ = hx;
          data.dy_ = hy;
          data.dz_ = hz;

          return model::ModelStruct<FloatType, ScalarType, Order>(data);
       }
    };
}
