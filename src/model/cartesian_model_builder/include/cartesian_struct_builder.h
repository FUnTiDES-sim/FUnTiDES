#pragma once

#include "model_builder_api.h"
#include <model_struct.h>

namespace model_builder {
    template<typename FloatType, typename ScalarType>
    class CartesianStructBuilder: public ModelBuilderBase<FloatType, ScalarType> {
     public:
       CartesianStructBuilder() = default;
       ~CartesianStructBuilder() = default;

       model::ModelStruct<FloatType, ScalarType> getModel(ScalarType e, FloatType h, int order) const {
          model::ModelStructData<FloatType, ScalarType> data;
          data.ex_ = e;
          data.ey_ = e;
          data.ez_ = e;

          data.dx_ = h;
          data.dy_ = h;
          data.dz_ = h;

          data.order_ = order;

          return model::ModelStruct<FloatType, ScalarType>(data);
       }
    };
}
