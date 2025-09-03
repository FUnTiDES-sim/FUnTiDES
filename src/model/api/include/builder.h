#pragma once

#include <model.h>

namespace model {
    template<typename FloatType, typename ScalarType>
    class ModelBuilderBase {
     public:
      ModelBuilderBase() = default;
      ~ModelBuilderBase() = default;

      model::ModelApi<FloatType, ScalarType> getModel() const;
    };
}
