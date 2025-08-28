#pragma once

#include <model_api.h>

namespace model_builder {
    template<typename FloatType, typename ScalarType>
    class ModelBuilderBase {
     public:
      ModelBuilderBase() = default;
      ~ModelBuilderBase() = default;

      model::ModelApi<FloatType, ScalarType> getModel() const;
    };
}
