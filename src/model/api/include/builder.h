#pragma once

#include <model.h>

namespace model {
    template<typename FloatType, typename ScalarType>
    class ModelBuilderBase {
     public:
      ModelBuilderBase() = default;
      ~ModelBuilderBase() = default;

      std::unique_ptr<model::ModelApi<FloatType, ScalarType>> getModel() const ;
    };
}
