#pragma once

#include <model.h>
#include <memory>

namespace model {
    template<typename FloatType, typename ScalarType>
    class ModelBuilderBase {
     public:
      ModelBuilderBase() = default;
      ~ModelBuilderBase() = default;

      virtual std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel() const = 0;
    };
}
