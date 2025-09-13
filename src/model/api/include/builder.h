#pragma once

#include <model.h>
#include <memory>

namespace model {
    template<typename FloatType, typename ScalarType>
    class ModelBuilderBase {
     public:
      ModelBuilderBase() = default;
      ~ModelBuilderBase() = default;

      static constexpr int MAX_ORDER = 5;

      std::unique_ptr<model::ModelApi<FloatType, ScalarType>> getModel() const;
    };
}
