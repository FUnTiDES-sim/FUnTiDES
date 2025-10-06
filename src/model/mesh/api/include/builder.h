#pragma once

#include <model.h>

#include <memory>

namespace model
{
namespace mesh
{
template <typename FloatType, typename ScalarType>
class ModelBuilderBase
{
 public:
  ModelBuilderBase() = default;
  ~ModelBuilderBase() = default;

  static constexpr int MAX_ORDER = 5;

  virtual std::shared_ptr<model::mesh::ModelApi<FloatType, ScalarType>> getModel()
      const = 0;
};
}  // namespace mesh
}  // namespace model
