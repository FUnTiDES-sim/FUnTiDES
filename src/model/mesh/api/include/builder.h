#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_

#pragma once

#include <model.h>

#include <memory>

namespace model
{
template <typename FloatType, typename ScalarType>
class ModelBuilderBase
{
 public:
  ModelBuilderBase() = default;
  ~ModelBuilderBase() = default;

  static constexpr int MAX_ORDER = 5;

  virtual std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const = 0;
};
}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_
