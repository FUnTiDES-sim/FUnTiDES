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

  /**
   * @brief Get the model instance.
   * @param free_surface_on_top Indicates if the free surface is on top, else we use damping on the top boundary.
   * @return A shared pointer to the model instance.
   */
  virtual std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel(
      bool free_surface_on_top) const = 0;
};
}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_
