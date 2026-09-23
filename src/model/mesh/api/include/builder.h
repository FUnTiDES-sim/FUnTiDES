#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_

#pragma once

#include <model.h>

#include <memory>

namespace model {
/**
 * @brief Abstract factory of ModelApi instances.
 *
 * A builder is configured with its mesh and material parameters at construction; getModel()
 * then creates the model. Host only.
 *
 * @tparam FloatType Floating-point type of the built model.
 * @tparam ScalarType Integer type of the built model.
 */
template <typename FloatType, typename ScalarType>
class ModelBuilderBase {
 public:
  ModelBuilderBase() = default;
  ~ModelBuilderBase() = default;

  static constexpr int MAX_ORDER = 9;  ///< Highest supported polynomial order.

  /**
   * @brief Build the model.
   * @param[in] free_surface_on_top true to flag the nodes of the z-max global face as Surface,
   * false to flag them as Damping like the other global boundary faces.
   * @return The new model.
   */
  virtual std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel(bool free_surface_on_top) const = 0;
};
}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_BUILDER_H_
