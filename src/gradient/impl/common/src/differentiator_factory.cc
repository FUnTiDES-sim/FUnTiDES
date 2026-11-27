#include "differentiator_factory.h"

#include "Integrals.h"
#include "differentiator_acoustic.h"
#include "differentiator_elastic.h"
#include "model_struct.h"
#include "model_unstruct.h"

namespace gradient {

using physicType = utils::enums::physicType;
using meshType = utils::enums::meshType;
using modelLocationType = utils::enums::modelLocationType;
using implemType = utils::enums::implemType;

/**
 * @brief Maps a runtime polynomial order to a compile-time constant.
 *
 * Tries CurrentOrder, then recurses down to 1, and calls @p func with the matching
 * std::integral_constant<int, order>.
 *
 * @tparam CurrentOrder Highest order that can be dispatched.
 * @param[in] order Runtime polynomial order.
 * @param[in] func Callable taking a std::integral_constant<int, order> and returning a
 *                 std::unique_ptr<Differentiator>.
 * @return The differentiator returned by @p func.
 * @throws std::runtime_error If @p order is not in [1, CurrentOrder].
 */
template <int CurrentOrder, typename FUNC>
std::unique_ptr<Differentiator> orderDispatch(int const order, FUNC&& func) {
  if (order == CurrentOrder) {
    return func(std::integral_constant<int, CurrentOrder>{});
  }

  if constexpr (CurrentOrder > 1) {
    return orderDispatch<CurrentOrder - 1>(order, std::forward<FUNC>(func));
  } else {
    throw std::runtime_error("Unsupported polynomial order: " + std::to_string(order));
  }
}

/**
 * @brief Creates a differentiator of physics PHYSIC for one mesh type.
 *
 * @tparam ImplTag Integral back-end tag passed to IntegralTypeSelector.
 * @tparam PHYSIC kAcoustic or kElastic.
 * @tparam ORDER Polynomial order.
 * @tparam MeshT Mesh type; structured meshes carry the polynomial order as a template argument.
 * @param[in] isModelOnNodes True if the model is defined on nodes, false if per element.
 */
template <auto ImplTag, physicType PHYSIC, int ORDER, typename MeshT>
std::unique_ptr<Differentiator> makeDifferentiatorForMesh(bool isModelOnNodes) {
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if constexpr (PHYSIC == physicType::kAcoustic) {
    if (isModelOnNodes) return std::make_unique<DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
    return std::make_unique<DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
  } else {
    if (isModelOnNodes) return std::make_unique<DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, true>>();
    return std::make_unique<DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, false>>();
  }
}

/**
 * @brief Creates a differentiator of physics PHYSIC and order ORDER, selecting the mesh type.
 */
template <auto ImplTag, physicType PHYSIC, int ORDER>
std::unique_ptr<Differentiator> makeDifferentiatorForOrder(bool isModelOnNodes, meshType mesh) {
  if (mesh == meshType::kStruct) {
    return makeDifferentiatorForMesh<ImplTag, PHYSIC, ORDER, model::ModelStruct<float, int, ORDER>>(isModelOnNodes);
  }
  return makeDifferentiatorForMesh<ImplTag, PHYSIC, ORDER, model::ModelUnstruct<float, int>>(isModelOnNodes);
}

/**
 * @brief Creates a differentiator for one integral back-end, from runtime options.
 *
 * Only the physics enabled by the COMPILE_ACOUSTIC and COMPILE_ELASTIC macros are available.
 *
 * @tparam ImplTag Integral back-end tag.
 * @param[in] order Polynomial order, limited by the per-physics maximum order.
 * @param[in] mesh Structured or unstructured mesh.
 * @param[in] modelLocation Whether the model is defined on nodes or per element.
 * @param[in] physic Acoustic or elastic.
 * @throws std::runtime_error If the order is unsupported, or the physics is unknown or not compiled.
 */
template <auto ImplTag>
std::unique_ptr<Differentiator> makeDifferentiatorSem(int order, meshType mesh, modelLocationType modelLocation,
                                                      physicType physic) {
  bool const isModelOnNodes = (modelLocation == modelLocationType::kOnNodes);

// Fallback maximum orders, used when the build system does not define them.
#ifndef MAX_DIFFERENTIATOR_ACOUSTIC_ORDER
#define MAX_DIFFERENTIATOR_ACOUSTIC_ORDER 3
#endif
#ifndef MAX_DIFFERENTIATOR_ELASTIC_ORDER
#define MAX_DIFFERENTIATOR_ELASTIC_ORDER 3
#endif

  // Dispatch on the physics first, so that each physics uses its own maximum order.
  switch (physic) {
    case physicType::kAcoustic:
#ifdef COMPILE_ACOUSTIC
      return orderDispatch<MAX_DIFFERENTIATOR_ACOUSTIC_ORDER>(order, [&](auto orderIC) {
        return makeDifferentiatorForOrder<ImplTag, physicType::kAcoustic, decltype(orderIC)::value>(isModelOnNodes,
                                                                                                    mesh);
      });
#else
      throw std::runtime_error("Differentiator: acoustic physics not compiled in this build (COMPILE_ACOUSTIC=OFF)");
#endif
    case physicType::kElastic:
#ifdef COMPILE_ELASTIC
      return orderDispatch<MAX_DIFFERENTIATOR_ELASTIC_ORDER>(order, [&](auto orderIC) {
        return makeDifferentiatorForOrder<ImplTag, physicType::kElastic, decltype(orderIC)::value>(isModelOnNodes,
                                                                                                   mesh);
      });
#else
      throw std::runtime_error("Differentiator: elastic physics not compiled in this build (COMPILE_ELASTIC=OFF)");
#endif
    default:
      throw std::runtime_error("Unknown physics type");
  }
}

/**
 * @brief Creates the differentiator matching the given runtime options.
 *
 * @param[in] implemType Integral implementation.
 * @param[in] mesh Structured or unstructured mesh.
 * @param[in] modelLocation Whether the model is defined on nodes or per element.
 * @param[in] physicType Acoustic or elastic.
 * @param[in] order Polynomial order.
 * @return Newly allocated differentiator, owned by the caller.
 * @throws std::runtime_error If the implementation, physics or order is unsupported.
 */
std::unique_ptr<Differentiator> createDifferentiator(implemType const implemType, meshType const mesh,
                                                     modelLocationType const modelLocation, physicType const physicType,
                                                     int const order) {
  switch (implemType) {
    case implemType::kMakutu:
      return makeDifferentiatorSem<IntegralType::MAKUTU>(order, mesh, modelLocation, physicType);
    default:
      throw std::runtime_error("Unknown implementation type: " + std::to_string(static_cast<int>(implemType)));
  }
}

}  // namespace gradient
