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
 * @brief Creates a differentiator on a structured mesh.
 *
 * @tparam ImplTag Integral back-end tag passed to IntegralTypeSelector.
 * @tparam ORDER Polynomial order.
 * @param[in] isModelOnNodes True if the model is defined on nodes, false if per element.
 * @param[in] physic Physics of the differentiator; any value other than kAcoustic
 *                   selects the elastic one.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Differentiator> makeDifferentiatorStruct(bool isModelOnNodes, physicType physic) {
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == physicType::kAcoustic) {
    if (isModelOnNodes) {
      return std::make_unique<DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
    } else {
      return std::make_unique<DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  } else  // kElastic
  {
    if (isModelOnNodes) {
      return std::make_unique<DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, true>>();
    } else {
      return std::make_unique<DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
}

/**
 * @brief Creates a differentiator on an unstructured mesh.
 *
 * @tparam ImplTag Integral back-end tag passed to IntegralTypeSelector.
 * @tparam ORDER Polynomial order.
 * @param[in] isModelOnNodes True if the model is defined on nodes, false if per element.
 * @param[in] physic Physics of the differentiator; any value other than kAcoustic
 *                   selects the elastic one.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Differentiator> makeDifferentiatorUnstruct(bool isModelOnNodes, physicType physic) {
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == physicType::kAcoustic) {
    if (isModelOnNodes) {
      return std::make_unique<DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
    } else {
      return std::make_unique<DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  } else  // kElastic
  {
    if (isModelOnNodes) {
      return std::make_unique<DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, true>>();
    } else {
      return std::make_unique<DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
}

/**
 * @brief Creates a differentiator for one integral back-end, from runtime options.
 *
 * @tparam ImplTag Integral back-end tag.
 * @param[in] order Polynomial order, limited by the per-physics maximum order.
 * @param[in] mesh Structured or unstructured mesh.
 * @param[in] modelLocation Whether the model is defined on nodes or per element.
 * @param[in] physic Acoustic or elastic.
 * @throws std::runtime_error If the order is unsupported or the physics is unknown.
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
  if (physic == physicType::kAcoustic) {
    return orderDispatch<MAX_DIFFERENTIATOR_ACOUSTIC_ORDER>(order, [&](auto orderIC) {
      constexpr int ORDER = decltype(orderIC)::value;
      return (mesh == meshType::kStruct) ? makeDifferentiatorStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                         : makeDifferentiatorUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
    });
  } else if (physic == physicType::kElastic) {
    return orderDispatch<MAX_DIFFERENTIATOR_ELASTIC_ORDER>(order, [&](auto orderIC) {
      constexpr int ORDER = decltype(orderIC)::value;
      return (mesh == meshType::kStruct) ? makeDifferentiatorStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                         : makeDifferentiatorUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
    });
  }

  throw std::runtime_error("Unknown physics type");
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
