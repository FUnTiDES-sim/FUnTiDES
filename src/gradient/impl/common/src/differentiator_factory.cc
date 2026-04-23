#include "differentiator_factory.h"

#include <Integrals.h>
#include <differentiator_acoustic.h>
#include <differentiator_elastic.h>
#include <model_struct.h>
#include <model_unstruct.h>

namespace gradient {

using physicType = utils::enums::physicType;
using meshType = utils::enums::meshType;
using modelLocationType = utils::enums::modelLocationType;
using implemType = utils::enums::implemType;

/**
 * @brief Template récursif C++17 qui remplace le switch codé en dur.
 * Il s'arrête exactement à CurrentOrder (défini par CMake) et descend
 * jusqu'à 1.
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
 * @brief Creates differentiator for structured mesh.
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
 * @brief Creates differentiator for unstructured mesh.
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
 * @brief Creates a SEM solver with the specified integral implementation.
 */
template <auto ImplTag>
std::unique_ptr<Differentiator> makeDifferentiatorSem(int order, meshType mesh, modelLocationType modelLocation,
                                                      physicType physic) {
  bool const isModelOnNodes = (modelLocation == modelLocationType::kOnNodes);

// Définir une macro de fallback au cas où CMake ne passe pas les variables
#ifndef MAX_DIFFERENTIATOR_ACOUSTIC_ORDER
#define MAX_DIFFERENTIATOR_ACOUSTIC_ORDER 3
#endif
#ifndef MAX_DIFFERENTIATOR_ELASTIC_ORDER
#define MAX_DIFFERENTIATOR_ELASTIC_ORDER 3
#endif

  // On dispatch selon la physique D'ABORD, pour injecter la bonne limite
  // d'ordre max.
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
