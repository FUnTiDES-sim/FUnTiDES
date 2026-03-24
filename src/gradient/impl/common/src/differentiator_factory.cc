#include "differentiator_factory.h"

#include <differentiator_acoustic.h>
#include <differentiator_elastic.h>
#include <model_struct.h>
#include <model_unstruct.h>

#include <fe/Integrals.hpp>

namespace gradient
{

using physicType = utils::enums::physicType;
using meshType = utils::enums::meshType;
using modelLocationType = utils::enums::modelLocationType;
using implemType = utils::enums::implemType;
using methodType = utils::enums::methodType;

/**
 * @brief Dispatches to the correct template instantiation based on runtime
 * order.
 */
template <typename FUNC>
std::unique_ptr<Differentiator> orderDispatch(int const order, FUNC&& func)
{
  switch (order)
  {
    case 1:
      return func(std::integral_constant<int, 1>{});
    case 2:
      return func(std::integral_constant<int, 2>{});
    case 3:
      return func(std::integral_constant<int, 3>{});
    // case 4:
    //   return func(std::integral_constant<int, 4>{});
    default:
      throw std::runtime_error("Unsupported polynomial order: " +
                               std::to_string(order));
  }
}

/**
 * @brief Creates differentiator for structured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Differentiator> makeDifferentiatorStruct(bool isModelOnNodes,
                                                         physicType physic)
{
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == physicType::kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
    }
    else
    {
      return std::make_unique<
          DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, true>>();
    }
    else
    {
      return std::make_unique<
          DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
}

/**
 * @brief Creates differentiator for unstructured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Differentiator> makeDifferentiatorUnstruct(bool isModelOnNodes,
                                                           physicType physic)
{
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == physicType::kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
    }
    else
    {
      return std::make_unique<
          DifferentiatorAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, true>>();
    }
    else
    {
      return std::make_unique<
          DifferentiatorElastic<ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
}

/**
 * @brief Creates a SEM solver with the specified integral implementation.
 */
template <auto ImplTag>
std::unique_ptr<Differentiator> makeDifferentiatorSem(
    int order, meshType mesh, modelLocationType modelLocation,
    physicType physic)
{
  bool const isModelOnNodes = (modelLocation == modelLocationType::kOnNodes);

  switch (mesh)
  {
    case meshType::kStruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeDifferentiatorStruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    case meshType::kUnstruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeDifferentiatorUnstruct<ImplTag, ORDER>(isModelOnNodes,
                                                          physic);
      });

    default:
      throw std::runtime_error("Unknown mesh type");
  }
}

std::unique_ptr<Differentiator> createDifferentiator(
    methodType const methodType, implemType const implemType,
    meshType const mesh, modelLocationType const modelLocation,
    physicType const physicType, int const order)
{
  if (methodType == methodType::kSem)
  {
    switch (implemType)
    {
      case implemType::kMakutu:
        return makeDifferentiatorSem<IntegralType::MAKUTU>(
            order, mesh, modelLocation, physicType);
      default:
        throw std::runtime_error("Unknown implementation type: " +
                                 to_string(implemType));
    }
  }

  // Add DG or other methods as needed
  throw std::runtime_error(
      "Unsupported solver configuration: methodType=" + to_string(methodType) +
      ", implemType=" + to_string(implemType) +
      ", physicType=" + to_string(physicType));
}

}  // namespace gradient
