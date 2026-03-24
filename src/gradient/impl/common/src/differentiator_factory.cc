#include "differentiator_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "differentiator.h"

namespace gradient
{


namespace feenum = solver::fe::enums;

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
                                                 feenum::physicType physic)
{
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kAcoustic>>();
    }
    else
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kAcoustic>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kElastic>>();
    }
    else
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kElastic>>();
    }
  }
}

/**
 * @brief Creates differentiator for unstructured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Differentiator> makeDifferentiatorUnstruct(bool isModelOnNodes,
                                           feenum::physicType physic)
{
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kAcoustic>>();
    }
    else
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kAcoustic>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kElastic>>();
    }
    else
    {
      return std::make_unique<
          Differentiator<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kElastic>>();
    }
  }
}

/**
 * @brief Creates a SEM solver with the specified integral implementation.
 */
template <auto ImplTag>
std::unique_ptr<Differentiator> makeDifferentiatorSem(int order, feenum::meshType mesh,
                                      feenum::modelLocationType modelLocation,
                                      feenum::physicType physic)
{
  bool const isModelOnNodes =
      (modelLocation == feenum::modelLocationType::kOnNodes);

  switch (mesh)
  {
    case feenum::meshType::kStruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeDifferentiatorStruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    case feenum::meshType::kUnstruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeDifferentiatorUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    default:
      throw std::runtime_error("Unknown mesh type");
  }
}

std::unique_ptr<Differentiator> createDifferentiator(
    feenum::methodType const methodType, feenum::implemType const implemType,
    feenum::meshType const mesh, feenum::modelLocationType const modelLocation,
    feenum::physicType const physicType, int const order)
{
  if (methodType == feenum::methodType::kSem)
  {
    switch (implemType)
    {
      case feenum::implemType::kMakutu:
        return makeDifferentiatorSem<IntegralType::MAKUTU>(order, mesh, modelLocation,
                                                   physicType);
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

}  // namespace gradient
