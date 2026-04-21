#include "solver_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver.h"
#include "sem_solver_acoustoelastic.h"

namespace solver
{

namespace fe
{

namespace solver_factory
{

namespace feenum = utils::enums;

/**
 * @brief Dispatches to the correct template instantiation based on runtime
 * order.
 */
template <typename FUNC>
std::unique_ptr<Solver> orderDispatch(int const order, FUNC&& func)
{
  switch (order)
  {
    case 1:
      return func(std::integral_constant<int, 1>{});
    case 2:
      return func(std::integral_constant<int, 2>{});
    case 3:
      return func(std::integral_constant<int, 3>{});
    case 4:
      return func(std::integral_constant<int, 4>{});
    case 5:
      return func(std::integral_constant<int, 5>{});
    case 6:
      return func(std::integral_constant<int, 6>{});
    default:
      throw std::runtime_error("Unsupported polynomial order: " +
                               std::to_string(order));
  }
}

/**
 * @brief Creates solver for structured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeSolverStruct(bool isModelOnNodes,
                                         feenum::physicType physic)
{
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kAcoustic>>();
    }
    else
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kAcoustic>>();
    }
  }
  else if (physic == feenum::physicType::kAcoustoElastic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<
          ORDER, SelectedIntegral, MeshT, true>>();
    }
    else
    {
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<
          ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kElastic>>();
    }
    else
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kElastic>>();
    }
  }
}

/**
 * @brief Creates solver for unstructured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeSolverUnstruct(bool isModelOnNodes,
                                           feenum::physicType physic)
{
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kAcoustic>>();
    }
    else
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kAcoustic>>();
    }
  }
  else if (physic == feenum::physicType::kAcoustoElastic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<
          ORDER, SelectedIntegral, MeshT, true>>();
    }
    else
    {
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<
          ORDER, SelectedIntegral, MeshT, false>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true,
                                feenum::physicType::kElastic>>();
    }
    else
    {
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false,
                                feenum::physicType::kElastic>>();
    }
  }
}

/**
 * @brief Creates a SEM solver with the specified integral implementation.
 */
template <auto ImplTag>
std::unique_ptr<Solver> makeSemSolver(int order, feenum::meshType mesh,
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
        return makeSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    case feenum::meshType::kUnstruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    default:
      throw std::runtime_error("Unknown mesh type");
  }
}

std::unique_ptr<Solver> createSolver(
    feenum::methodType const methodType, feenum::implemType const implemType,
    feenum::meshType const mesh, feenum::modelLocationType const modelLocation,
    feenum::physicType const physicType, int const order)
{
  if (methodType == feenum::methodType::kSem)
  {
    switch (implemType)
    {
      case feenum::implemType::kMakutu:
        return makeSemSolver<IntegralType::MAKUTU>(order, mesh, modelLocation,
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
}  // namespace solver_factory
}  // namespace solver_factory
}  // namespace fe
}  // namespace solver
