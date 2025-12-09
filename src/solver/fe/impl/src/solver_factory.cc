#include "solver_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver.h"

namespace SolverFactory
{

namespace
{

/**
 * @brief Dispatches to the correct template instantiation based on runtime
 * order.
 */
template <typename FUNC>
std::unique_ptr<SEMSolverBase> orderDispatch(int const order, FUNC&& func)
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
 * @brief Creates solver for structured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<SEMSolverBase> makeSolverStruct(bool isModelOnNodes,
                                                physicType physic)
{
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, true, kAcoustic>>();
    }
    else
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, false, kAcoustic>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, true, kElastic>>();
    }
    else
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, false, kElastic>>();
    }
  }
}

/**
 * @brief Creates solver for unstructured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<SEMSolverBase> makeSolverUnstruct(bool isModelOnNodes,
                                                  physicType physic)
{
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == kAcoustic)
  {
    if (isModelOnNodes)
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, true, kAcoustic>>();
    }
    else
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, false, kAcoustic>>();
    }
  }
  else  // kElastic
  {
    if (isModelOnNodes)
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, true, kElastic>>();
    }
    else
    {
      return std::make_unique<solver::fe::SEMsolver<ORDER, SelectedIntegral,
                                                    MeshT, false, kElastic>>();
    }
  }
}

/**
 * @brief Creates a SEM solver with the specified integral implementation.
 */
template <auto ImplTag>
std::unique_ptr<SEMSolverBase> makeSemSolver(int order, meshType mesh,
                                             modelLocationType modelLocation,
                                             physicType physic)
{
  bool const isModelOnNodes = (modelLocation == kOnNodes);

  switch (mesh)
  {
    case kStruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    case kUnstruct:
      return orderDispatch(order, [&](auto orderIC) {
        constexpr int ORDER = decltype(orderIC)::value;
        return makeSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
      });

    default:
      throw std::runtime_error("Unknown mesh type");
  }
}

}  // anonymous namespace

std::unique_ptr<SEMSolverBase> createSolver(
    methodType const methodType, implemType const implemType,
    meshType const mesh, modelLocationType const modelLocation,
    physicType const physicType, int const order)
{
  if (methodType == kSem)
  {
    switch (implemType)
    {
      case kMakutu:
        return makeSemSolver<IntegralType::MAKUTU>(order, mesh, modelLocation,
                                                   physicType);
      case kShiva:
#ifdef ENABLE_Shiva
        return makeSemSolver<IntegralType::SHIVA>(order, mesh, modelLocation,
                                                  physicType);
#else
        throw std::runtime_error(
            "Shiva is not compiled. Cannot be used, please compile with "
            "ENABLE_Shiva=ON");
#endif  // ENABLE_Shiva

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

}  // namespace SolverFactory
