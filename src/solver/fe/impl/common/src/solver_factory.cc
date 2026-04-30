#include "solver_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver.h"
#include "sem_solver_acoustoelastic.h"

namespace solver {
namespace fe {
namespace solver_factory {

namespace feenum = utils::enums;

template <int CurrentOrder, typename FUNC>
std::unique_ptr<Solver> orderDispatch(int const order, FUNC&& func) {
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
 * @brief Creates solver for structured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeSolverStruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  } else if (physic == feenum::physicType::kAcoustoElastic) {
    if (isModelOnNodes)
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<ORDER, SelectedIntegral, MeshT, true>>();
    else
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<ORDER, SelectedIntegral, MeshT, false>>();
  } else  // kElastic
  {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kElastic>>();
    else
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kElastic>>();
  }
}

/**
 * @brief Creates solver for unstructured mesh.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeSolverUnstruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  } else if (physic == feenum::physicType::kAcoustoElastic) {
    if (isModelOnNodes)
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<ORDER, SelectedIntegral, MeshT, true>>();
    else
      return std::make_unique<solver::fe::SEMsolverAcoustoElastic<ORDER, SelectedIntegral, MeshT, false>>();
  } else  // kElastic
  {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kElastic>>();
    else
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kElastic>>();
  }
}

/**
 * @brief Creates a SEM solver with the specified integral implementation.
 */
template <auto ImplTag>
std::unique_ptr<Solver> makeSemSolver(int order, feenum::meshType mesh, feenum::modelLocationType modelLocation,
                                      feenum::physicType physic) {
  bool const isModelOnNodes = (modelLocation == feenum::modelLocationType::kOnNodes);

  if (physic == feenum::physicType::kAcoustic) {
    return orderDispatch<MAX_SOLVER_ACOUSTIC_ORDER>(order, [&](auto orderIC) {
      constexpr int ORDER = decltype(orderIC)::value;
      return (mesh == feenum::meshType::kStruct) ? makeSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                                 : makeSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
    });
  } else if (physic == feenum::physicType::kElastic) {
    return orderDispatch<MAX_SOLVER_ELASTIC_ORDER>(order, [&](auto orderIC) {
      constexpr int ORDER = decltype(orderIC)::value;
      return (mesh == feenum::meshType::kStruct) ? makeSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                                 : makeSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
    });
  } else if (physic == feenum::physicType::kAcoustoElastic) {
    return orderDispatch<MAX_SOLVER_ELASTOACOUSTIC_ORDER>(order, [&](auto orderIC) {
      constexpr int ORDER = decltype(orderIC)::value;
      return (mesh == feenum::meshType::kStruct) ? makeSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                                 : makeSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
    });
  }

  throw std::runtime_error("Unknown physics type");
}

std::unique_ptr<Solver> createSolver(feenum::methodType const methodType, feenum::implemType const implemType,
                                     feenum::meshType const mesh, feenum::modelLocationType const modelLocation,
                                     feenum::physicType const physicType, int const order) {
  if (methodType == feenum::methodType::kSem) {
    switch (implemType) {
      case feenum::implemType::kMakutu:
        return makeSemSolver<IntegralType::MAKUTU>(order, mesh, modelLocation, physicType);
      default:
        throw std::runtime_error("Unknown implementation type: " + std::to_string(static_cast<int>(implemType)));
    }
  }

  // Add DG or other methods as needed
  throw std::runtime_error(
      "Unsupported solver configuration: methodType=" + std::to_string(static_cast<int>(methodType)) + ", implemType=" +
      std::to_string(static_cast<int>(implemType)) + ", physicType=" + std::to_string(static_cast<int>(physicType)));
}

}  // namespace solver_factory
}  // namespace fe
}  // namespace solver
