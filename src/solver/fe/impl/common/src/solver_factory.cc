#include "solver_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver.h"
#include "sem_solver_acoustoelastic.h"
#ifdef COMPILE_DG
#include "dg_solver.h"
#endif
#ifdef COMPILE_DG_SEM
#include "dg-sem_solver.h"
#endif
#ifdef COMPILE_DG_PADAPTIVE
#include "dg_padaptive_solver.h"
#endif

namespace solver {
namespace fe {
namespace solver_factory {

namespace feenum = utils::enums;

/**
 * @brief Maps a runtime polynomial order to a compile-time one.
 *
 * Tries CurrentOrder, then CurrentOrder - 1, down to 1.
 *
 * @tparam CurrentOrder Highest order tried first; must be at least 1.
 * @param[in] order Runtime polynomial order.
 * @param[in] func Callable taking a std::integral_constant<int, ORDER> and returning a solver.
 * @return The solver built by func for the matching order.
 * @throws std::runtime_error If order is not in [1, CurrentOrder].
 */
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
 * @brief Creates a SEM-type solver (acoustic, elastic or acousto-elastic) for a structured mesh.
 *
 * Any physic other than kAcoustic and kAcoustoElastic selects the elastic solver.
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
  } else  // any other value selects elastic
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
 * @brief Creates a SEM-type solver (acoustic, elastic or acousto-elastic) for an unstructured mesh.
 *
 * Any physic other than kAcoustic and kAcoustoElastic selects the elastic solver.
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
  } else  // any other value selects elastic
  {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kElastic>>();
    else
      return std::make_unique<
          solver::fe::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kElastic>>();
  }
}

#ifdef COMPILE_DG
/**
 * @brief Creates a DG solver for a structured mesh.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeDgSolverStruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::DGsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::fe::DGsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG: unsupported physics type");
}

/**
 * @brief Creates a DG solver for an unstructured mesh.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeDgSolverUnstruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::DGsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::fe::DGsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG: unsupported physics type");
}

/**
 * @brief Creates a DG solver, dispatching the runtime order and mesh type to template arguments.
 * @throws std::runtime_error If the order or the physics is unsupported.
 */
template <auto ImplTag>
std::unique_ptr<Solver> makeDgSolver(int order, feenum::meshType mesh, feenum::modelLocationType modelLocation,
                                     feenum::physicType physic) {
  bool const isModelOnNodes = (modelLocation == feenum::modelLocationType::kOnNodes);
  return orderDispatch<MAX_DG_SOLVER_ACOUSTIC_ORDER>(order, [&](auto orderIC) {
    constexpr int ORDER = decltype(orderIC)::value;
    return (mesh == feenum::meshType::kStruct) ? makeDgSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                               : makeDgSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
  });
}
#endif

#ifdef COMPILE_DG_SEM
/**
 * @brief Creates a coupled DG-SEM solver for a structured mesh.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeDgSemSolverStruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelStruct<float, int, ORDER>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::DGSEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::fe::DGSEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG-SEM: unsupported physics type");
}

/**
 * @brief Creates a coupled DG-SEM solver for an unstructured mesh.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeDgSemSolverUnstruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelUnstruct<float, int>;
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::fe::DGSEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::fe::DGSEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG-SEM: unsupported physics type");
}

/**
 * @brief Creates a coupled DG-SEM solver, dispatching the runtime order and mesh type to template arguments.
 * @throws std::runtime_error If the order or the physics is unsupported.
 */
template <auto ImplTag>
std::unique_ptr<Solver> makeDgSemSolver(int order, feenum::meshType mesh, feenum::modelLocationType modelLocation,
                                        feenum::physicType physic) {
  bool const isModelOnNodes = (modelLocation == feenum::modelLocationType::kOnNodes);
  return orderDispatch<MAX_DG_SEM_SOLVER_ACOUSTIC_ORDER>(order, [&](auto orderIC) {
    constexpr int ORDER = decltype(orderIC)::value;
    return (mesh == feenum::meshType::kStruct) ? makeDgSemSolverStruct<ImplTag, ORDER>(isModelOnNodes, physic)
                                               : makeDgSemSolverUnstruct<ImplTag, ORDER>(isModelOnNodes, physic);
  });
}
#endif

#ifdef COMPILE_DG_PADAPTIVE
/**
 * @brief Creates a DG p-adaptive solver for a structured mesh.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER_MIN, int ORDER_MAX>
std::unique_ptr<Solver> makeDgPAdaptiveSolverStruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelStruct<float, int, ORDER_MAX>;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<solver::fe::DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, IntegralTypeSelector, ImplTag, MeshT,
                                                            true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<solver::fe::DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, IntegralTypeSelector, ImplTag, MeshT,
                                                            false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG p-adaptive: unsupported physics type");
}

/**
 * @brief Creates a DG p-adaptive solver for an unstructured mesh.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER_MIN, int ORDER_MAX>
std::unique_ptr<Solver> makeDgPAdaptiveSolverUnstruct(bool isModelOnNodes, feenum::physicType physic) {
  using MeshT = model::ModelUnstruct<float, int>;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<solver::fe::DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, IntegralTypeSelector, ImplTag, MeshT,
                                                            true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<solver::fe::DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, IntegralTypeSelector, ImplTag, MeshT,
                                                            false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG p-adaptive: unsupported physics type");
}

/**
 * @brief Creates a DG p-adaptive solver, dispatching the runtime orders and mesh type to template arguments.
 *
 * ORDER_MAX is dispatched over [1, MAX_DG_PADAPTIVE_SOLVER_ACOUSTIC_ORDER], then ORDER_MIN over
 * [1, ORDER_MAX - 1]. One explicit instantiation must exist per ordered pair (see
 * generate_padaptive_solver_implementations()).
 *
 * @throws std::runtime_error If the orders or the physics are unsupported.
 */
template <auto ImplTag>
std::unique_ptr<Solver> makeDgPAdaptiveSolver(int order_min, int order_max, feenum::meshType mesh,
                                              feenum::modelLocationType modelLocation, feenum::physicType physic) {
  bool const isModelOnNodes = (modelLocation == feenum::modelLocationType::kOnNodes);
  return orderDispatch<MAX_DG_PADAPTIVE_SOLVER_ACOUSTIC_ORDER>(
      order_max, [&](auto orderMaxIC) -> std::unique_ptr<Solver> {
        constexpr int ORDER_MAX = decltype(orderMaxIC)::value;

        if constexpr (ORDER_MAX > 1) {
          return orderDispatch<ORDER_MAX - 1>(order_min, [&](auto orderMinIC) {
            constexpr int ORDER_MIN = decltype(orderMinIC)::value;
            return (mesh == feenum::meshType::kStruct)
                       ? makeDgPAdaptiveSolverStruct<ImplTag, ORDER_MIN, ORDER_MAX>(isModelOnNodes, physic)
                       : makeDgPAdaptiveSolverUnstruct<ImplTag, ORDER_MIN, ORDER_MAX>(isModelOnNodes, physic);
          });
        } else {
          throw std::runtime_error("DG p-adaptive requires order_max > 1");
        }
      });
}
#endif

/**
 * @brief Creates a SEM-type solver, dispatching the runtime order, mesh type and physics to template arguments.
 *
 * The maximum supported order depends on the physics.
 *
 * @throws std::runtime_error If the order or the physics is unsupported.
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

/**
 * @brief Creates the solver matching a run configuration.
 *
 * DG, DG-SEM and p-adaptive solvers are only available when the corresponding COMPILE_* macro is defined.
 *
 * @param[in] order Polynomial order; for the p-adaptive method, the maximum order.
 * @param[in] order_min Minimum polynomial order; only used by the p-adaptive method, where it must satisfy
 *            0 < order_min < order.
 * @return The newly created solver.
 * @throws std::runtime_error If the configuration is not supported.
 */
std::unique_ptr<Solver> createSolver(feenum::methodType const methodType, feenum::implemType const implemType,
                                     feenum::meshType const mesh, feenum::modelLocationType const modelLocation,
                                     feenum::physicType const physicType, int const order, int const order_min) {
  if (methodType == feenum::methodType::kSem) {
    switch (implemType) {
      case feenum::implemType::kMakutu:
        return makeSemSolver<IntegralType::MAKUTU>(order, mesh, modelLocation, physicType);
      default:
        throw std::runtime_error("Unknown implementation type: " + std::to_string(static_cast<int>(implemType)));
    }
  }

#ifdef COMPILE_DG
  if (methodType == feenum::methodType::kDg) {
    switch (implemType) {
      case feenum::implemType::kMakutu:
        return makeDgSolver<IntegralType::MAKUTU>(order, mesh, modelLocation, physicType);
      default:
        throw std::runtime_error("Unknown DG implementation type: " + std::to_string(static_cast<int>(implemType)));
    }
  }
#endif

#ifdef COMPILE_DG_SEM
  if (methodType == feenum::methodType::kDgSem) {
    switch (implemType) {
      case feenum::implemType::kMakutu:
        return makeDgSemSolver<IntegralType::MAKUTU>(order, mesh, modelLocation, physicType);
      default:
        throw std::runtime_error("Unknown DG-SEM implementation type: " + std::to_string(static_cast<int>(implemType)));
    }
  }
#endif

#ifdef COMPILE_DG_PADAPTIVE
  if (methodType == feenum::methodType::kDgPAdaptive) {
    int const order_max = order;
    switch (implemType) {
      case feenum::implemType::kMakutu:
        if (order_min <= 0 || order_min >= order_max)
          throw std::runtime_error("DG p-adaptive requires 0 < order_min < order_max");
        return makeDgPAdaptiveSolver<IntegralType::MAKUTU>(order_min, order_max, mesh, modelLocation, physicType);
      default:
        throw std::runtime_error("Unknown DG p-adaptive implementation type: " +
                                 std::to_string(static_cast<int>(implemType)));
    }
  }
#endif

  throw std::runtime_error(
      "Unsupported solver configuration: methodType=" + std::to_string(static_cast<int>(methodType)) + ", implemType=" +
      std::to_string(static_cast<int>(implemType)) + ", physicType=" + std::to_string(static_cast<int>(physicType)));
}

}  // namespace solver_factory
}  // namespace fe
}  // namespace solver
