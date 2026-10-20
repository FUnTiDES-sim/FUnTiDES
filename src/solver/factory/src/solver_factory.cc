#include "solver_factory.h"

#include <Integrals.h>
#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver.h"
#include "sem_solver_acoustoelastic.h"
#ifdef COMPILE_DG
#include "dg_solver.h"
#endif
#ifdef COMPILE_DG_SEM
#include "dgsem_solver.h"
#endif
#ifdef COMPILE_DG_PADAPTIVE
#include "dg_padaptive_solver.h"
#endif

namespace solver {
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

// ============================================================================
// SEM solver (acoustic, elastic or acousto-elastic)
// ============================================================================

/**
 * @brief Creates a SEM-type solver for one mesh type.
 *
 * Any physic other than kAcoustic and kAcoustoElastic selects the elastic solver.
 *
 * @tparam MeshT Mesh type; structured meshes carry the polynomial order as a
 *         template argument, unstructured meshes do not.
 */
template <auto ImplTag, int ORDER, typename MeshT>
std::unique_ptr<Solver> makeSemSolverForMesh(bool isModelOnNodes, feenum::physicType physic) {
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  } else if (physic == feenum::physicType::kAcoustoElastic) {
    if (isModelOnNodes)
      return std::make_unique<solver::SEMsolverAcoustoElastic<ORDER, SelectedIntegral, MeshT, true>>();
    else
      return std::make_unique<solver::SEMsolverAcoustoElastic<ORDER, SelectedIntegral, MeshT, false>>();
  } else  // any other value selects elastic
  {
    if (isModelOnNodes)
      return std::make_unique<
          solver::SEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kElastic>>();
    else
      return std::make_unique<
          solver::SEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kElastic>>();
  }
}

/**
 * @brief Creates a SEM-type solver of order ORDER, selecting the mesh type.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeSemSolverForOrder(bool isModelOnNodes, feenum::meshType mesh, feenum::physicType physic) {
  if (mesh == feenum::meshType::kStruct) {
    return makeSemSolverForMesh<ImplTag, ORDER, model::ModelStruct<float, int, ORDER>>(isModelOnNodes, physic);
  }
  return makeSemSolverForMesh<ImplTag, ORDER, model::ModelUnstruct<float, int>>(isModelOnNodes, physic);
}

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
      return makeSemSolverForOrder<ImplTag, decltype(orderIC)::value>(isModelOnNodes, mesh, physic);
    });
  } else if (physic == feenum::physicType::kElastic) {
    return orderDispatch<MAX_SOLVER_ELASTIC_ORDER>(order, [&](auto orderIC) {
      return makeSemSolverForOrder<ImplTag, decltype(orderIC)::value>(isModelOnNodes, mesh, physic);
    });
  } else if (physic == feenum::physicType::kAcoustoElastic) {
    return orderDispatch<MAX_SOLVER_ELASTOACOUSTIC_ORDER>(order, [&](auto orderIC) {
      return makeSemSolverForOrder<ImplTag, decltype(orderIC)::value>(isModelOnNodes, mesh, physic);
    });
  }

  throw std::runtime_error("Unknown physics type");
}

#ifdef COMPILE_DG
// ============================================================================
// DG solver (acoustic only)
// ============================================================================

/**
 * @brief Creates a DG solver for one mesh type.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER, typename MeshT>
std::unique_ptr<Solver> makeDgSolverForMesh(bool isModelOnNodes, feenum::physicType physic) {
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::DGsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::DGsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG: unsupported physics type");
}

/**
 * @brief Creates a DG solver of order ORDER, selecting the mesh type.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeDgSolverForOrder(bool isModelOnNodes, feenum::meshType mesh, feenum::physicType physic) {
  if (mesh == feenum::meshType::kStruct) {
    return makeDgSolverForMesh<ImplTag, ORDER, model::ModelStruct<float, int, ORDER>>(isModelOnNodes, physic);
  }
  return makeDgSolverForMesh<ImplTag, ORDER, model::ModelUnstruct<float, int>>(isModelOnNodes, physic);
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
    return makeDgSolverForOrder<ImplTag, decltype(orderIC)::value>(isModelOnNodes, mesh, physic);
  });
}
#endif

#ifdef COMPILE_DG_SEM
// ============================================================================
// Coupled DG-SEM solver (acoustic only)
// ============================================================================

/**
 * @brief Creates a coupled DG-SEM solver for one mesh type.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER, typename MeshT>
std::unique_ptr<Solver> makeDgSemSolverForMesh(bool isModelOnNodes, feenum::physicType physic) {
  using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;

  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<
          solver::DGSEMsolver<ORDER, SelectedIntegral, MeshT, true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<
          solver::DGSEMsolver<ORDER, SelectedIntegral, MeshT, false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG-SEM: unsupported physics type");
}

/**
 * @brief Creates a coupled DG-SEM solver of order ORDER, selecting the mesh type.
 */
template <auto ImplTag, int ORDER>
std::unique_ptr<Solver> makeDgSemSolverForOrder(bool isModelOnNodes, feenum::meshType mesh, feenum::physicType physic) {
  if (mesh == feenum::meshType::kStruct) {
    return makeDgSemSolverForMesh<ImplTag, ORDER, model::ModelStruct<float, int, ORDER>>(isModelOnNodes, physic);
  }
  return makeDgSemSolverForMesh<ImplTag, ORDER, model::ModelUnstruct<float, int>>(isModelOnNodes, physic);
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
    return makeDgSemSolverForOrder<ImplTag, decltype(orderIC)::value>(isModelOnNodes, mesh, physic);
  });
}
#endif

#ifdef COMPILE_DG_PADAPTIVE
// ============================================================================
// DG p-adaptive solver (acoustic only)
// ============================================================================

/**
 * @brief Creates a DG p-adaptive solver for one mesh type.
 *
 * The structured mesh type carries ORDER_MAX as its polynomial order.
 * @throws std::runtime_error If physic is not kAcoustic.
 */
template <auto ImplTag, int ORDER_MIN, int ORDER_MAX, typename MeshT>
std::unique_ptr<Solver> makeDgPAdaptiveSolverForMesh(bool isModelOnNodes, feenum::physicType physic) {
  if (physic == feenum::physicType::kAcoustic) {
    if (isModelOnNodes)
      return std::make_unique<solver::DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, IntegralTypeSelector, ImplTag, MeshT,
                                                            true, feenum::physicType::kAcoustic>>();
    else
      return std::make_unique<solver::DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, IntegralTypeSelector, ImplTag, MeshT,
                                                            false, feenum::physicType::kAcoustic>>();
  }
  throw std::runtime_error("DG p-adaptive: unsupported physics type");
}

/**
 * @brief Creates a DG p-adaptive solver of orders (ORDER_MIN, ORDER_MAX), selecting the mesh type.
 */
template <auto ImplTag, int ORDER_MIN, int ORDER_MAX>
std::unique_ptr<Solver> makeDgPAdaptiveSolverForOrders(bool isModelOnNodes, feenum::meshType mesh,
                                                        feenum::physicType physic) {
  if (mesh == feenum::meshType::kStruct) {
    return makeDgPAdaptiveSolverForMesh<ImplTag, ORDER_MIN, ORDER_MAX, model::ModelStruct<float, int, ORDER_MAX>>(
        isModelOnNodes, physic);
  }
  return makeDgPAdaptiveSolverForMesh<ImplTag, ORDER_MIN, ORDER_MAX, model::ModelUnstruct<float, int>>(
      isModelOnNodes, physic);
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
            return makeDgPAdaptiveSolverForOrders<ImplTag, decltype(orderMinIC)::value, ORDER_MAX>(isModelOnNodes, mesh,
                                                                                                   physic);
          });
        } else {
          throw std::runtime_error("DG p-adaptive requires order_max > 1");
        }
      });
}
#endif

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
      "Unsupported solver configuration: methodType=" + std::to_string(static_cast<int>(methodType)) +
      ", implemType=" + std::to_string(static_cast<int>(implemType)) + ", physicType=" +
      std::to_string(static_cast<int>(physicType)));
}

}  // namespace solver_factory
}  // namespace solver
