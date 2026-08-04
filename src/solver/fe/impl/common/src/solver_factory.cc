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

#ifdef COMPILE_DG
/**
 * @brief Creates DG solver for structured mesh (acoustic only).
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
 * @brief Creates DG solver for unstructured mesh (acoustic only).
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
 * @brief Creates a DG solver with the specified integral implementation.
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
 * @brief Creates DG-SEM coupled solver for structured mesh (acoustic only).
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
 * @brief Creates DG-SEM coupled solver for unstructured mesh (acoustic only).
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
 * @brief Creates a DG-SEM coupled solver with the specified integral implementation.
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
 * @brief Creates DG p-adaptive solver for structured mesh (acoustic only).
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
 * @brief Creates DG p-adaptive solver for unstructured mesh (acoustic only).
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
 * @brief Creates a DG p-adaptive solver with the specified integral implementation.
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
