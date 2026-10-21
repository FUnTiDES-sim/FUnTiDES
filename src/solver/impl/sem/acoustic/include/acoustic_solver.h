#ifndef FUNTIDES_SOLVER_IMPL_SEM_ACOUSTIC_INCLUDE_ACOUSTIC_SOLVER_H_
#define FUNTIDES_SOLVER_IMPL_SEM_ACOUSTIC_INCLUDE_ACOUSTIC_SOLVER_H_

/**
 * @file acoustic_solver.h
 * @brief Named acoustic solver type.
 *
 * Thin final subclass of SEMsolver with the physics bound to kAcoustic;
 * order, integral back-end, mesh type and model location stay template
 * parameters. The name exists for readability at call sites and in
 * diagnostics; the factory remains the only construction path.
 *
 * Include from the umbrella (sem_solver_impl.h), never directly: the
 * definition of SEMsolver must be complete first.
 */

#include "sem_solver.h"

namespace solver {

/**
 * @brief The acoustic SEM solver, as a named type.
 *
 * @tparam ORDER           polynomial order of the SEM elements.
 * @tparam INTEGRAL_TYPE   integral back-end (sum-factorization, GEMM, ...).
 * @tparam MESH_TYPE       structured or unstructured mesh model.
 * @tparam IS_MODEL_ON_NODES whether the medium properties live on nodes.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
class AcousticSolver final : public SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
                                              utils::enums::physicType::kAcoustic> {
 public:
  AcousticSolver() = default;
};

static_assert(std::is_base_of_v<SEMsolver<1, int, int, false, utils::enums::physicType::kAcoustic>,
                                AcousticSolver<1, int, int, false>>,
              "AcousticSolver must derive from the acoustic SEMsolver");

}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_SEM_ACOUSTIC_INCLUDE_ACOUSTIC_SOLVER_H_
