#ifndef FUNTIDES_SOLVER_FE_IMPL_PHYSICS_ELASTIC_INCLUDE_ELASTIC_SOLVER_H_
#define FUNTIDES_SOLVER_FE_IMPL_PHYSICS_ELASTIC_INCLUDE_ELASTIC_SOLVER_H_

/**
 * @file elastic_solver.h
 * @brief Named elastic solver type (physics-policy migration, step 3).
 *
 * Thin final subclass of SEMsolver with the physics bound to kElastic;
 * order, integral back-end, mesh type and model location stay template
 * parameters. The anisotropy (Iso/Vti/Tti) remains a runtime choice of
 * the configured model, exactly as with the underlying template. The name
 * exists for readability at call sites and in diagnostics; the factory
 * remains the only construction path.
 *
 * Include from the umbrella (sem_solver_impl.h), never directly: the
 * definition of SEMsolver must be complete first.
 */

#include "sem_solver.h"

namespace solver {
namespace fe {

/**
 * @brief The elastic SEM solver, as a named type.
 *
 * @tparam ORDER           polynomial order of the SEM elements.
 * @tparam INTEGRAL_TYPE   integral back-end (sum-factorization, GEMM, ...).
 * @tparam MESH_TYPE       structured or unstructured mesh model.
 * @tparam IS_MODEL_ON_NODES whether the medium properties live on nodes.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
class ElasticSolver final : public SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
                                             utils::enums::physicType::kElastic> {
 public:
  ElasticSolver() = default;
};

static_assert(std::is_base_of_v<SEMsolver<1, int, int, false, utils::enums::physicType::kElastic>,
                                ElasticSolver<1, int, int, false>>,
              "ElasticSolver must derive from the elastic SEMsolver");

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_PHYSICS_ELASTIC_INCLUDE_ELASTIC_SOLVER_H_
