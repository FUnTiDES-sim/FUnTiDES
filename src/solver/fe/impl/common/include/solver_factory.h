#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SOLVER_FACTORY_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SOLVER_FACTORY_H_
#include <Integrals.h>

#include <memory>

#include "sem_enums.h"
#include "solver.h"

namespace solver {
namespace fe {
namespace solver_factory {

/**
 * @brief Creates the solver matching a runtime configuration.
 *
 * @param[in] methodType Numerical method.
 * @param[in] implemType Integral back-end.
 * @param[in] meshType Mesh kind (structured or unstructured).
 * @param[in] modelLocation Whether model parameters are stored on nodes or on elements.
 * @param[in] physicType Physics of the equation solved.
 * @param[in] order Polynomial order of the elements (the higher order for the p-adaptive method).
 * @param[in] order_min Lower polynomial order, used by the p-adaptive method only
 *            (0 < order_min < order); ignored by the other methods.
 * @return Owning pointer to the created solver.
 * @throws std::runtime_error if the configuration is unsupported.
 */
std::unique_ptr<Solver> createSolver(utils::enums::methodType methodType, utils::enums::implemType implemType,
                                     utils::enums::meshType meshType, utils::enums::modelLocationType modelLocation,
                                     utils::enums::physicType physicType, int const order, int const order_min = 0);
}  // namespace solver_factory
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SOLVER_FACTORY_H_
