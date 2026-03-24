#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SOLVER_FACTORY_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SOLVER_FACTORY_H_
#include <fe/Integrals.hpp>
#include <memory>

#include "sem_enums.h"
#include "solver.h"

namespace solver
{
namespace fe
{
namespace solver_factory
{

/**
 * @brief Creates a SEM solver instance based on the specified configuration.
 *
 * This factory function creates the appropriate solver type based on the
 * method type, implementation type, mesh type, model location, physics type,
 * and polynomial order.
 *
 * @param methodType The numerical method (SEM or DG)
 * @param implemType The implementation backend (Makutu)
 * @param meshType The mesh type (Struct or Unstruct)
 * @param modelLocation Where model parameters are stored (OnNodes or
 * OnElements)
 * @param physicType The physics type (Acoustic or Elastic)
 * @param order The polynomial order of spectral elements
 * @return A unique pointer to the created solver
 * @throws std::runtime_error if the configuration is unsupported
 */
std::unique_ptr<Solver> createSolver(utils::enums::methodType methodType,
                                     utils::enums::implemType implemType,
                                     utils::enums::meshType meshType,
                                     utils::enums::modelLocationType modelLocation,
                                     utils::enums::physicType physicType, int order);
}  // namespace solver_factory
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SOLVER_FACTORY_H_
