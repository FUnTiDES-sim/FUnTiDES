#ifndef SRC_SOLVER_FE_IMPL_INCLUDE_SOLVERFACTORY_H_
#define SRC_SOLVER_FE_IMPL_INCLUDE_SOLVERFACTORY_H_

#include <fe/Integrals.hpp>
#include <memory>

#include "sem_enums.h"
#include "sem_solver_base.h"

using namespace solver::fe::enums;

namespace SolverFactory
{

/**
 * @brief Creates a SEM solver instance based on the specified configuration.
 *
 * This factory function creates the appropriate solver type based on the
 * method type, implementation type, mesh type, model location, physics type,
 * and polynomial order.
 *
 * @param methodType The numerical method (SEM or DG)
 * @param implemType The implementation backend (Makutu or Shiva)
 * @param meshType The mesh type (Struct or Unstruct)
 * @param modelLocation Where model parameters are stored (OnNodes or
 * OnElements)
 * @param physicType The physics type (Acoustic or Elastic)
 * @param order The polynomial order of spectral elements
 * @return A unique pointer to the created solver
 * @throws std::runtime_error if the configuration is unsupported
 */
std::unique_ptr<SEMSolverBase> createSolver(methodType methodType,
                                            implemType implemType,
                                            meshType meshType,
                                            modelLocationType modelLocation,
                                            physicType physicType, int order);

}  // namespace SolverFactory
#endif  // SRC_SOLVER_FE_IMPL_INCLUDE_SOLVERFACTORY_H_
