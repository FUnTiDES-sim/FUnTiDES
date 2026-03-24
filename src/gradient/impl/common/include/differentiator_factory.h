#ifndef FUNTIDES_GRADIENT_IMPL_COMMON_INCLUDE_DIFFERENTIATOR_FACTORY_H_
#define FUNTIDES_GRADIENT_IMPL_COMMON_INCLUDE_DIFFERENTIATOR_FACTORY_H_

#include "sem_enums.h"
#include "differentiator.h"

namespace gradient
{

namespace feenum = solver::fe::enums;

/**
 * @brief Factory function to create a Differentiator instance with specified parameters.
 *
 * This function instantiates and returns a Differentiator object configured according to the
 * provided implementation, mesh, model location, and physics type specifications.
 *
 * @param implemType The implementation type that determines the computational backend or algorithm.
 * @param meshType The type of mesh used in the finite element analysis.
 * @param modelLocation The location within the model where computations are performed.
 * @param physicType The physics type governing the differentiation behavior.
 * @param order The order of differentiation or accuracy order for the numerical scheme.
 *
 * @return std::unique_ptr<Differentiator> A unique pointer to the created Differentiator object.
 *         The caller takes ownership of the returned object.
 *
 * @throws Potentially throws exceptions if invalid parameter combinations are provided or
 *         if memory allocation fails.
 */
std::unique_ptr<Differentiator> createDifferentiator(feenum::implemType implemType,
                                                     feenum::meshType meshType,
                                                     feenum::modelLocationType modelLocation,
                                                     feenum::physicType physicType, int order);

}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_COMMON_INCLUDE_DIFFERENTIATOR_FACTORY_H_
