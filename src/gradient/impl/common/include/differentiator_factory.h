#ifndef FUNTIDES_GRADIENT_IMPL_COMMON_INCLUDE_DIFFERENTIATOR_FACTORY_H_
#define FUNTIDES_GRADIENT_IMPL_COMMON_INCLUDE_DIFFERENTIATOR_FACTORY_H_

#include "differentiator.h"
#include "sem_enums.h"

namespace gradient {

/**
 * @brief Creates the Differentiator matching the requested configuration.
 *
 * @param implemType Back-end selector.
 * @param meshType Mesh kind (structured or unstructured).
 * @param modelLocation Where the model parameters are stored.
 * @param physicType Physics of the differentiator (selects the wavefield and gradient types).
 * @param order Polynomial order of the spectral elements.
 * @todo VERIFY: is `order` the polynomial order of the elements, and which orders are supported?
 *
 * @return Newly created differentiator; the caller owns it.
 *
 * @throws @todo VERIFY: which exception type, and for which unsupported combinations of arguments?
 */
std::unique_ptr<Differentiator> createDifferentiator(utils::enums::implemType implemType,
                                                     utils::enums::meshType meshType,
                                                     utils::enums::modelLocationType modelLocation,
                                                     utils::enums::physicType physicType, int order);

}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_COMMON_INCLUDE_DIFFERENTIATOR_FACTORY_H_
