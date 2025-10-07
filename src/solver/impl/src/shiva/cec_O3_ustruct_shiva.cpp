#include <model_unstruct.h>

#include "compute_element_contributions_impl.h"

#include "Integrals.hpp"

namespace solver
{
namespace kernels
{
constexpr int ORDER = 3;
template void 
computeElementContributions< ORDER,
                             IntegralTypeSelector<ORDER, false, IntegralType::SHIVA>::type,
                             model::ModelUnstruct<float, int> > 
                            ( int i2, 
                               const ARRAY_REAL_VIEW &pnGlobal, 
                               const model::ModelUnstruct<float, int> &m_mesh,
                               VECTOR_REAL_VIEW & yGlobal );

} // namespace solver
} // namespace kernels