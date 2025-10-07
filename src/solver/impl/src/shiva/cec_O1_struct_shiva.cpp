#include <model_struct.h>

#include "compute_element_contributions_impl.h"

#include "Integrals.hpp"

namespace solver
{
namespace kernels
{
constexpr int ORDER = 1;
template void 
computeElementContributions< ORDER,
                             IntegralTypeSelector<ORDER, true, IntegralType::SHIVA>::type,
                             model::ModelStruct<float, int, ORDER> > 
                            ( int i2, 
                               const ARRAY_REAL_VIEW &pnGlobal, 
                               const model::ModelStruct<float, int, ORDER> &m_mesh,
                               VECTOR_REAL_VIEW & yGlobal );

} // namespace solver
} // namespace kernels