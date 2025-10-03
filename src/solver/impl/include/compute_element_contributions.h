#pragma once

#include <data_type.h>
#include <model.h>

namespace solver
{
namespace kernels
{


template < int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE >
void computeElementContributions( int i2, 
                                  const ARRAY_REAL_VIEW &pnGlobal, 
                                  const MESH_TYPE &m_mesh,
                                  VECTOR_REAL_VIEW & yGlobal );

} // namespace kernels
} // namespace solver