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
                                  VECTOR_REAL_VIEW & yGlobal )
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;
  static constexpr int nPointsElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float pnLocal[nPointsElement] = {0};
  float Y[nPointsElement] = {0};

  constexpr int dim = ORDER + 1;
  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const globalIdx = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    pnLocal[i] = pnGlobal(globalIdx, i2);
  }

  typename INTEGRAL_TYPE::TransformType transformData;
  INTEGRAL_TYPE::gatherCoordinates(elementNumber, m_mesh, transformData);

  // Stiffness term
  INTEGRAL_TYPE::computeStiffnessTerm(
      transformData, [&](const int i, const int j, const real_t val) {
        float localIncrement = val * pnLocal[j];
        Y[i] += localIncrement;
      });

  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    ATOMICADD(yGlobal[gIndex], Y[i]);
  }

  MAINLOOPEND
}


} // namespace kernels
} // namespace solver