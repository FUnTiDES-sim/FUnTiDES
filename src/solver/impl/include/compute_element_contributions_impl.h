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

  Kokkos::parallel_for( Kokkos::RangePolicy<Kokkos::LaunchBounds<32, 1>>(0, m_mesh.getNumberOfElements()), \
                        KOKKOS_LAMBDA ( const int elementNumber )
  {

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  constexpr int N = ORDER + 1;
  constexpr int N2 = N * N;
  
  typename MESH_TYPE::IndexType elemIndex;
  m_mesh.elementIndex(elementNumber, elemIndex);


  static constexpr int nPointsElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  float pnLocal[nPointsElement] = {0};
  float Y[nPointsElement] = {0};

  for( int k=0; k<N; ++k )
  {
    for( int j=0 ; j<N; ++j )
    {
      for( int i=0 ; i<N; ++i )
      {
        int const globalNodeIndex = m_mesh.globalNodeIndex( elemIndex, i, j, k);
        pnLocal[i + j * N + k * N2] = pnGlobal(globalNodeIndex, i2);
      }
    }
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
    int x = i % N;
    int z = (i / N) % N;
    int y = i / (N * N);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    ATOMICADD(yGlobal[gIndex], Y[i]);
  }

  });
}


} // namespace kernels
} // namespace solver