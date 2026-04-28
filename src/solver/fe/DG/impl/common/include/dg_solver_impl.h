#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#include <data_type.h>

#include <array>
#include <cstdlib>

#include "Integrals.h"
#include "dg_solver.h"

namespace solver {
namespace fe {

//============================================================================
// Update Solution
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolution(const float& dt,
                                                                                            Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  updateFields(dt, myData);
  FENCE
}

//============================================================================
// outputSolutionValues - Output field values for diagnostics
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const VECTOR_REAL_VIEW& fieldGlobal, const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << fieldGlobal(m_mesh.globalNodeIndex(e, 0, 0, 0)) << endl;
}


//============================================================================
// applyRHSTerm
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyRHSTerm(int timeSample, float dt,
                                                                                          const DataType& data) {
  int nb_rhs_element = data.getRhsElement().extent(0);
  auto mesh_local = m_mesh;  // Capture mesh for lambda

  static constexpr const char* rhsNames[1] = {"rhs"};
  for (int f = 0; f < kNumRhs; ++f) {
    rhsTermGlobal[f] = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), rhsNames[f]);
  }

  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int i) {
        for (int z = 0; z < ORDER + 1; z++) {
          for (int y = 0; y < ORDER + 1; y++) {
            for (int x = 0; x < ORDER + 1; x++) {
              int localNodeId = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
              int nodeRHS = mesh_local.globalNodeIndex(data.getRhsElement()[i], x, y, z);

              for (int f = 0; f < kNumRhs; ++f) {
                float source = data.getRhsTerm(f)(i, timeSample) * data.getRhsWeights()(i, localNodeId);
                rhsTermGlobal[f](nodeRHS) -= source;
              }
            }
          }
        }
      });
}


//============================================================================
// updateFields - Time integration update
//============================================================================
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFields(float dt,
                                                                                          const DataType& data) {
  // Extract scalar constants to local variables
  float const dt_local = dt;
  float const dt2_local = dt * dt;

  int const kNumElem = mesh_local.getNumberOfElements();

  auto mesh_local = m_mesh;

  allocateArray2D(kNumElem, kPointsPerElement, "current_field");
  allocateArray2D(kNumElem, kPointsPerElement, "prev_field")
  current_field = data.getCurrentField(0);
  prev_field = data.getPreviousField(0);

  applyRHSTerm(timeSample, dt, data);
  
  Kokkos::parallel_for(
      "Solver DG Update Field Acoustic", kNumElem, KOKKOS_LAMBDA(const int e) {
    float massMatrixLocal[kPointsPerElement] = {0};
    float stiffnessMatrixLocal[kPointsPerElement] = {0};
    float elementCoords[8][3];
    {
      auto const eIdx = mesh_local.elementIndex(e);
      int I = 0;
      for (int kv = 0; kv < 2; ++kv)
        for (int jv = 0; jv < 2; ++jv)
          for (int iv = 0; iv < 2; ++iv)
            mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), elementCoords[I++]);
    }
    computeMassTerm(elementCoords, [&](const int j, const real_t val) {massMatrixLocal[j] += val; });
    computeStiffnessTerm(elementCoords, [&](const int, const int, const int){}, [&](const int i, const int j, const real_t val){
      stiffnessMatrixLocal[i] += val*prev_field[e][j];});
    
    for(int faceId=0; i<6; ++i) {

      if (!mesh_local.isBoundaryFace(f)) {

        real_t const gamma = 0.5;

        int const f = mesh_local.getGlobalFace(e, static_cast<model::CubicFace>(faceId)); 

        // Get corner coordinates of the face for integration
        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) {
            faceCoords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
          }
        }

        float normal[3];
        mesh_local.faceNormal(e, static_cast<model::CubicFace>(faceId), normal);

        int const neighbor_e = 0;

        computeInterfaceFluxTerm(faceId, faceCoords, [&](const int i, const int j, const int k, const real_t val) { 
          int const ei = mesh_local.getGlobalNodeFromFace(f, i);
          int const ej = mesh_local.getGlobalNodeFromFace(f, j);
          int const ej_perm = 0;
          int const ei_perm = 0;
          stiffnessMatrixLocal[ei] += 0.5*val*prev_field[e][ej]*normal[k] + 0.5*val*prev_field[neighbor_e][ej_perm]*normal[k];
          stiffnessMatrixLocal[ej] += 0.5*val*prev_field[e][ei]*normal[k] + 0.5*val*prev_field[neighbor_e][ei_perm]*normal[k];
        });

        for(int i=0; i<knumNodesPerFace; ++i) {
          int const ei = mesh_local.getGlobalNodeFromFace(f, i);// node index in element
          int const ei_perm = 0;
          stiffnessMatrixLocal[ei] += gamma * ( computeDampingTerm(i, faceCoords)*prev_field[e][ei]
           - computeDampingTerm(i, faceCoords)*prev_field[neighbor_e][ei_perm] );
        }

      }
    }

    for (int z = 0; z < ORDER + 1; z++) {
      for (int y = 0; y < ORDER + 1; y++) {
        for (int x = 0; x < ORDER + 1; x++) {
          int i = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
          int nodeGlobalIndex = mesh_local.globalNodeIndex(e, x, y, z);
          stiffnessMatrixLocal[i] += rhsTermGlobal[0](nodeGlobalIndex);

          current_field[e][i] = 1 / massMatrixLocal[i] * ( 2 * massMatrixLocal[i] * prev_field[e][i]
              - dt2_local * stiffnessMatrixLocal[i] - massMatrixLocal[i] * current_field[e][i]); // add damping    
        }
      }
    }
    
  });

}


}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_


