#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_IMPL_H_
#include <data_type.h>

#include <array>
#include <cstdlib>

#include "Integrals.h"
#include "dg_penalty.h"
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const ARRAY_REAL_VIEW& fieldGlobal, const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << fieldGlobal(e, 0) << endl;
}


//============================================================================
// computeFEInit - Initialize mesh and face connectivity
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& /*sponge_size*/,
    const bool /*surface_sponge*/, const float /*taper_delta*/) {
  if (auto* typed_mesh = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh = *typed_mesh;
  } else {
    throw std::runtime_error("Incompatible mesh type in DG solver");
  }
  m_face_connectivity_.build(m_mesh);
}

//============================================================================
// Compute Forces (Phase 1)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeForces(const float& dt,
                                                                                           const int& timeSample,
                                                                                           Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  applyRHSTerm(timeSample, dt, myData);
  FENCE
}

//============================================================================
// applyRHSTerm
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyRHSTerm(int timeSample, float dt,
                                                                                          const DataType& data) {
  int nb_rhs_element = data.getRhsElement().extent(0);
  auto mesh_local = m_mesh;  // Capture mesh for lambda

  static constexpr const char* rhsNames[1] = {"rhs"};
  for (int f = 0; f < kNumRhs; ++f) {
    rhsTermGlobal[f] = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), rhsNames[f]);
  }

  auto rhs_views = rhsTermGlobal;

  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int i) {
        for (int z = 0; z < ORDER + 1; z++) {
          for (int y = 0; y < ORDER + 1; y++) {
            for (int x = 0; x < ORDER + 1; x++) {
              int localNodeId = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
              int nodeRHS = mesh_local.globalNodeIndex(data.getRhsElement()[i], x, y, z);

              for (int f = 0; f < kNumRhs; ++f) {
                float source = data.getRhsTerm(f)(i, timeSample) * data.getRhsWeights()(i, localNodeId);
                rhs_views[f](nodeRHS) -= source;
              }
            }
          }
        }
      });
}


//============================================================================
// updateFields - Time integration update
//============================================================================
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFields(float dt,
                                                                                          const DataType& data) {
  float const dt2_local = dt * dt;

  auto mesh_local = m_mesh;
  auto face_connectivity_local = m_face_connectivity_;
  VECTOR_REAL_VIEW rhs_global_0 = rhsTermGlobal[0];
  real_t const penalty_local = m_penalty_factor_;

  int const kNumElem = mesh_local.getNumberOfElements();
  ARRAY_REAL_VIEW current_field = data.getCurrentField(0);
  ARRAY_REAL_VIEW prev_field = data.getPreviousField(0);

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
    INTEGRAL_TYPE::computeMassTerm(elementCoords, [&](const int j, const real_t val) { massMatrixLocal[j] += val; });
    INTEGRAL_TYPE::computeStiffnessTerm(elementCoords, [&](const int, const int, const int) {},
                                        [&](const int i, const int j, const real_t val) {
                                          stiffnessMatrixLocal[i] += val * prev_field(e, j);
                                        });

    for (int faceId = 0; faceId < 6; ++faceId) {
      int const f = face_connectivity_local.getGlobalFace(e, static_cast<model::CubicFace>(faceId));

      if (!face_connectivity_local.isBoundaryFace(f)) {
        // Get corner coordinates of the face for integration
        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const globalNodeIndex = face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) {
            faceCoords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
          }
        }

        float normal[3];
        mesh_local.faceNormal(e, static_cast<model::CubicFace>(faceId), normal);

        int const neighbor_e = face_connectivity_local.elemNeighbor(f);

        model::CubicFace const neighbor_local_face =
            static_cast<model::CubicFace>(face_connectivity_local.localFaceNeighbor(f));

        real_t const gamma = computeSIPGPenalty<ORDER>(faceCoords, elementCoords, penalty_local);

        INTEGRAL_TYPE::computeInterfaceFluxTerm(faceCoords, elementCoords, faceId, [&](const int i, const int j, const int k, const real_t val) {
          int const ei = model::faceLocalToElemLocal(static_cast<model::CubicFace>(faceId), i, ORDER);
          int const ej = model::faceLocalToElemLocal(static_cast<model::CubicFace>(faceId), j, ORDER);
          int const ej_perm =
              model::faceLocalToElemLocal(neighbor_local_face, face_connectivity_local.getNeighborFaceDof(f, j), ORDER);
          int const ei_perm =
              model::faceLocalToElemLocal(neighbor_local_face, face_connectivity_local.getNeighborFaceDof(f, i), ORDER);
          stiffnessMatrixLocal[ei] +=
              -0.5f * val * prev_field(e, ej) * normal[k] + 0.5f * val * prev_field(neighbor_e, ej_perm) * normal[k];
          stiffnessMatrixLocal[ej] +=
              -0.5f * val * prev_field(e, ei) * normal[k] - 0.5f * val * prev_field(neighbor_e, ei_perm) * normal[k];
        });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei = model::faceLocalToElemLocal(static_cast<model::CubicFace>(faceId), i, ORDER);
          int const ei_perm =
              model::faceLocalToElemLocal(neighbor_local_face, face_connectivity_local.getNeighborFaceDof(f, i), ORDER);
          stiffnessMatrixLocal[ei] += gamma * INTEGRAL_TYPE::computeDampingTerm(i,faceCoords)*(prev_field(e, ei) - prev_field(neighbor_e, ei_perm));
        }
      }
    }

    for (int z = 0; z < ORDER + 1; z++) {
      for (int y = 0; y < ORDER + 1; y++) {
        for (int x = 0; x < ORDER + 1; x++) {
          int i = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
          int nodeGlobalIndex = mesh_local.globalNodeIndex(e, x, y, z);
          stiffnessMatrixLocal[i] += rhs_global_0(nodeGlobalIndex);

          current_field(e, i) = (2.0f * massMatrixLocal[i] * prev_field(e, i) - dt2_local * stiffnessMatrixLocal[i] -
                               massMatrixLocal[i] * current_field(e, i)) /
                              massMatrixLocal[i];  // add damping and source term
        }
      }
    }
    
  });

}


}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_


