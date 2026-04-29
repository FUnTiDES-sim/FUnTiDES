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

  int const kNumElem = m_mesh.getNumberOfElements();
  m_rhs_elem_ = allocateArray2D<ARRAY_REAL_VIEW>(kNumElem, kPointsPerElement, "rhsElem");

  auto rhs_elem_view = m_rhs_elem_;
  auto rhs_element_view = data.getRhsElement();
  auto rhs_term_view = data.getRhsTerm(0);
  auto rhs_weights_view = data.getRhsWeights();

  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int s) {
        int const src_elem = rhs_element_view[s];
        float const wavelet_val = rhs_term_view(s, timeSample);
        for (int dof = 0; dof < kPointsPerElement; ++dof) {
          rhs_elem_view(src_elem, dof) -= wavelet_val * rhs_weights_view(s, dof);
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
  ARRAY_REAL_VIEW rhs_elem_local = m_rhs_elem_;
  real_t const penalty_local = m_penalty_factor_;

  int const kNumElem = mesh_local.getNumberOfElements();
  // SEM convention: current_field = p^n, prev_field = p^{n-1}; result written into prev_field
  ARRAY_REAL_VIEW current_field =  data.getCurrentField(0);
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

    real_t const vp = mesh_local.getModelVpOnElement(e);
    real_t const rho = mesh_local.getModelRhoOnElement(e);
    real_t const inv_model_factor = 1.0f / (vp * vp * rho);
    real_t const inv_rho = 1.0f / rho;

    INTEGRAL_TYPE::computeMassTerm(elementCoords,
                                   [&](const int j, const real_t val) { massMatrixLocal[j] += inv_model_factor * val; });
    // Stiffness uses p^n = current_field (SEM convention); scaled by 1/rho for acoustic wave equation
    INTEGRAL_TYPE::computeStiffnessTerm(elementCoords, [&](const int, const int, const int) {},
                                        [&](const int i, const int j, const real_t val) {
                                          stiffnessMatrixLocal[i] += inv_rho * val * current_field(e, j);
                                        });

    for (int faceId = 0; faceId < 6; ++faceId) {
      int const f = face_connectivity_local.getGlobalFace(e, static_cast<model::CubicFace>(faceId));

      if (!face_connectivity_local.isBoundaryFace(f)) {
        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni =
              face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float normal[3];
        mesh_local.faceNormal(e, static_cast<model::CubicFace>(faceId), normal);

        // Determine which element is "self" and which is "other": the
        // face_connectivity stores one fixed owner/neighbor assignment that is
        // independent of the current loop element e.
        bool const e_is_owner = (face_connectivity_local.elemOwner(f) == e);
        int const neighbor_e = e_is_owner ? face_connectivity_local.elemNeighbor(f)
                                           : face_connectivity_local.elemOwner(f);
        model::CubicFace const neighbor_local_face = static_cast<model::CubicFace>(
            e_is_owner ? face_connectivity_local.localFaceNeighbor(f)
                       : face_connectivity_local.localFaceOwner(f));
        real_t const gamma = computeSIPGPenalty<ORDER>(faceCoords, elementCoords, penalty_local);

        INTEGRAL_TYPE::computeInterfaceFluxTerm(faceCoords, elementCoords, faceId,
            [&](const int i, const int j, const int k, const real_t val) {
              int const ei = model::faceLocalToElemLocal(static_cast<model::CubicFace>(faceId), i, ORDER);
              int const ej = model::faceLocalToElemLocal(static_cast<model::CubicFace>(faceId), j, ORDER);
              int const ej_perm = model::faceLocalToElemLocal(
                  neighbor_local_face, face_connectivity_local.getNeighborFaceDof(f, j), ORDER);
              int const ei_perm = model::faceLocalToElemLocal(
                  neighbor_local_face, face_connectivity_local.getNeighborFaceDof(f, i), ORDER);
              // Symmetry term (stiffnessLocal[ei]) and consistency term (stiffnessLocal[ej])
              // both require 1/rho scaling for the acoustic wave equation.
              stiffnessMatrixLocal[ei] += inv_rho * (-0.5f * val * current_field(e, ej) * normal[k] +
                                                      0.5f * val * current_field(neighbor_e, ej_perm) * normal[k]);
              stiffnessMatrixLocal[ej] += inv_rho * (-0.5f * val * current_field(e, ei) * normal[k] +
                                                      0.5f * val * current_field(neighbor_e, ei_perm) * normal[k]);
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei = model::faceLocalToElemLocal(static_cast<model::CubicFace>(faceId), i, ORDER);
          int const ei_perm = model::faceLocalToElemLocal(
              neighbor_local_face, face_connectivity_local.getNeighborFaceDof(f, i), ORDER);
          stiffnessMatrixLocal[ei] +=
              gamma * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
              (current_field(e, ei) - current_field(neighbor_e, ei_perm));
        }
      }
    }

    // Verlet (SEM convention): write p^{n+1} into prev_field; source element-indexed
    for (int i = 0; i < kPointsPerElement; ++i) {
      stiffnessMatrixLocal[i] += rhs_elem_local(e, i);
      prev_field(e, i) = (2.0f * massMatrixLocal[i] * current_field(e, i) - dt2_local * stiffnessMatrixLocal[i] -
                          massMatrixLocal[i] * prev_field(e, i)) /
                         massMatrixLocal[i];
    }
  });

}


}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_


