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
// outputSolutionValues
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const ARRAY_REAL_VIEW& fieldGlobal, const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << fieldGlobal(e, 0) << endl;
}

//============================================================================
// computeFEInit - Initialize mesh, face connectivity, and persistent arrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& /*sponge_size*/, const bool /*surface_sponge*/,
    const float /*taper_delta*/) {
  if (auto* typed_mesh = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh = *typed_mesh;
  } else {
    throw std::runtime_error("Incompatible mesh type in DG solver");
  }
  m_face_connectivity_.build(m_mesh);
  int const kNumElem = m_mesh.getNumberOfElements();
  m_rhs_elem_    = allocateArray2D<ARRAY_REAL_VIEW>(kNumElem, kPointsPerElement, "rhsElem");
  m_mass_local_  = allocateArray2D<ARRAY_REAL_VIEW>(kNumElem, kPointsPerElement, "massLocal");
  m_stiff_local_ = allocateArray2D<ARRAY_REAL_VIEW>(kNumElem, kPointsPerElement, "stiffLocal");
  m_damp_local_  = allocateArray2D<ARRAY_REAL_VIEW>(kNumElem, kPointsPerElement, "dampLocal");
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyRHSTerm(int timeSample, float /*dt*/,
                                                                                         const DataType& data) {
  int const nb_rhs_element  = data.getRhsElement().extent(0);
  auto rhs_elem_view        = m_rhs_elem_;
  auto rhs_element_view     = data.getRhsElement();
  auto rhs_term_view        = data.getRhsTerm(0);
  auto rhs_weights_view     = data.getRhsWeights();

  // Direct assignment (not +=): m_rhs_elem_ is zero-initialized at allocation and
  // non-source entries are never touched, so only source entries need overwriting.
  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int s) {
        int const src_elem  = rhs_element_view[s];
        float const wavelet_val = rhs_term_view(s, timeSample);
        for (int dof = 0; dof < kPointsPerElement; ++dof) {
          rhs_elem_view(src_elem, dof) = -wavelet_val * rhs_weights_view(s, dof);
        }
      });
}

//============================================================================
// computeVolumeAndBoundary - Kernel 1
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeVolumeAndBoundary(
    int kNumElem, ARRAY_REAL_VIEW current_field) {
  auto mesh_local              = m_mesh;
  auto face_connectivity_local = m_face_connectivity_;
  auto const face_to_elem_dof  = kFaceToElemDof;  // local copy for device capture
  ARRAY_REAL_VIEW mass_local_view  = m_mass_local_;
  ARRAY_REAL_VIEW stiff_local_view = m_stiff_local_;
  ARRAY_REAL_VIEW damp_local_view  = m_damp_local_;

  Kokkos::parallel_for(
      "DG Volume+Boundary", kNumElem, KOKKOS_LAMBDA(const int e) {
        float massLocal[kPointsPerElement]  = {0};
        float stiffLocal[kPointsPerElement] = {0};
        float dampLocal[kPointsPerElement]  = {0};
        float elementCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(e);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), elementCoords[I++]);
        }

        real_t const vp               = mesh_local.getModelVpOnElement(e);
        real_t const rho              = mesh_local.getModelRhoOnElement(e);
        real_t const inv_model_factor = 1.0f / (vp * vp * rho);
        real_t const inv_rho          = 1.0f / rho;
        real_t const inv_vp           = 1.0f / vp;

        INTEGRAL_TYPE::computeMassTerm(
            elementCoords, [&](const int j, const real_t val) { massLocal[j] += inv_model_factor * val; });

        real_t p_local[kPointsPerElement];
        for (int i = 0; i < kPointsPerElement; ++i) p_local[i] = current_field(e, i);
        INTEGRAL_TYPE::computeStiffnessTermSumFact(
            elementCoords, p_local, stiffLocal,
            [&](const int , const int , const int ) -> real_t { return inv_rho; });

        for (int faceId = 0; faceId < 6; ++faceId) {
          int const f = face_connectivity_local.getGlobalFace(e, static_cast<model::CubicFace>(faceId));
          if (!face_connectivity_local.isBoundaryFace(f)) continue;
          float faceCoords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const gni =
                face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
          }
          for (int i = 0; i < knumNodesPerFace; ++i) {
            int const ei = face_to_elem_dof[faceId][i];
            dampLocal[ei] += inv_vp * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords);
          }
        }

        for (int i = 0; i < kPointsPerElement; ++i) {
          mass_local_view(e, i)  = massLocal[i];
          stiff_local_view(e, i) = stiffLocal[i];
          damp_local_view(e, i)  = dampLocal[i];
        }
      });
}

//============================================================================
// computeInterfaceFlux - Kernel 2
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeInterfaceFlux(
    int kNumElem, ARRAY_REAL_VIEW current_field) {
  auto mesh_local              = m_mesh;
  auto face_connectivity_local = m_face_connectivity_;
  auto const face_to_elem_dof  = kFaceToElemDof;  // local copy for device capture
  ARRAY_REAL_VIEW stiff_local_view = m_stiff_local_;
  real_t const penalty_local       = m_penalty_factor_;

  Kokkos::parallel_for(
      "DG Interface Flux", kNumElem, KOKKOS_LAMBDA(const int e) {
        float elementCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(e);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), elementCoords[I++]);
        }

        real_t const rho     = mesh_local.getModelRhoOnElement(e);
        real_t const inv_rho = 1.0f / rho;

        for (int faceId = 0; faceId < 6; ++faceId) {
          int const f = face_connectivity_local.getGlobalFace(e, static_cast<model::CubicFace>(faceId));
          if (face_connectivity_local.isBoundaryFace(f)) continue;

          float faceCoords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const gni =
                face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
          }

          float normal[3];
          mesh_local.faceNormal(e, static_cast<model::CubicFace>(faceId), normal);

          bool const e_is_owner = (face_connectivity_local.elemOwner(f) == e);
          int const neighbor_e  = e_is_owner ? face_connectivity_local.elemNeighbor(f)
                                              : face_connectivity_local.elemOwner(f);
          model::CubicFace const neighbor_local_face = static_cast<model::CubicFace>(
              e_is_owner ? face_connectivity_local.localFaceNeighbor(f)
                         : face_connectivity_local.localFaceOwner(f));
          real_t const gamma = computeSIPGPenalty<ORDER>(faceCoords, elementCoords, penalty_local);

          INTEGRAL_TYPE::computeInterfaceFluxTerm(
              faceCoords, elementCoords, faceId,
              [&](const int i, const int j, const int k, const real_t val) {
                int const ei      = face_to_elem_dof[faceId][i];
                int const ej      = face_to_elem_dof[faceId][j];
                int const ej_perm = face_to_elem_dof[static_cast<int>(neighbor_local_face)]
                                                  [face_connectivity_local.getNeighborFaceDof(f, j)];
                int const ei_perm = face_to_elem_dof[static_cast<int>(neighbor_local_face)]
                                                  [face_connectivity_local.getNeighborFaceDof(f, i)];
                stiff_local_view(e, ei) +=
                    inv_rho * (-0.5f * val * current_field(e, ej) * normal[k] +
                                0.5f * val * current_field(neighbor_e, ej_perm) * normal[k]);
                stiff_local_view(e, ej) +=
                    inv_rho * (-0.5f * val * current_field(e, ei) * normal[k] +
                                0.5f * val * current_field(neighbor_e, ei_perm) * normal[k]);
              });

          for (int i = 0; i < knumNodesPerFace; ++i) {
            int const ei      = face_to_elem_dof[faceId][i];
            int const ei_perm = face_to_elem_dof[static_cast<int>(neighbor_local_face)]
                                              [face_connectivity_local.getNeighborFaceDof(f, i)];
            stiff_local_view(e, ei) += gamma * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
                                       (current_field(e, ei) - current_field(neighbor_e, ei_perm));
          }
        }
      });
}

//============================================================================
// applyVerlet - Kernel 3
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyVerlet(
    int kNumElem, float dt, ARRAY_REAL_VIEW current_field, ARRAY_REAL_VIEW prev_field) {
  float const dt_local  = dt;
  float const dt2_local = dt * dt;
  ARRAY_REAL_VIEW mass_local_view  = m_mass_local_;
  ARRAY_REAL_VIEW stiff_local_view = m_stiff_local_;
  ARRAY_REAL_VIEW damp_local_view  = m_damp_local_;
  ARRAY_REAL_VIEW rhs_elem_local   = m_rhs_elem_;

  Kokkos::parallel_for(
      "DG Verlet", kNumElem, KOKKOS_LAMBDA(const int e) {
        for (int i = 0; i < kPointsPerElement; ++i) {
          float const M = mass_local_view(e, i);
          float const K = stiff_local_view(e, i) + rhs_elem_local(e, i);
          float const D = damp_local_view(e, i);
          prev_field(e, i) =
              (2.0f * M * current_field(e, i) - dt2_local * K - (M - 0.5f * dt_local * D) * prev_field(e, i)) /
              (M + 0.5f * dt_local * D);
        }
      });
}

//============================================================================
// updateFields - Orchestrates the 3 kernels with fences between them
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFields(float dt,
                                                                                         const DataType& data) {
  int const kNumElem = m_mesh.getNumberOfElements();
  // SEM convention: current_field = p^n, prev_field = p^{n-1}; result written into prev_field
  ARRAY_REAL_VIEW current_field = data.getCurrentField(0);
  ARRAY_REAL_VIEW prev_field    = data.getPreviousField(0);

  computeVolumeAndBoundary(kNumElem, current_field);
  FENCE
  computeInterfaceFlux(kNumElem, current_field);
  FENCE
  applyVerlet(kNumElem, dt, current_field, prev_field);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
