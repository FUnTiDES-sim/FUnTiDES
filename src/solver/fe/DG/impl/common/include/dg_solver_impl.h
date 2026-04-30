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
  ARRAY_REAL_VIEW mass_local_view  = m_mass_local_;
  ARRAY_REAL_VIEW stiff_local_view = m_stiff_local_;
  ARRAY_REAL_VIEW damp_local_view  = m_damp_local_;

  Kokkos::parallel_for(
      "DG Volume+Boundary", kNumElem, KOKKOS_LAMBDA(const int e) {
        float massLocal[kPointsPerElement]  = {0};
        float stiffLocal[kPointsPerElement] = {0};
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

        INTEGRAL_TYPE::computeMassTerm(
            elementCoords, [&](const int j, const real_t val) { massLocal[j] += inv_model_factor * val; });

        real_t p_local[kPointsPerElement];
        for (int i = 0; i < kPointsPerElement; ++i) p_local[i] = current_field(e, i);
        INTEGRAL_TYPE::computeStiffnessTermSumFact(
            elementCoords, p_local, stiffLocal,
            [&](const int , const int , const int ) -> real_t { return inv_rho; });

        for (int i = 0; i < kPointsPerElement; ++i) {
          mass_local_view(e, i)  = massLocal[i];
          stiff_local_view(e, i) = stiffLocal[i];
          damp_local_view(e, i)  = 0.0f;  // zeroed here; filled by computeBoundaryDamping
        }
      });
}

//============================================================================
// computeBoundaryDamping - Kernel 1b (face-loop, boundary faces only)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeBoundaryDamping(
    int kNumFaces) {
  auto mesh_local              = m_mesh;
  auto face_connectivity_local = m_face_connectivity_;
  auto const face_to_elem_dof  = kFaceToElemDof;  // local copy for device capture
  ARRAY_REAL_VIEW damp_local_view = m_damp_local_;

  Kokkos::parallel_for(
      "DG Boundary Damping", kNumFaces, KOKKOS_LAMBDA(const int f) {
        if (!face_connectivity_local.isBoundaryFace(f)) return;

        int const e      = face_connectivity_local.elemOwner(f);
        int const faceId = face_connectivity_local.localFaceOwner(f);

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni =
              face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        real_t const inv_vp = 1.0f / mesh_local.getModelVpOnElement(e);

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei = face_to_elem_dof[faceId][i];
          ATOMICADD(damp_local_view(e, ei), inv_vp * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords));
        }
      });
}

//============================================================================
// computeInterfaceFlux - Kernel 2 (face-loop: each interior face processed once)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeInterfaceFlux(
    int /*kNumElem*/, ARRAY_REAL_VIEW current_field) {
  auto mesh_local              = m_mesh;
  auto face_connectivity_local = m_face_connectivity_;
  auto const face_to_elem_dof  = kFaceToElemDof;  // local copy for device capture
  ARRAY_REAL_VIEW stiff_local_view = m_stiff_local_;
  real_t const penalty_local       = m_penalty_factor_;

  int const kNumFaces = face_connectivity_local.getNumberOfFaces();

  Kokkos::parallel_for(
      "DG Interface Flux", kNumFaces, KOKKOS_LAMBDA(const int f) {
        if (face_connectivity_local.isBoundaryFace(f)) return;

        int const owner_e    = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o      = face_connectivity_local.localFaceOwner(f);
        int const fid_n      = face_connectivity_local.localFaceNeighbor(f);

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni =
              face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float owner_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(owner_e);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), owner_coords[I++]);
        }

        float neighbor_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(neighbor_e);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), neighbor_coords[I++]);
        }

        real_t const inv_rho_o = 1.0f / mesh_local.getModelRhoOnElement(owner_e);
        real_t const inv_rho_n = 1.0f / mesh_local.getModelRhoOnElement(neighbor_e);

        float normal[3];
        mesh_local.faceNormal(owner_e, static_cast<model::CubicFace>(fid_o), normal);

        real_t const gamma_o = computeSIPGPenalty<ORDER>(faceCoords, owner_coords, penalty_local);
        real_t const gamma_n = computeSIPGPenalty<ORDER>(faceCoords, neighbor_coords, penalty_local);

        float stiff_o[kPointsPerElement] = {0};
        float stiff_n[kPointsPerElement] = {0};

        // --- Owner side (outward normal = normal[]) ---
        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, owner_coords, fid_o,
            [&](const int i, const int j, const int k, const real_t val) {
              int const ei      = face_to_elem_dof[fid_o][i];
              int const ej      = face_to_elem_dof[fid_o][j];
              int const ej_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, j)];
              int const ei_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, i)];
              float const nk    = normal[k];
              stiff_o[ei] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) * nk +
                                           0.5f * val * current_field(neighbor_e, ej_perm) * nk);
              stiff_o[ej] += inv_rho_o * (-0.5f * val * current_field(owner_e, ei) * nk +
                                           0.5f * val * current_field(neighbor_e, ei_perm) * nk);
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei      = face_to_elem_dof[fid_o][i];
          int const ei_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, i)];
          stiff_o[ei] += gamma_o * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
                         (current_field(owner_e, ei) - current_field(neighbor_e, ei_perm));
        }

        // --- Neighbor side (outward normal = -normal[]) ---
        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, neighbor_coords, fid_n,
            [&](const int i, const int j, const int k, const real_t val) {
              int const ei      = face_to_elem_dof[fid_n][i];
              int const ej      = face_to_elem_dof[fid_n][j];
              int const ej_perm = face_to_elem_dof[fid_o][face_connectivity_local.getNeighborFaceDof(f, j)];
              int const ei_perm = face_to_elem_dof[fid_o][face_connectivity_local.getNeighborFaceDof(f, i)];
              float const nk    = -normal[k];
              stiff_n[ei] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) * nk +
                                           0.5f * val * current_field(owner_e, ej_perm) * nk);
              stiff_n[ej] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ei) * nk +
                                           0.5f * val * current_field(owner_e, ei_perm) * nk);
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei      = face_to_elem_dof[fid_n][i];
          int const ei_perm = face_to_elem_dof[fid_o][face_connectivity_local.getNeighborFaceDof(f, i)];
          stiff_n[ei] += gamma_n * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
                         (current_field(neighbor_e, ei) - current_field(owner_e, ei_perm));
        }

        // Atomic write-back: multiple faces can share the same element
        for (int i = 0; i < kPointsPerElement; ++i) {
          ATOMICADD(stiff_local_view(owner_e, i), stiff_o[i]);
          ATOMICADD(stiff_local_view(neighbor_e, i), stiff_n[i]);
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
  int const kNumElem  = m_mesh.getNumberOfElements();
  int const kNumFaces = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  // SEM convention: current_field = p^n, prev_field = p^{n-1}; result written into prev_field
  ARRAY_REAL_VIEW current_field = data.getCurrentField(0);
  ARRAY_REAL_VIEW prev_field    = data.getPreviousField(0);

  computeVolumeAndBoundary(kNumElem, current_field);
  FENCE
  computeBoundaryDamping(kNumFaces);
  FENCE
  computeInterfaceFlux(kNumElem, current_field);
  FENCE
  applyVerlet(kNumElem, dt, current_field, prev_field);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
