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
// Update Solution Forward (Phase 2)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolutionForward(
    const float& dt, Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  updateFieldsForward(dt, myData);
  FENCE
}

//============================================================================
// Update Solution Backward (Phase 2 - Adjoint Mode)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolutionBackward(
    const float& dt, Solver::DataStruct& data) {
  throw std::runtime_error(
      "DGsolver::updateSolutionBackward not yet implemented. "
      "DG backward mode requires 3-buffer wavefield support.");
}

//============================================================================
// outputSolutionValues
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const arrayReal& fieldGlobal, const char* fieldName) {
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
  // Skipped when a caller already injected one via setFaceConnectivity().
  if (m_face_connectivity_.getNumberOfFaces() == 0) m_face_connectivity_.build(m_mesh);
  int const kNumElem = m_mesh.getNumberOfElements();
  m_rhs_elem_ = allocateArray2D<arrayReal>(kNumElem, kPointsPerElement, "rhsElem");
  m_mass_local_ = allocateArray2D<arrayReal>(kNumElem, kPointsPerElement, "massLocal");
  m_stiff_local_ = allocateArray2D<arrayReal>(kNumElem, kPointsPerElement, "stiffLocal");
  m_damp_local_ = allocateArray2D<arrayReal>(kNumElem, kPointsPerElement, "dampLocal");
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
  int const nb_rhs_element = data.getRhsElement().extent(0);
  auto rhs_elem_view = m_rhs_elem_;
  auto rhs_element_view = data.getRhsElement();
  auto rhs_term_view = data.getRhsTerm(0);
  auto rhs_weights_view = data.getRhsWeights();

  // Direct assignment (not +=): m_rhs_elem_ is zero-initialized at allocation and
  // non-source entries are never touched, so only source entries need overwriting.
  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int s) {
        int const src_elem = rhs_element_view[s];
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
    int kNumElem, arrayReal current_field, int timeSample, const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : kNumElem;

  arrayReal mass_local_view = m_mass_local_;
  arrayReal stiff_local_view = m_stiff_local_;
  arrayReal damp_local_view = m_damp_local_;
  arrayReal rhs_elem_view = m_rhs_elem_;

  // RHS source data, folded into this kernel to avoid a separate launch.
  int const nb_rhs_element = data.getRhsElement().extent(0);
  auto rhs_element_view = data.getRhsElement();
  auto rhs_term_view = data.getRhsTerm(0);
  auto rhs_weights_view = data.getRhsWeights();

  Kokkos::parallel_for(
      "DG Volume+Boundary", n_iter, KOKKOS_LAMBDA(const int _loop_idx) {
        int const e = list_on ? list_local[_loop_idx] : _loop_idx;
        float massLocal[kPointsPerElement] = {0};
        float stiffLocal[kPointsPerElement] = {0};
        float elementCoords[8][3];
        auto const eIdx = mesh_local.elementIndex(e);
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv)
              mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                      elementCoords[iv + 2 * jv + 4 * kv]);

        real_t const vp = mesh_local.getModelVpOnElement(e);
        real_t const rho = mesh_local.getModelRhoOnElement(e);
        real_t const inv_model_factor = 1.0f / (vp * vp * rho);
        real_t const inv_rho = 1.0f / rho;

        INTEGRAL_TYPE::computeMassTerm(elementCoords,
                                       [&](const int j, const real_t val) { massLocal[j] += inv_model_factor * val; });

        real_t p_local[kPointsPerElement];
        for (int i = 0; i < kPointsPerElement; ++i) p_local[i] = current_field(e, i);
        INTEGRAL_TYPE::computeStiffnessTermSumFact(elementCoords, p_local, stiffLocal,
                                                   [&](const int, const int, const int) -> real_t { return inv_rho; });

        // Fold in the external RHS forcing: write m_rhs_elem_ for DG source
        // elements (same values applyRHSTerm would produce), so the separate
        // applyRHSTerm launch can be skipped in the DG-SEM path. timeSample < 0
        // means "do not fold" (standalone DG path fills m_rhs_elem_ separately).
        if (timeSample >= 0) {
          for (int s = 0; s < nb_rhs_element; ++s) {
            if (rhs_element_view[s] == e) {
              float const wavelet_val = rhs_term_view(s, timeSample);
              for (int i = 0; i < kPointsPerElement; ++i)
                rhs_elem_view(e, i) = -wavelet_val * rhs_weights_view(s, i);
            }
          }
        }

        for (int i = 0; i < kPointsPerElement; ++i) {
          mass_local_view(e, i) = massLocal[i];
          stiff_local_view(e, i) = stiffLocal[i];
          damp_local_view(e, i) = 0.0f;  // zeroed here; filled by computeBoundaryDampingAndInterfaceFlux
        }
      });
}

//============================================================================
// computeBoundaryDampingAndInterfaceFlux - Kernel 1b+2, fused (face-loop: boundary faces take
// the damping branch, interior faces take the SIPG flux branch; mutually exclusive per face,
// disjoint accumulators)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeBoundaryDampingAndInterfaceFlux(
    int kNumFaces, arrayReal current_field, const vectorReal* p_sem, const vectorReal* work_sem,
    const vectorInt* element_type) {
  auto mesh_local = m_mesh;

  bool const list_on = m_list_mode_;
  auto list_local = m_face_list_;
  int const n_iter = list_on ? m_n_face_list_ : kNumFaces;

  auto face_connectivity_local = m_face_connectivity_;
  auto const face_to_elem_dof = kFaceToElemDof;  // local copy for device capture
  auto const face_to_elem_dof_depth = kFaceToElemDofAtDepth;
  arrayReal damp_local_view = m_damp_local_;
  arrayReal stiff_local_view = m_stiff_local_;
  real_t const penalty_local = m_penalty_factor_;

  // DG-SEM coupling state (null when running standalone DG).
  bool const has_coupling = (p_sem != nullptr) && (work_sem != nullptr) && (element_type != nullptr);
  auto p_sem_view = has_coupling ? *p_sem : vectorReal();
  auto work_sem_view = has_coupling ? *work_sem : vectorReal();
  auto element_type_view = has_coupling ? *element_type : vectorInt();

  Kokkos::parallel_for(
      "DG BoundaryDamping+InterfaceFlux", n_iter, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = list_on ? list_local[_loop_idx] : _loop_idx;

        if (face_connectivity_local.isBoundaryFace(f)) {
          int const e = face_connectivity_local.elemOwner(f);
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
          return;
        }

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        // DG-SEM interface face: owner and neighbor belong to different domains.
        // Take the SIPG coupling branch (SEM pressure ↔ DG stiff_local_ and work_sem).
        // kElementTypeDG == 1 (defined in dg-sem_solver.h); the DG solver uses the literal
        // to avoid a cross-header dependency.
        if (has_coupling && element_type_view(owner_e) != element_type_view(neighbor_e)) {
          bool const owner_is_dg = (element_type_view(owner_e) == 1);
          int const dg_e = owner_is_dg ? owner_e : neighbor_e;
          int const sem_e = owner_is_dg ? neighbor_e : owner_e;
          int const fid_dg = owner_is_dg ? fid_o : fid_n;
          int const fid_sem = owner_is_dg ? fid_n : fid_o;

          auto dg_to_sem = [&](int i) {
            return owner_is_dg ? face_connectivity_local.getNeighborFaceDof(f, i)
                               : face_connectivity_local.getOwnerFaceDof(f, i);
          };
          auto sem_to_dg = [&](int i) {
            return owner_is_dg ? face_connectivity_local.getOwnerFaceDof(f, i)
                               : face_connectivity_local.getNeighborFaceDof(f, i);
          };

          float faceCoords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const gni =
                face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
          }

          float dg_coords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(dg_e);
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                          dg_coords[iv + 2 * jv + 4 * kv]);
          }

          float sem_coords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(sem_e);
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                          sem_coords[iv + 2 * jv + 4 * kv]);
          }

          float normal_dg[3];
          mesh_local.faceNormal(dg_e, static_cast<model::CubicFace>(fid_dg), normal_dg);

          real_t const face_area = computeFaceArea(faceCoords);
          real_t const inv_rho_dg = 1.0f / mesh_local.getModelRhoOnElement(dg_e);
          real_t const gamma_dg = computeSIPGPenaltyFromArea<ORDER>(face_area, dg_coords, penalty_local);
          real_t const inv_rho_sem = 1.0f / mesh_local.getModelRhoOnElement(sem_e);
          real_t const gamma_sem = computeSIPGPenaltyFromArea<ORDER>(face_area, sem_coords, penalty_local);

          float stiff_dg_local[knumNodesPerFace] = {0};
          float work_sem_local[knumNodesPerFace] = {0};

          INTEGRAL_TYPE::computeInterfaceFluxTerm(
              faceCoords, dg_coords, fid_dg, [&](const int i, const int j, const int k, const real_t val) {
                int const ei = face_to_elem_dof[fid_dg][i];
                int const ej = face_to_elem_dof[fid_dg][j];
                int const sd_j = dg_to_sem(j);
                int const gn_j = face_connectivity_local.getGlobalNodeFromFace(f, sd_j);
                float const nk = normal_dg[k];
                stiff_dg_local[i] += inv_rho_dg * nk * (-0.5f * val * current_field(dg_e, ej) + 0.5f * val * p_sem_view(gn_j));
                stiff_dg_local[j] += inv_rho_dg * nk * (-0.5f * val * current_field(dg_e, ei));
                work_sem_local[sd_j] += inv_rho_dg * nk * (0.5f * val * current_field(dg_e, ei));
              });

          for (int i = 0; i < knumNodesPerFace; ++i) {
            int const ei = face_to_elem_dof[fid_dg][i];
            int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, dg_to_sem(i));
            stiff_dg_local[i] +=
                gamma_dg * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) * (current_field(dg_e, ei) - p_sem_view(gn_i));
          }

          INTEGRAL_TYPE::computeInterfaceFluxTerm(
              faceCoords, sem_coords, fid_sem, [&](const int i, const int j, const int k, const real_t val) {
                int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, i);
                int const gn_j = face_connectivity_local.getGlobalNodeFromFace(f, j);
                int const sd_j = sem_to_dg(j);
                int const ej_perm = face_to_elem_dof[fid_dg][sd_j];
                float const nk = -normal_dg[k];  // SEM outward = -DG outward
                stiff_dg_local[sd_j] += inv_rho_sem * nk * (0.5f * val * p_sem_view(gn_i));
                work_sem_local[i] += inv_rho_sem * nk * (-0.5f * val * p_sem_view(gn_j) + 0.5f * val * current_field(dg_e, ej_perm));
                work_sem_local[j] += inv_rho_sem * nk * (-0.5f * val * p_sem_view(gn_i));
              });

          for (int i = 0; i < knumNodesPerFace; ++i) {
            int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, i);
            int const ei_perm = face_to_elem_dof[fid_dg][sem_to_dg(i)];
            work_sem_local[i] += inv_rho_sem * gamma_sem * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
                                 (p_sem_view(gn_i) - current_field(dg_e, ei_perm));
          }

          for (int i = 0; i < knumNodesPerFace; ++i) {
            ATOMICADD(stiff_local_view(dg_e, face_to_elem_dof[fid_dg][i]), stiff_dg_local[i]);
            ATOMICADD(work_sem_view(face_connectivity_local.getGlobalNodeFromFace(f, i)), work_sem_local[i]);
          }
          return;
        }

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni = face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float owner_coords[8][3];
        auto const eIdx_o = mesh_local.elementIndex(owner_e);
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv)
              mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx_o, iv, jv, kv),
                                      owner_coords[iv + 2 * jv + 4 * kv]);

        float neighbor_coords[8][3];
        auto const eIdx_n = mesh_local.elementIndex(neighbor_e);
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv)
              mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx_n, iv, jv, kv),
                                      neighbor_coords[iv + 2 * jv + 4 * kv]);

        real_t const inv_rho_o = 1.0f / mesh_local.getModelRhoOnElement(owner_e);
        real_t const inv_rho_n = 1.0f / mesh_local.getModelRhoOnElement(neighbor_e);

        float normal[3];
        mesh_local.faceNormal(owner_e, static_cast<model::CubicFace>(fid_o), normal);

        real_t const face_area = computeFaceArea(faceCoords);
        real_t const gamma_o = computeSIPGPenaltyFromArea<ORDER>(face_area, owner_coords, penalty_local);
        real_t const gamma_n = computeSIPGPenaltyFromArea<ORDER>(face_area, neighbor_coords, penalty_local);

        float stiff_o[kPointsPerElement] = {0};
        float stiff_n[kPointsPerElement] = {0};

        // --- Owner side (outward normal = normal[]) ---
        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, owner_coords, fid_o,
            [&](const int i, const int j, const int k, const real_t val) {
              int const ei = face_to_elem_dof[fid_o][i];
              int const ej = face_to_elem_dof[fid_o][j];
              int const ej_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, j)];
              float const nk = normal[k];
              stiff_o[ei] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) * nk +
                                          0.5f * val * current_field(neighbor_e, ej_perm) * nk);
              stiff_o[ej] += inv_rho_o * (-0.5f * val * current_field(owner_e, ei) * nk);
              stiff_n[ej_perm] += inv_rho_o * (0.5f * val * current_field(owner_e, ei) * nk);
            },
            [&](const int m, const int j, const int k, const real_t val) {
              int const em = face_to_elem_dof_depth[fid_o][j][m];
              int const ej = face_to_elem_dof[fid_o][j];
              int const ej_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, j)];
              float const nk = normal[k];
              stiff_o[em] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) * nk +
                                          0.5f * val * current_field(neighbor_e, ej_perm) * nk);
              stiff_o[ej] += inv_rho_o * (-0.5f * val * current_field(owner_e, em) * nk);
              stiff_n[ej_perm] += inv_rho_o * (0.5f * val * current_field(owner_e, em) * nk);
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei = face_to_elem_dof[fid_o][i];
          int const ei_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, i)];
          stiff_o[ei] += gamma_o * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
                         (current_field(owner_e, ei) - current_field(neighbor_e, ei_perm));
        }

        // --- Neighbor side (outward normal = -normal[]) ---
        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, neighbor_coords, fid_n,
            [&](const int i, const int j, const int k, const real_t val) {
              int const ei = face_to_elem_dof[fid_n][i];
              int const ej = face_to_elem_dof[fid_n][j];
              int const ej_perm = face_to_elem_dof[fid_o][face_connectivity_local.getOwnerFaceDof(f, j)];
              float const nk = -normal[k];
              stiff_n[ei] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) * nk +
                                          0.5f * val * current_field(owner_e, ej_perm) * nk);
              stiff_n[ej] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ei) * nk);
              stiff_o[ej_perm] += inv_rho_n * (0.5f * val * current_field(neighbor_e, ei) * nk);
            },
            [&](const int m, const int j, const int k, const real_t val) {
              int const em = face_to_elem_dof_depth[fid_n][j][m];
              int const ej = face_to_elem_dof[fid_n][j];
              int const ej_perm = face_to_elem_dof[fid_o][face_connectivity_local.getOwnerFaceDof(f, j)];
              float const nk = -normal[k];
              stiff_n[em] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) * nk +
                                          0.5f * val * current_field(owner_e, ej_perm) * nk);
              stiff_n[ej] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, em) * nk);
              stiff_o[ej_perm] += inv_rho_n * (0.5f * val * current_field(neighbor_e, em) * nk);
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei = face_to_elem_dof[fid_n][i];
          int const ei_perm = face_to_elem_dof[fid_o][face_connectivity_local.getOwnerFaceDof(f, i)];
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
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyVerlet(int kNumElem, float dt,
                                                                                        arrayReal current_field,
                                                                                        arrayReal prev_field) {
  float const dt_local = dt;
  float const dt2_local = dt * dt;

  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : kNumElem;

  arrayReal mass_local_view = m_mass_local_;
  arrayReal stiff_local_view = m_stiff_local_;
  arrayReal damp_local_view = m_damp_local_;
  arrayReal rhs_elem_local = m_rhs_elem_;

  Kokkos::parallel_for(
      "DG Verlet", n_iter, KOKKOS_LAMBDA(const int _loop_idx) {
        int const e = list_on ? list_local[_loop_idx] : _loop_idx;
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
//============================================================================
// updateFieldsForward - Orchestrates the 3 kernels with fences between them
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsForward(float dt,
                                                                                                const DataType& data) {
  int const kNumElem = m_mesh.getNumberOfElements();
  int const kNumFaces = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  // SEM convention: current_field = p^n, prev_field = p^{n-1}; result written into prev_field
  arrayReal current_field = data.getCurrentField(0);
  arrayReal prev_field = data.getPreviousField(0);

  computeVolumeAndBoundary(kNumElem, current_field, -1, data);
  FENCE
  computeBoundaryDampingAndInterfaceFlux(kNumFaces, current_field);
  FENCE
  applyVerlet(kNumElem, dt, current_field, prev_field);
}

//============================================================================
// updateFieldsBackward - Backward/adjoint mode (not yet fully implemented for DG)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsBackward(float dt,
                                                                                                 const DataType& data) {
  throw std::runtime_error(
      "DGsolver::updateFieldsBackward not yet implemented. "
      "DG backward mode requires 3-buffer wavefield support.");
}

//============================================================================
// updateFieldsFromListForward - Verlet update restricted to a compact element list (forward mode)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForward(
    float dt, const DataType& data, const vectorInt& elem_list, int n_elems) {
  m_elem_list_ = elem_list;
  m_n_elem_list_ = n_elems;
  faceListFromElementList();
  m_list_mode_ = true;
  updateFieldsForward(dt, data);
  m_list_mode_ = false;
}

//============================================================================
// updateFieldsFromListBackward - Verlet update restricted to a compact element list (backward mode)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListBackward(
    float dt, const DataType& data, const vectorInt& elem_list, int n_elems) {
  throw std::runtime_error(
      "DGsolver::updateFieldsFromListBackward not yet implemented. "
      "DG backward mode requires 3-buffer wavefield support.");
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
