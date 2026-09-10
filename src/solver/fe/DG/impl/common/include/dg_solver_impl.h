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
    int kNumElem, arrayReal current_field) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : kNumElem;

  arrayReal mass_local_view = m_mass_local_;
  arrayReal stiff_local_view = m_stiff_local_;
  arrayReal damp_local_view = m_damp_local_;

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
    int kNumFaces, arrayReal current_field) {
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

        // Face-sized accumulators, indexed by face dof in each side's own face numbering: the
        // interface flux only ever touches the shared face's (ORDER+1)^2 dofs, so element-sized
        // arrays cost a 7x larger per-thread footprint (which spills out of registers) and flush an
        // atomic zero to every off-face dof. Same reasoning as the DG-SEM coupling kernel.
        float stiff_o[knumNodesPerFace] = {0};
        float stiff_n[knumNodesPerFace] = {0};

        // Driving one quadrature point at a time is what keeps the accumulators face-sized: over a
        // whole face the normal channel spans every element dof, but for a single point it only
        // reaches the ORDER+1 dofs of the line through it.
        for (int q = 0; q < knumNodesPerFace; ++q) {
          // Face-normal accumulators. These dofs are off-face, so they cannot go through the
          // face-sized rows. The lines of two quadrature points are disjoint, hence plain adds.
          float norm_o[ORDER + 1] = {0};
          float norm_n[ORDER + 1] = {0};

          // --- Owner side (outward normal = normal[]) ---
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, owner_coords, fid_o,
              [&](const int i, const int j, const int k, const real_t val) {
                int const nfd_j = face_connectivity_local.getNeighborFaceDof(f, j);
                int const ei = face_to_elem_dof[fid_o][i];
                int const ej = face_to_elem_dof[fid_o][j];
                int const ej_perm = face_to_elem_dof[fid_n][nfd_j];
                float const nk = normal[k];
                stiff_o[i] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) * nk +
                                           0.5f * val * current_field(neighbor_e, ej_perm) * nk);
                stiff_o[j] += inv_rho_o * (-0.5f * val * current_field(owner_e, ei) * nk);
                stiff_n[nfd_j] += inv_rho_o * (0.5f * val * current_field(owner_e, ei) * nk);
              },
              [&](const int m, const int j, const int k, const real_t val) {
                int const nfd_j = face_connectivity_local.getNeighborFaceDof(f, j);
                int const em = face_to_elem_dof_depth[fid_o][j][m];
                int const ej = face_to_elem_dof[fid_o][j];
                int const ej_perm = face_to_elem_dof[fid_n][nfd_j];
                float const nk = normal[k];
                norm_o[m] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) * nk +
                                          0.5f * val * current_field(neighbor_e, ej_perm) * nk);
                stiff_o[j] += inv_rho_o * (-0.5f * val * current_field(owner_e, em) * nk);
                stiff_n[nfd_j] += inv_rho_o * (0.5f * val * current_field(owner_e, em) * nk);
              });

          // --- Neighbor side (outward normal = -normal[]) ---
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, neighbor_coords, fid_n,
              [&](const int i, const int j, const int k, const real_t val) {
                int const ofd_j = face_connectivity_local.getOwnerFaceDof(f, j);
                int const ei = face_to_elem_dof[fid_n][i];
                int const ej = face_to_elem_dof[fid_n][j];
                int const ej_perm = face_to_elem_dof[fid_o][ofd_j];
                float const nk = -normal[k];
                stiff_n[i] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) * nk +
                                           0.5f * val * current_field(owner_e, ej_perm) * nk);
                stiff_n[j] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ei) * nk);
                stiff_o[ofd_j] += inv_rho_n * (0.5f * val * current_field(neighbor_e, ei) * nk);
              },
              [&](const int m, const int j, const int k, const real_t val) {
                int const ofd_j = face_connectivity_local.getOwnerFaceDof(f, j);
                int const em = face_to_elem_dof_depth[fid_n][j][m];
                int const ej = face_to_elem_dof[fid_n][j];
                int const ej_perm = face_to_elem_dof[fid_o][ofd_j];
                float const nk = -normal[k];
                norm_n[m] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) * nk +
                                          0.5f * val * current_field(owner_e, ej_perm) * nk);
                stiff_n[j] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, em) * nk);
                stiff_o[ofd_j] += inv_rho_n * (0.5f * val * current_field(neighbor_e, em) * nk);
              });

          // The normal lines are off-face, so they bypass the face-sized flush below. Both
          // callbacks always fire with j == q, so the line is the one through face dof q on each
          // side, in that side's own face numbering.
          for (int m = 0; m <= ORDER; ++m) {
            ATOMICADD(stiff_local_view(owner_e, face_to_elem_dof_depth[fid_o][q][m]), norm_o[m]);
            ATOMICADD(stiff_local_view(neighbor_e, face_to_elem_dof_depth[fid_n][q][m]), norm_n[m]);
          }
        }

        // SIPG penalty and atomic write-back, fused: both sides use the same damping term at face
        // dof i, so it is computed once here instead of once per side. Atomic because several faces
        // share an element.
        for (int i = 0; i < knumNodesPerFace; ++i) {
          real_t const damping_i = INTEGRAL_TYPE::computeDampingTerm(i, faceCoords);

          int const ei_o = face_to_elem_dof[fid_o][i];
          int const ei_o_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, i)];
          stiff_o[i] += gamma_o * damping_i * (current_field(owner_e, ei_o) - current_field(neighbor_e, ei_o_perm));

          int const ei_n = face_to_elem_dof[fid_n][i];
          int const ei_n_perm = face_to_elem_dof[fid_o][face_connectivity_local.getOwnerFaceDof(f, i)];
          stiff_n[i] += gamma_n * damping_i * (current_field(neighbor_e, ei_n) - current_field(owner_e, ei_n_perm));

          ATOMICADD(stiff_local_view(owner_e, ei_o), stiff_o[i]);
          ATOMICADD(stiff_local_view(neighbor_e, ei_n), stiff_n[i]);
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

  computeVolumeAndBoundary(kNumElem, current_field);
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
