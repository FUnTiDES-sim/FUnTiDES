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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolutionForward(
    const float& dt, Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  updateFieldsForward(dt, myData);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolutionBackward(
    const float& dt, Solver::DataStruct& data) {
  throw std::runtime_error(
      "DGsolver::updateSolutionBackward not yet implemented. "
      "DG backward mode requires 3-buffer wavefield support.");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const arrayReal& fieldGlobal, const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << fieldGlobal(e, 0) << endl;
}

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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeForces(const float& dt,
                                                                                          const int& timeSample,
                                                                                          Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  applyRHSTerm(timeSample, dt, myData);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyRHSTerm(int timeSample, float /*dt*/,
                                                                                         const DataType& data) {
  int const nb_rhs_element = data.getRhsElement().extent(0);
  auto rhs_elem_view = m_rhs_elem_;
  auto rhs_element_view = data.getRhsElement();
  auto rhs_term_view = data.getRhsTerm(0);
  auto rhs_weights_view = data.getRhsWeights();

  // Assignment, not +=: m_rhs_elem_ is zero-initialized at allocation and non-source
  // entries are never touched, so only source entries need overwriting.
  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int s) {
        int const src_elem = rhs_element_view[s];
        float const wavelet_val = rhs_term_view(s, timeSample);
        for (int dof = 0; dof < kPointsPerElement; ++dof) {
          rhs_elem_view(src_elem, dof) = -wavelet_val * rhs_weights_view(s, dof);
        }
      });
}

/// Kernel 1: per-element mass and stiffness (volume term); resets the damping accumulator.
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
          damp_local_view(e, i) = 0.0f;  // accumulated later by computeBoundaryDampingAndInterfaceFlux
        }
      });
}

/// Kernels 1b and 2, fused: loop over faces. Boundary faces add absorbing damping, interior faces
/// add the SIPG interface flux. The two branches are exclusive per face and write disjoint
/// accumulators (damping vs stiffness).
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
        float const neg_normal[3] = {-normal[0], -normal[1], -normal[2]};

        real_t const face_area = computeFaceArea(faceCoords);
        real_t const gamma_o = computeSIPGPenaltyFromArea<ORDER>(face_area, owner_coords, penalty_local);
        real_t const gamma_n = computeSIPGPenaltyFromArea<ORDER>(face_area, neighbor_coords, penalty_local);

        // Face-sized accumulators, indexed by face dof in each side's own face numbering. They are
        // indexed dynamically, so they live in local memory: element-sized arrays would put
        // 2*(ORDER+1)^3 floats per thread there (2.7 kB at order 6), and every read-modify-write in
        // the quadrature loop goes through L1. Face-sized arrays cut that by (ORDER+1).
        float stiff_o[knumNodesPerFace] = {0};
        float stiff_n[knumNodesPerFace] = {0};

        // One quadrature point at a time keeps the accumulators face-sized: over a whole face the
        // normal channel spans every element dof, but for a single point it only reaches the
        // ORDER+1 dofs of the line through it, which fit in registers.
        for (int q = 0; q < knumNodesPerFace; ++q) {
          // Face-normal accumulators for the off-face dofs, so they bypass the face-sized arrays.
          // The lines of two quadrature points are disjoint, hence plain adds.
          float norm_o[ORDER + 1] = {0};
          float norm_n[ORDER + 1] = {0};

          // Both callbacks always fire with j == q, so everything derived from j is invariant over
          // the point and hoisted here. This includes the two face-dof permutations, which are
          // virtual calls.
          int const nfd_q = face_connectivity_local.getNeighborFaceDof(f, q);
          int const ofd_q = face_connectivity_local.getOwnerFaceDof(f, q);
          int const ej_o = face_to_elem_dof[fid_o][q];
          int const ej_o_perm = face_to_elem_dof[fid_n][nfd_q];
          int const ej_n = face_to_elem_dof[fid_n][q];
          int const ej_n_perm = face_to_elem_dof[fid_o][ofd_q];

          real_t const half_o = 0.5f * inv_rho_o;
          real_t const half_n = 0.5f * inv_rho_n;

          // Pressure jump across the face at q, seen from each side. The first contribution of
          // both callbacks is exactly val times this.
          real_t const dp_o = half_o * (current_field(neighbor_e, ej_o_perm) - current_field(owner_e, ej_o));
          real_t const dp_n = half_n * (current_field(owner_e, ej_n_perm) - current_field(neighbor_e, ej_n));

          // The two other contributions of each side land on fixed slots (j and its image on the
          // opposite side) with the same magnitude and opposite signs. One register carries both
          // and is flushed once below, which avoids many dynamically indexed read-modify-writes.
          float acc_o = 0.0f;
          float acc_n = 0.0f;

          // Owner side, outward normal = normal[]. The callbacks are normal-contracted (the
          // discretization folds sum_k C_ijk * n_k), so each contribution fires once instead of
          // once per physical direction.
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, owner_coords, fid_o, normal,
              [&](const int i, const int, const real_t val) {
                stiff_o[i] += val * dp_o;
                acc_o -= half_o * val * current_field(owner_e, face_to_elem_dof[fid_o][i]);
              },
              [&](const int m, const int, const real_t val) {
                norm_o[m] += val * dp_o;
                acc_o -= half_o * val * current_field(owner_e, face_to_elem_dof_depth[fid_o][q][m]);
              });

          // Neighbor side, outward normal = -normal[].
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, neighbor_coords, fid_n, neg_normal,
              [&](const int i, const int, const real_t val) {
                stiff_n[i] += val * dp_n;
                acc_n -= half_n * val * current_field(neighbor_e, face_to_elem_dof[fid_n][i]);
              },
              [&](const int m, const int, const real_t val) {
                norm_n[m] += val * dp_n;
                acc_n -= half_n * val * current_field(neighbor_e, face_to_elem_dof_depth[fid_n][q][m]);
              });

          // Mirror each register onto the opposite side (negation is exact in IEEE-754, so this
          // equals what the callbacks would have accumulated there).
          stiff_o[q] += acc_o;
          stiff_n[nfd_q] -= acc_o;
          stiff_n[q] += acc_n;
          stiff_o[ofd_q] -= acc_n;

          // Flush the off-face normal lines: with j == q, the line is the one through face dof q
          // on each side, in that side's own face numbering.
          for (int m = 0; m <= ORDER; ++m) {
            ATOMICADD(stiff_local_view(owner_e, face_to_elem_dof_depth[fid_o][q][m]), norm_o[m]);
            ATOMICADD(stiff_local_view(neighbor_e, face_to_elem_dof_depth[fid_n][q][m]), norm_n[m]);
          }
        }

        // SIPG penalty and write-back, fused: both sides use the same damping term at face dof i,
        // so it is computed once. Atomic because several faces share an element.
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

/// Kernel 3: explicit Verlet step with damping. Reads current_field (step n) and overwrites
/// prev_field (step n-1) with the field at step n+1.
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

/// Runs the three kernels in order with a fence between them: volume terms, face terms, Verlet.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsForward(float dt,
                                                                                                const DataType& data) {
  int const kNumElem = m_mesh.getNumberOfElements();
  int const kNumFaces = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  // current_field = p^n, prev_field = p^{n-1}; the result is written into prev_field.
  arrayReal current_field = data.getCurrentField(0);
  arrayReal prev_field = data.getPreviousField(0);

  computeVolumeAndBoundary(kNumElem, current_field);
  FENCE
  computeBoundaryDampingAndInterfaceFlux(kNumFaces, current_field);
  FENCE
  applyVerlet(kNumElem, dt, current_field, prev_field);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsBackward(float dt,
                                                                                                 const DataType& data) {
  throw std::runtime_error(
      "DGsolver::updateFieldsBackward not yet implemented. "
      "DG backward mode requires 3-buffer wavefield support.");
}

/// Same as updateFieldsForward, restricted to the first n_elems elements of elem_list (and the
/// faces derived from them).
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListBackward(
    float dt, const DataType& data, const vectorInt& elem_list, int n_elems) {
  throw std::runtime_error(
      "DGsolver::updateFieldsFromListBackward not yet implemented. "
      "DG backward mode requires 3-buffer wavefield support.");
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_SOLVER_IMPL_H_
