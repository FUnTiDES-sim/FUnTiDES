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
  // A hard if/else, so the branch not taken is never instantiated: nvcc arbitrates registers per
  // compilation unit, and an unused second kernel there still costs occupancy to the live one.
  if constexpr (ORDER <= kMaxOrderForFlatFace)
    computeBoundaryDampingAndInterfaceFlux_Flat(kNumFaces, current_field);
  else
    computeBoundaryDampingAndInterfaceFlux_Team(kNumFaces, current_field);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeBoundaryDampingAndInterfaceFlux_Flat(
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

        // Face-sized accumulators, indexed by face dof in each side's own face numbering. These are
        // indexed dynamically, so they never live in registers: element-sized arrays put
        // 2*(ORDER+1)^3 floats per thread in local memory (2.7 kB at order 6) and every read-modify-
        // write in the quadrature loop goes through L1, which ncu measures as the kernel's busiest
        // pipe. Face-sized cuts that by (ORDER+1).
        float stiff_o[knumNodesPerFace] = {0};
        float stiff_n[knumNodesPerFace] = {0};

        // Driving one quadrature point at a time is what keeps the accumulators face-sized: over a
        // whole face the normal channel spans every element dof, but for a single point it only
        // reaches the ORDER+1 dofs of the line through it, which fit in registers.
        for (int q = 0; q < knumNodesPerFace; ++q) {
          // Face-normal accumulators. Off-face dofs, so they bypass the face-sized rows; the lines
          // of two quadrature points are disjoint, hence plain adds.
          float norm_o[ORDER + 1] = {0};
          float norm_n[ORDER + 1] = {0};

          // --- Owner side (outward normal = normal[]) ---
          // Normal-contracted callbacks: the discretization folds sum_k C_ijk * n_k, so each
          // contribution fires once instead of once per physical direction, dividing the local-memory
          // traffic on stiff_o/stiff_n by three.
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, owner_coords, fid_o, normal,
              [&](const int i, const int j, const real_t val) {
                int const nfd_j = face_connectivity_local.getNeighborFaceDof(f, j);
                int const ei = face_to_elem_dof[fid_o][i];
                int const ej = face_to_elem_dof[fid_o][j];
                int const ej_perm = face_to_elem_dof[fid_n][nfd_j];
                stiff_o[i] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) +
                                           0.5f * val * current_field(neighbor_e, ej_perm));
                stiff_o[j] += inv_rho_o * (-0.5f * val * current_field(owner_e, ei));
                stiff_n[nfd_j] += inv_rho_o * (0.5f * val * current_field(owner_e, ei));
              },
              [&](const int m, const int j, const real_t val) {
                int const nfd_j = face_connectivity_local.getNeighborFaceDof(f, j);
                int const em = face_to_elem_dof_depth[fid_o][j][m];
                int const ej = face_to_elem_dof[fid_o][j];
                int const ej_perm = face_to_elem_dof[fid_n][nfd_j];
                norm_o[m] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) +
                                          0.5f * val * current_field(neighbor_e, ej_perm));
                stiff_o[j] += inv_rho_o * (-0.5f * val * current_field(owner_e, em));
                stiff_n[nfd_j] += inv_rho_o * (0.5f * val * current_field(owner_e, em));
              });

          // --- Neighbor side (outward normal = -normal[]) ---
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, neighbor_coords, fid_n, neg_normal,
              [&](const int i, const int j, const real_t val) {
                int const ofd_j = face_connectivity_local.getOwnerFaceDof(f, j);
                int const ei = face_to_elem_dof[fid_n][i];
                int const ej = face_to_elem_dof[fid_n][j];
                int const ej_perm = face_to_elem_dof[fid_o][ofd_j];
                stiff_n[i] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) +
                                           0.5f * val * current_field(owner_e, ej_perm));
                stiff_n[j] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ei));
                stiff_o[ofd_j] += inv_rho_n * (0.5f * val * current_field(neighbor_e, ei));
              },
              [&](const int m, const int j, const real_t val) {
                int const ofd_j = face_connectivity_local.getOwnerFaceDof(f, j);
                int const em = face_to_elem_dof_depth[fid_n][j][m];
                int const ej = face_to_elem_dof[fid_n][j];
                int const ej_perm = face_to_elem_dof[fid_o][ofd_j];
                norm_n[m] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) +
                                          0.5f * val * current_field(owner_e, ej_perm));
                stiff_n[j] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, em));
                stiff_o[ofd_j] += inv_rho_n * (0.5f * val * current_field(neighbor_e, em));
              });

          // The normal lines are off-face, so they bypass the face-sized flush below. Both callbacks
          // always fire with j == q, so the line is the one through face dof q on each side, in that
          // side's own face numbering.
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeBoundaryDampingAndInterfaceFlux_Team(
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

  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMember = typename TeamPolicyType::member_type;
  using ScratchView1D =
      Kokkos::View<float *, Kokkos::LayoutRight, ExecSpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ScratchView2D = Kokkos::View<float **, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // One team (one warp) per face, its threads spread over the face's (ORDER+1)^2 quadrature
  // points. A TeamThreadRange is correct whatever team size the backend hands back.
  TeamPolicyType policy(n_iter, kFaceTeamSize);
  policy.set_scratch_size(0, Kokkos::PerTeam(ScratchView2D::shmem_size(kFaceScratchRows, knumNodesPerFace) * 2 +
                                            ScratchView1D::shmem_size(kFaceGeomFloats) +
                                            ScratchView1D::shmem_size(2)));

  Kokkos::parallel_for(
      "DG BoundaryDamping+InterfaceFlux Team", policy, KOKKOS_LAMBDA(const TeamMember &team) {
        int const f = list_on ? list_local[team.league_rank()] : team.league_rank();

        // Allocated unconditionally, before the branch below, so every thread of the team walks the
        // same scratch bump-allocation sequence.
        ScratchView2D priv_o(team.team_scratch(0), kFaceScratchRows, knumNodesPerFace);
        ScratchView2D priv_n(team.team_scratch(0), kFaceScratchRows, knumNodesPerFace);
        ScratchView1D geom(team.team_scratch(0), kFaceGeomFloats);
        ScratchView1D penalty(team.team_scratch(0), 2);

        // The discretization primitives take C-array references, so the geometry buffer is viewed as
        // the fixed-size arrays they expect. Layout: faceCoords | owner_coords | neighbor_coords.
        auto &faceCoords = *reinterpret_cast<float(*)[4][3]>(geom.data());
        auto &owner_coords = *reinterpret_cast<float(*)[8][3]>(geom.data() + 12);
        auto &neighbor_coords = *reinterpret_cast<float(*)[8][3]>(geom.data() + 36);

        // f is a per-team value, so this branch is uniform across the team: the whole team takes the
        // same side and the early return leaves nobody behind at a later barrier.
        if (face_connectivity_local.isBoundaryFace(f)) {
          int const e = face_connectivity_local.elemOwner(f);
          int const faceId = face_connectivity_local.localFaceOwner(f);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 4), [&](const int j) {
            int const gni =
                face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
          });
          team.team_barrier();

          real_t const inv_vp = 1.0f / mesh_local.getModelVpOnElement(e);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, knumNodesPerFace), [&](const int i) {
            int const ei = face_to_elem_dof[faceId][i];
            ATOMICADD(damp_local_view(e, ei), inv_vp * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords));
          });
          return;
        }

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        auto const eIdx_o = mesh_local.elementIndex(owner_e);
        auto const eIdx_n = mesh_local.elementIndex(neighbor_e);

        // 20 items: the 4 face corner nodes, then the owner and neighbor element vertices. Spreading
        // the gather over the team is what makes it coalesced.
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 20), [&](const int p) {
          if (p < 4) {
            int const gni =
                face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(p));
            for (int d = 0; d < 3; ++d) faceCoords[p][d] = mesh_local.nodeCoord(gni, d);
          } else if (p < 12) {
            int const v = p - 4;
            mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx_o, v % 2, (v / 2) % 2, v / 4), owner_coords[v]);
          } else {
            int const v = p - 12;
            mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx_n, v % 2, (v / 2) % 2, v / 4), neighbor_coords[v]);
          }
        });

        // Every row is zeroed exactly once, by whichever thread the range hands it to, so the
        // reduction below is well defined even when the backend gave the team fewer threads than
        // kFaceTeamSize.
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kFaceScratchRows), [&](const int r) {
          for (int i = 0; i < knumNodesPerFace; ++i) {
            priv_o(r, i) = 0.0f;
            priv_n(r, i) = 0.0f;
          }
        });
        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 2), [&](const int s) {
          real_t const face_area = computeFaceArea(faceCoords);
          penalty(s) = (s == 0) ? computeSIPGPenaltyFromArea<ORDER>(face_area, owner_coords, penalty_local)
                                : computeSIPGPenaltyFromArea<ORDER>(face_area, neighbor_coords, penalty_local);
        });
        team.team_barrier();

        real_t const inv_rho_o = 1.0f / mesh_local.getModelRhoOnElement(owner_e);
        real_t const inv_rho_n = 1.0f / mesh_local.getModelRhoOnElement(neighbor_e);

        float normal[3];
        mesh_local.faceNormal(owner_e, static_cast<model::CubicFace>(fid_o), normal);
        float const neg_normal[3] = {-normal[0], -normal[1], -normal[2]};

        // Groups of kFaceGroupSize lane-adjacent threads share one accumulator row: halves the
        // scratch, which is the binding occupancy constraint, for the price of a shared-memory
        // atomic instead of a plain add.
        int const group_id = team.team_rank() / kFaceGroupSize;

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, knumNodesPerFace), [&](const int q) {
          // Face-normal accumulators, private to the thread: these dofs are off-face, so they cannot
          // go through the face-sized rows. Only ORDER+1 floats each.
          float norm_o[ORDER + 1] = {0};
          float norm_n[ORDER + 1] = {0};

          // Normal-contracted callbacks: the discretization folds sum_k C_ijk * n_k, so each
          // contribution fires once instead of once per physical direction. That divides the
          // shared-memory atomics below -- the kernel's dominant cost -- by three.
          // --- Owner side (outward normal = normal[]) ---
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, owner_coords, fid_o, normal,
              [&](const int i, const int j, const real_t val) {
                int const nfd_j = face_connectivity_local.getNeighborFaceDof(f, j);
                int const ei = face_to_elem_dof[fid_o][i];
                int const ej = face_to_elem_dof[fid_o][j];
                int const ej_perm = face_to_elem_dof[fid_n][nfd_j];
                ATOMICADD(priv_o(group_id, i), inv_rho_o * (-0.5f * val * current_field(owner_e, ej) +
                                                            0.5f * val * current_field(neighbor_e, ej_perm)));
                ATOMICADD(priv_o(group_id, j), inv_rho_o * (-0.5f * val * current_field(owner_e, ei)));
                ATOMICADD(priv_n(group_id, nfd_j), inv_rho_o * (0.5f * val * current_field(owner_e, ei)));
              },
              [&](const int m, const int j, const real_t val) {
                int const nfd_j = face_connectivity_local.getNeighborFaceDof(f, j);
                int const em = face_to_elem_dof_depth[fid_o][j][m];
                int const ej = face_to_elem_dof[fid_o][j];
                int const ej_perm = face_to_elem_dof[fid_n][nfd_j];
                norm_o[m] += inv_rho_o * (-0.5f * val * current_field(owner_e, ej) +
                                          0.5f * val * current_field(neighbor_e, ej_perm));
                ATOMICADD(priv_o(group_id, j), inv_rho_o * (-0.5f * val * current_field(owner_e, em)));
                ATOMICADD(priv_n(group_id, nfd_j), inv_rho_o * (0.5f * val * current_field(owner_e, em)));
              });

          // --- Neighbor side (outward normal = -normal[]) ---
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, neighbor_coords, fid_n, neg_normal,
              [&](const int i, const int j, const real_t val) {
                int const ofd_j = face_connectivity_local.getOwnerFaceDof(f, j);
                int const ei = face_to_elem_dof[fid_n][i];
                int const ej = face_to_elem_dof[fid_n][j];
                int const ej_perm = face_to_elem_dof[fid_o][ofd_j];
                ATOMICADD(priv_n(group_id, i), inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) +
                                                            0.5f * val * current_field(owner_e, ej_perm)));
                ATOMICADD(priv_n(group_id, j), inv_rho_n * (-0.5f * val * current_field(neighbor_e, ei)));
                ATOMICADD(priv_o(group_id, ofd_j), inv_rho_n * (0.5f * val * current_field(neighbor_e, ei)));
              },
              [&](const int m, const int j, const real_t val) {
                int const ofd_j = face_connectivity_local.getOwnerFaceDof(f, j);
                int const em = face_to_elem_dof_depth[fid_n][j][m];
                int const ej = face_to_elem_dof[fid_n][j];
                int const ej_perm = face_to_elem_dof[fid_o][ofd_j];
                norm_n[m] += inv_rho_n * (-0.5f * val * current_field(neighbor_e, ej) +
                                          0.5f * val * current_field(owner_e, ej_perm));
                ATOMICADD(priv_n(group_id, j), inv_rho_n * (-0.5f * val * current_field(neighbor_e, em)));
                ATOMICADD(priv_o(group_id, ofd_j), inv_rho_n * (0.5f * val * current_field(neighbor_e, em)));
              });

          // The normal lines are off-face, so they bypass the column reduction below. Both callbacks
          // always fire with j == q, so the line is the one through face dof q on each side, in that
          // side's own face numbering.
          for (int m = 0; m <= ORDER; ++m) {
            ATOMICADD(stiff_local_view(owner_e, face_to_elem_dof_depth[fid_o][q][m]), norm_o[m]);
            ATOMICADD(stiff_local_view(neighbor_e, face_to_elem_dof_depth[fid_n][q][m]), norm_n[m]);
          }
        });
        team.team_barrier();

        // Thread i owns face dof i: it sums column i across the team's rows, so neither the
        // reduction nor the penalty term needs an atomic. Only the write-back does, because several
        // faces share an element.
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, knumNodesPerFace), [&](const int i) {
          float sum_o = 0.0f;
          float sum_n = 0.0f;
          for (int t = 0; t < kFaceScratchRows; ++t) {
            sum_o += priv_o(t, i);
            sum_n += priv_n(t, i);
          }

          real_t const damping_i = INTEGRAL_TYPE::computeDampingTerm(i, faceCoords);

          int const ei_o = face_to_elem_dof[fid_o][i];
          int const ei_o_perm = face_to_elem_dof[fid_n][face_connectivity_local.getNeighborFaceDof(f, i)];
          sum_o += penalty(0) * damping_i * (current_field(owner_e, ei_o) - current_field(neighbor_e, ei_o_perm));

          int const ei_n = face_to_elem_dof[fid_n][i];
          int const ei_n_perm = face_to_elem_dof[fid_o][face_connectivity_local.getOwnerFaceDof(f, i)];
          sum_n += penalty(1) * damping_i * (current_field(neighbor_e, ei_n) - current_field(owner_e, ei_n_perm));

          ATOMICADD(stiff_local_view(owner_e, ei_o), sum_o);
          ATOMICADD(stiff_local_view(neighbor_e, ei_n), sum_n);
        });
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
