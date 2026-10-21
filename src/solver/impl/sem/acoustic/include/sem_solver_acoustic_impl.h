#ifndef FUNTIDES_SOLVER_IMPL_SEM_ACOUSTIC_INCLUDE_SEM_SOLVER_ACOUSTIC_IMPL_H_
#define FUNTIDES_SOLVER_IMPL_SEM_ACOUSTIC_INCLUDE_SEM_SOLVER_ACOUSTIC_IMPL_H_

/**
 * @file sem_solver_acoustic_impl.h
 * @brief Acoustic element kernels of SEMsolver.
 *
 * Out-of-line member definitions of the acoustic element kernels:
 * computeElementContributions_Acoustic{,_Gemm,_Flat,_FlatZ,_TeamZ} and
 * computeAttenuationContributionsAcoustic.
 *
 * Fragment contract: include after sem_solver_kernels.h, which defines the
 * detail:: capability probes and the shared z-deformed geometry helpers this
 * code uses. The explicit-instantiation units include the umbrella
 * (sem_solver_impl.h), never this file directly.
 */

#include <data_type.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>

#include "Integrals.h"
#include "mesh_type_traits.h"
#include "sem_solver.h"

namespace solver {

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_Gemm(
    const DataType& data) {
  if constexpr (detail::has_team_gemm<INTEGRAL_TYPE>::value) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    constexpr int kStride = INTEGRAL_TYPE::numNodes * 6;

    auto mesh_local = m_mesh;
    int const dim = mesh_local.getOrder() + 1;
    int const pointsPerElem = dim * dim * dim;
    int const nElemsFull = mesh_local.getNumberOfElements();

    // W = w * alpha * B is precomputed once for every element, kStride values per element.
    // It runs on the first call (the benchmark warmup), so it stays out of the timed loop.
    if (!gemmMetricsReady_) {
      gemmMetrics_ = allocateVector<vectorReal>(static_cast<size_t>(nElemsFull) * kStride, "gemmMetrics");
      auto W_global = gemmMetrics_;
      auto mesh_pc = m_mesh;
      Kokkos::parallel_for(
          "Gemm Precompute Metrics", Kokkos::RangePolicy<ExecSpace>(0, nElemsFull), KOKKOS_LAMBDA(const int e) {
            float cornerCoords[8][3];
            auto const eIdx = mesh_pc.elementIndex(e);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_pc.vertexCoords(mesh_pc.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);

            real_t inv_density = 0.0f;
            if constexpr (!IS_MODEL_ON_NODES) {
              inv_density = 1.0f / mesh_pc.getModelRhoOnElement(e);
            }

            INTEGRAL_TYPE::computeElementMetrics(
                cornerCoords,
                [&](const int qa, const int qb, const int qc) -> real_t {
                  if constexpr (IS_MODEL_ON_NODES) {
                    int const g = mesh_pc.globalNodeIndex(e, qa, qb, qc);
                    return 1.0f / mesh_pc.getModelRhoOnNodes(g);
                  } else {
                    return inv_density;
                  }
                },
                W_global.data() + static_cast<size_t>(e) * kStride);
          });
      Kokkos::fence();
      gemmMetricsReady_ = true;
    }

    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : nElemsFull;

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) local_workVectorsGlobal[f] = workVectorsGlobal_[f];

    auto W_global = gemmMetrics_;

    // 1D derivative operator, built on host, captured by value (no device View).
    Kokkos::Array<real_t, INTEGRAL_TYPE::num1dNodes * INTEGRAL_TYPE::num1dNodes> D_arr;
    INTEGRAL_TYPE::fillDerivativeMatrix(D_arr.data());

    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView2D = Kokkos::View<float**, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    TeamPolicyType policy(n_iter, Kokkos::AUTO);
    // Scratch holds the gather/scatter buffers and the 12 [N][N^2] GEMM tensors.
    // W lives in global memory, so it needs no scratch: less shared memory than
    // the streaming variant, which also relaxes the occupancy limit slightly.
    size_t bytes_fields = ScratchView2D::shmem_size(kNumFields, pointsPerElem) * 2;
    size_t bytes_gemm = ScratchView2D::shmem_size(INTEGRAL_TYPE::num1dNodes, INTEGRAL_TYPE::numNodesPerFace) * 12;
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_gemm));

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic Gemm", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          ScratchView2D localFields(team.team_scratch(0), kNumFields, pointsPerElem);
          ScratchView2D localWork(team.team_scratch(0), kNumFields, pointsPerElem);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pointsPerElem), [&](const int localIdx) {
            int i = localIdx % dim;
            int j = (localIdx / dim) % dim;
            int k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            for (int f = 0; f < kNumFields; ++f) {
              localFields(f, localIdx) = data.getCurrentField(f)(globalIdx);
              localWork(f, localIdx) = 0.0f;
            }
          });
          team.team_barrier();

          real_t const* W_ptr = W_global.data() + static_cast<size_t>(elementNumber) * kStride;
          INTEGRAL_TYPE::computeStiffnessOperatorTeamVector(team, &localFields(0, 0), &localWork(0, 0), W_ptr,
                                                            D_arr.data());

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pointsPerElem), [&](const int localIdx) {
            int i = localIdx % dim;
            int j = (localIdx / dim) % dim;
            int k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork(f, localIdx));
            }
          });
        });
  } else {
    throw std::runtime_error(
        "computeElementContributions_Acoustic_Gemm: INTEGRAL_TYPE has no GEMM path "
        "(use Qk_Hexahedron_Tensorial).");
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic(
    const DataType& data) {
  if constexpr (detail::has_team_gemm<INTEGRAL_TYPE>::value) {
    computeElementContributions_Acoustic_Gemm(data);
  } else {
    if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
      prepareZDeformedGeometry();
      if (zDeformedEnabled_) {
        computeElementContributions_Acoustic_TeamZ(data);
        return;
      }
    }
    computeElementContributions_Acoustic_Flat(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_FlatZ(
    const DataType& data) {
  if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    auto mesh_local = m_mesh;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = mesh_local.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto geom = zDeformedGeom_;
    float const* invRhoNodes = zDeformedInvRhoNodeTable();
    int const* elemNodes = zDeformedNodeTable();
    auto force = workVectorsGlobal_[0];

    // Input and output arrays of one element live in registers: 3 blocks per SM leave room for both.
    using Policy = Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, 3>>;

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic FlatZ", detail::lightWeight(Policy(0, n_iter)),
        KOKKOS_LAMBDA(const int _loop_idx) {
          int const e = list_on ? list_local[_loop_idx] : _loop_idx;
          constexpr int dim = ORDER + 1;

          float u[kPointsPerElement];
          for (int k = 0; k < dim; ++k)
            for (int j = 0; j < dim; ++j)
              for (int i = 0; i < dim; ++i)
                u[i + j * dim + k * dim * dim] =
                    data.getCurrentField(0)(elemNodes[(i + j * dim + k * dim * dim) * n_iter + _loop_idx]);

          real_t const J00 = geom(e);
          real_t const J11 = geom(nElems + e);
          real_t Z[8];
          for (int k = 0; k < 8; ++k) Z[k] = geom((2 + k) * nElems + e);
          real_t const invRhoElem = geom((kZGeomStride - 1) * nElems + e);

          auto get_alpha = [&](const int q) -> real_t {
            if (IS_MODEL_ON_NODES) return invRhoNodes[elemNodes[q * n_iter + _loop_idx]];
            return invRhoElem;
          };

          float v[kPointsPerElement];
          INTEGRAL_TYPE::computeStiffnessTermSumFactZDeformedScatter(J00, J11, Z, u, v, get_alpha);

          for (int k = 0; k < dim; ++k)
            for (int j = 0; j < dim; ++j)
              for (int i = 0; i < dim; ++i)
                ATOMICADD(force[elemNodes[(i + j * dim + k * dim * dim) * n_iter + _loop_idx]],
                          v[i + j * dim + k * dim * dim]);
        });
  } else {
    throw std::runtime_error(
        "computeElementContributions_Acoustic_FlatZ: INTEGRAL_TYPE has no z-deformed sum-factorization path.");
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_TeamZ(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kAcoustic ||
                !detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    throw std::runtime_error("computeElementContributions_Acoustic_TeamZ: needs an acoustic z-deformed setup.");
  } else {
    constexpr int n = ORDER + 1;
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    constexpr int kMinTeamsPerSM = std::clamp(1024 / kPreferredTeamSize, 1, 16);

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<kPreferredTeamSize, kMinTeamsPerSM>>;
    using TeamMember = typename TeamPolicyType::member_type;

    // One thread per node and one per geometry entry are required. The last geometry entry, 1/rho of
    // the element, is only read when the model lives on elements.
    constexpr int kGeomLoaded = IS_MODEL_ON_NODES ? kZGeomStride - 1 : kZGeomStride;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    if (team_size < std::max(kPointsPerElement, kGeomLoaded)) {
      computeElementContributions_Acoustic_FlatZ(data);
      return;
    }

    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = m_mesh.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto zgeom = zDeformedGeom_;
    int const* elemNodes = zDeformedNodeTable(true);
    float const* invRhoNodes = zDeformedInvRhoNodeTable();
    float const* basisTabGlobal = zDeformedBasisTab_.data();
    auto force = workVectorsGlobal_[0];

    // Same choice as the TTI team kernel: up to 8 elements per team, at least four waves of teams.
    constexpr int kMaxElemsPerTeam = 8;
    constexpr int kMinWaves = 4;
    int const resident_teams = std::max<int>(1, ExecSpace::concurrency() / kPreferredTeamSize);
    int const elems_per_team = std::clamp(n_iter / (kMinWaves * resident_teams), 1, kMaxElemsPerTeam);
    TeamPolicyType policy((n_iter + elems_per_team - 1) / elems_per_team, team_size);

    // Shared scratch, every block on a 16-byte boundary: rows of D padded to a multiple of 4, element
    // geometry, nodal values, then the three fluxes.
    constexpr int kRowChunks = (n + 3) / 4;
    constexpr int kRowStride = 4 * kRowChunks;
    constexpr int kGeomPadded = ((kZGeomStride + 3) / 4) * 4;
    static_assert(kGeomPadded == 12, "the geometry is read as three 16-byte blocks");
    constexpr int kGeomOffset = n * kRowStride;
    constexpr int kValueOffset = kGeomOffset + kGeomPadded;
    constexpr int kFluxOffset = kValueOffset + kPointsPerElement;
    constexpr int kScratchFloats = kFluxOffset + 3 * kPointsPerElement;
    constexpr size_t kScratchBytes = kScratchFloats * sizeof(float);
    policy.set_scratch_size(0, Kokkos::PerTeam(kScratchBytes + alignof(detail::Float4)));

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic TeamZ", detail::lightWeight(policy),
        KOKKOS_LAMBDA(const TeamMember& team) {
          float* scratch =
              static_cast<float*>(team.team_scratch(0).get_shmem_aligned(kScratchBytes, alignof(detail::Float4)));
          // Captured here rather than first inside an if constexpr branch, which nvcc rejects.
          float const* inv_rho_nodes = invRhoNodes;
          auto const geom = zgeom;
          float* rows = scratch;
          float* geom_sh = scratch + kGeomOffset;
          float* u = scratch + kValueOffset;
          float* G = scratch + kFluxOffset;
          detail::Float4 const* rows4 = reinterpret_cast<detail::Float4 const*>(rows);
          detail::Float4 const* geom4 = reinterpret_cast<detail::Float4 const*>(geom_sh);

          // Thread q owns point q, which is also node q.
          int const q = team.team_rank();
          bool const owns = q < kPointsPerElement;
          int const qa = q % n;
          int const qb = (q / n) % n;
          int const qc = q / (n * n);

          // Filled before the first barrier, which covers it.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n * kRowStride), [&](const int idx) {
            int const r = idx / kRowStride;
            int const p = idx % kRowStride;
            rows[idx] = p < n ? basisTabGlobal[r * n + p] : 0.0f;
          });

          // The same for every element of the team: basis derivatives along the three lines through q,
          // the quadrature weight and the interpolation coordinates of q.
          float d_xi[n];
          float d_eta[n];
          float d_zeta[n];
          float weight = 0.0f;
          float alpha[3] = {0.0f, 0.0f, 0.0f};
          if (owns) {
            for (int i = 0; i < n; ++i) {
              d_xi[i] = basisTabGlobal[i * n + qa];
              d_eta[i] = basisTabGlobal[i * n + qb];
              d_zeta[i] = basisTabGlobal[i * n + qc];
            }
            weight = basisTabGlobal[n * n + qa] * basisTabGlobal[n * n + qb] * basisTabGlobal[n * n + qc];
            alpha[0] = basisTabGlobal[n * n + n + qa];
            alpha[1] = basisTabGlobal[n * n + n + qb];
            alpha[2] = basisTabGlobal[n * n + n + qc];
          }

          int const first_rank = team.league_rank() * elems_per_team;
          int const n_here = Kokkos::min(elems_per_team, n_iter - first_rank);

          // Two-stage pipeline held in registers: node index and element number two elements ahead, values
          // and one geometry entry one element ahead, so that a load never waits on its index.
          int ahead_node = 0;
          int ahead_e = 0;
          int next_node = 0;
          float next_u = 0.0f;
          float next_inv_rho = 0.0f;
          float next_geom = 0.0f;
          auto const load_index = [&](int rank) {
            if (owns) ahead_node = elemNodes[rank * kPointsPerElement + q];
            if (q < kGeomLoaded) ahead_e = list_on ? list_local[rank] : rank;
          };
          auto const load_values = [&]() {
            if (owns) {
              next_node = ahead_node;
              next_u = data.getCurrentField(0)(next_node);
              if constexpr (IS_MODEL_ON_NODES) next_inv_rho = inv_rho_nodes[next_node];
            }
            if (q < kGeomLoaded) next_geom = geom(q * nElems + ahead_e);
          };
          load_index(first_rank);
          load_values();
          if (n_here > 1) load_index(first_rank + 1);

          for (int k = 0; k < n_here; ++k) {
            int const node = next_node;
            float inv_rho = next_inv_rho;

            if (owns) u[q] = next_u;
            if (q < kGeomLoaded) geom_sh[q] = next_geom;
            team.team_barrier();

            if (k + 1 < n_here) load_values();
            if (k + 2 < n_here) load_index(first_rank + k + 2);

            // Entries J00, J11, the eight vertex z and 1/rho of the element, in that order.
            detail::Float4 const g0 = geom4[0];
            detail::Float4 const g1 = geom4[1];
            detail::Float4 const g2 = geom4[2];
            float const J00 = g0.v[0];
            float const J11 = g0.v[1];
            float const Z[8] = {g0.v[2], g0.v[3], g1.v[0], g1.v[1], g1.v[2], g1.v[3], g2.v[0], g2.v[1]};
            if constexpr (!IS_MODEL_ON_NODES) inv_rho = g2.v[2];

            if (owns) {
              float du[3] = {0.0f, 0.0f, 0.0f};
              for (int i = 0; i < n; ++i) {
                du[0] += d_xi[i] * u[i + qb * n + qc * n * n];
                du[1] += d_eta[i] * u[qa + i * n + qc * n * n];
                du[2] += d_zeta[i] * u[qa + qb * n + i * n * n];
              }

              float dZ[3][4];
              INTEGRAL_TYPE::zDeformedHalfDifferences(Z, dZ);
              float invJ[3][3];
              float const detJ = INTEGRAL_TYPE::invJacobianZDeformedBilinear(alpha, J00, J11, dZ, invJ);
              float const a = invJ[0][0];
              float const b = invJ[1][1];
              float const c = invJ[2][0];
              float const d = invJ[2][1];
              float const ez = invJ[2][2];
              // B = det(J) J^-1 J^-T, whose xy entry is zero.
              float const scale = weight * inv_rho * detJ;
              float const B0 = a * a;
              float const B1 = b * b;
              float const B2 = c * c + d * d + ez * ez;
              float const B3 = b * d;
              float const B4 = a * c;
              G[q] = scale * (B0 * du[0] + B4 * du[2]);
              G[kPointsPerElement + q] = scale * (B1 * du[1] + B3 * du[2]);
              G[2 * kPointsPerElement + q] = scale * (B4 * du[0] + B3 * du[1] + B2 * du[2]);
            }
            team.team_barrier();

            // Transpose: node (qa, qb, qc) gathers D(qa, p) G_xi(p, qb, qc), and so on.
            if (owns) {
              float v = 0.0f;
              for (int c = 0; c < kRowChunks; ++c) {
                detail::Float4 const ra = rows4[qa * kRowChunks + c];
                detail::Float4 const rb = rows4[qb * kRowChunks + c];
                detail::Float4 const rc = rows4[qc * kRowChunks + c];
                for (int j = 0; j < 4; ++j) {
                  int const p = 4 * c + j;
                  if (p < n) {
                    v += ra.v[j] * G[p + qb * n + qc * n * n];
                    v += rb.v[j] * G[kPointsPerElement + qa + p * n + qc * n * n];
                    v += rc.v[j] * G[2 * kPointsPerElement + qa + qb * n + p * n * n];
                  }
                }
              }
              ATOMICADD(force[node], v);
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  using Policy = Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>;

  Kokkos::parallel_for(
      "Solver Element Contribution Acoustic", Policy(0, n_iter), KOKKOS_LAMBDA(const int _loop_idx) {
        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        constexpr int dim = ORDER + 1;

        float localFields[kNumFields][kPointsPerElement];
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int k = 0; k < dim; ++k) {
          for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        real_t inv_density = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density = 1.0f / mesh_local.getModelRhoOnElement(elementNumber);
        }

        INTEGRAL_TYPE::computeStiffnessTermSumFact(
            cornerCoords, localFields[0], localWork[0], [&](const int qa, const int qb, const int qc) -> real_t {
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                return 1.0f / mesh_local.getModelRhoOnNodes(gIndex);
              } else {
                return inv_density;
              }
            });

        for (int k = 0; k < dim; ++k) {
          for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) {
                ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsAcoustic(
    const DataType& data) {
  auto mesh_local = m_mesh;
  Kokkos::Array<vectorReal, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Acoustic",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWorkA[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
            }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        real_t inv_density_q = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density_q =
              1.0f / (mesh_local.getModelRhoOnElement(elementNumber) * mesh_local.getModelQpOnElement(elementNumber));
        }

        INTEGRAL_TYPE::computeStiffnessTerm(
            cornerCoords,
            [&](const int qa, const int qb, const int qc) {
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                inv_density_q = 1.0f / (mesh_local.getModelRhoOnNodes(gIndex) * mesh_local.getModelQpOnNodes(gIndex));
              }
            },
            [&](const int i, const int j, const real_t val) {
              localWorkA[0][i] += inv_density_q * val * localFields[0][j];
            });

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
                ATOMICADD(local_attenuationWorkVectorsGlobal[f][globalIdx], localWorkA[f][localIdx]);
            }
      });
}

}  // namespace solver

#endif  // FUNTIDES_SOLVER_IMPL_SEM_ACOUSTIC_INCLUDE_SEM_SOLVER_ACOUSTIC_IMPL_H_
