#ifndef FUNTIDES_SOLVER_IMPL_PHYSICS_ELASTIC_INCLUDE_SEM_SOLVER_ELASTIC_KERNELS_H_
#define FUNTIDES_SOLVER_IMPL_PHYSICS_ELASTIC_INCLUDE_SEM_SOLVER_ELASTIC_KERNELS_H_

/**
 * @file sem_solver_elastic_kernels.h
 * @brief Elastic element kernels of SEMsolver (physics-policy migration, step 2).
 *
 * Out-of-line member definitions moved verbatim from sem_solver_kernels.h:
 * computeAttenuationContributionsElastic, computeElementContributions_
 * {Iso,Vti,Tti}{,_Flat,_Team,_TeamZ} and precomputeTtiTensorsOnNodes.
 * The split is positional (no line edited): kept+extracted is an exact
 * partition of the original file.
 *
 * Fragment contract: include after sem_solver_acoustic_kernels.h; the
 * explicit-instantiation units include the umbrella (sem_solver_impl.h),
 * never this file directly.
 */

#include <data_type.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>

#include "Integrals.h"
#include "elastic_flux.h"
#include "mesh_type_traits.h"
#include "sem_solver.h"

namespace solver {

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsElastic(
    const DataType& data) {
  if (anisotropyType_ != model::AnisotropyType::kIso) return;

  auto mesh_local = m_mesh;

  Kokkos::Array<vectorReal, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Elastic",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        constexpr int dim = ORDER + 1;
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

        // Per (p, r) pair of reference directions: the six distinct coefficients of the
        // isotropic stiffness contracted with the inverse Jacobian.
        // a0, a1, a2 multiply the diagonal displacement terms; a3, b0, b1 the xy, xz, yz couplings.
        struct CJPacked {
          float a0, a1, a2, a3;
          float b0, b1;
        };
        CJPacked CJflat[3 * 3];

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            cornerCoords,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              float vp, vs, rho, qp, qs;
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                vp = mesh_local.getModelVpOnNodes(gIndex);
                vs = mesh_local.getModelVsOnNodes(gIndex);
                rho = mesh_local.getModelRhoOnNodes(gIndex);
                qp = mesh_local.getModelQpOnNodes(gIndex);
                qs = mesh_local.getModelQsOnNodes(gIndex);
              } else {
                vp = mesh_local.getModelVpOnElement(elementNumber);
                vs = mesh_local.getModelVsOnElement(elementNumber);
                rho = mesh_local.getModelRhoOnElement(elementNumber);
                qp = mesh_local.getModelQpOnElement(elementNumber);
                qs = mesh_local.getModelQsOnElement(elementNumber);
              }

              float mu = rho * vs * vs;
              float lambda = rho * (vp * vp - 2.0f * vs * vs);
              float lambdap2mua = (lambda + 2.0f * mu) / qp;
              mu = mu / qs;
              lambda = lambdap2mua - 2.0f * mu;
              float lambda_plus_2mu = lambda + 2.0f * mu;

              for (int p = 0; p < 3; ++p) {
                float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
                for (int r = 0; r < 3; ++r) {
                  float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                  int const idx = p * 3 + r;
                  CJflat[idx].a0 = lambda_plus_2mu * Jp0 * Jr0 + mu * (Jp1 * Jr1 + Jp2 * Jr2);
                  CJflat[idx].a1 = mu * Jp0 * Jr0 + lambda_plus_2mu * Jp1 * Jr1 + mu * Jp2 * Jr2;
                  CJflat[idx].a2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) + lambda_plus_2mu * Jp2 * Jr2;
                  CJflat[idx].a3 = lambda * Jp0 * Jr1 + mu * Jp1 * Jr0;
                  CJflat[idx].b0 = lambda * Jp0 * Jr2 + mu * Jp2 * Jr0;
                  CJflat[idx].b1 = lambda * Jp1 * Jr2 + mu * Jp2 * Jr1;
                }
              }
            },
            [&](int i, int j, float val, const int p, const int r) {
              int const idx = p * 3 + r;
              float const uxj = localFields[0][j];
              float const uyj = localFields[1][j];
              float const uzj = localFields[2][j];
              localWorkA[0][i] += val * (CJflat[idx].a0 * uxj + CJflat[idx].a3 * uyj + CJflat[idx].b0 * uzj);
              localWorkA[1][i] += val * (CJflat[idx].a3 * uxj + CJflat[idx].a1 * uyj + CJflat[idx].b1 * uzj);
              localWorkA[2][i] += val * (CJflat[idx].b0 * uxj + CJflat[idx].b1 * uyj + CJflat[idx].a2 * uzj);
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Element Contribution Iso Flat",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;

        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float mu_e = 0.0f, lambda_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
          }

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork,
                                                        [&](int qa, int qb, int qc, float const(&J_inv)[3][3],
                                                            float const(&grad_u_ref)[3][3], float(&flux)[3][3]) {
                                                          float mu, lambda;
                                                          if constexpr (IS_MODEL_ON_NODES) {
                                                            int const gIndex =
                                                                mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                                                            float const vp = mesh_local.getModelVpOnNodes(gIndex);
                                                            float const vs = mesh_local.getModelVsOnNodes(gIndex);
                                                            float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                                                            mu = rho * vs * vs;
                                                            lambda = rho * (vp * vp - 2.0f * vs * vs);
                                                          } else {
                                                            mu = mu_e;
                                                            lambda = lambda_e;
                                                          }

                                                          flux::elasticFluxIso(J_inv, mu, lambda, grad_u_ref, flux);
                                                        });

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Iso_Flat(data);
  } else {
    computeElementContributions_Iso_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso_Team(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  // avoid extended __host__ __device__ lambda cannot first-capture variable in
  // constexpr-if context
  (void)mesh_local;
  (void)list_on;
  (void)list_local;
  (void)n_iter;
  (void)local_workVectorsGlobal;

  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;

    // A team owns one element: gather and scatter are spread over its threads,
    // which makes them coalesced, and the per-element buffers move from
    // per-thread local memory to shared memory.
    //
    // Every team range covers kPointsPerElement, so a larger team only adds
    // threads that idle at the barriers. Round up to a whole warp so no partial
    // warp idles where kPointsPerElement is not a multiple of 32, then clamp to
    // what the backend can actually provide: a host backend caps the team size
    // at its thread count, and asking for more aborts at launch. A
    // TeamThreadRange is correct at any team size.
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    // Constant-Jacobian meshes carry the element geometry (9 inverse-Jacobian
    // entries + the determinant) in a small per-team scratch buffer.
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;
    TeamPolicyType policy(n_iter, team_size);
    // One buffer carries the displacements in and the forces out: the stiffness
    // kernel overwrites it after the barrier that ends its read phase. Halving
    // this buffer is what lifts the occupancy off its shared-memory limit.
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Iso Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float mu_e = 0.0f, lambda_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
          }

          auto const iso_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            float mu, lambda;
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float const vp = mesh_local.getModelVpOnNodes(gIndex);
              float const vs = mesh_local.getModelVsOnNodes(gIndex);
              float const rho = mesh_local.getModelRhoOnNodes(gIndex);
              mu = rho * vs * vs;
              lambda = rho * (vp * vp - 2.0f * vs * vs);
            } else {
              mu = mu_e;
              lambda = lambda_e;
            }
            flux::elasticFluxIso(J_inv, mu, lambda, grad_u_ref, flux);
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              // invJacobianTransformation accumulates into J, so it must start at zero.
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), iso_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), iso_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Element Contribution Vti Flat",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;
        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float c11_e = 0, c12_e = 0, c13_e = 0, c33_e = 0, c44_e = 0, c66_e = 0;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            float const delta_e = mesh_local.getModelDeltaOnElement(elementNumber);
            float const epsilon_e = mesh_local.getModelEpsilonOnElement(elementNumber);
            float const gamma_e = mesh_local.getModelGammaOnElement(elementNumber);
            float const rho_vp2 = rho_e * vp_e * vp_e;
            float const rho_vs2 = rho_e * vs_e * vs_e;
            c33_e = rho_vp2;
            c44_e = rho_vs2;
            c11_e = rho_vp2 * (1.0f + 2.0f * epsilon_e);
            c66_e = rho_vs2 * (1.0f + 2.0f * gamma_e);
            float const vp2_vs2 = vp_e * vp_e - vs_e * vs_e;
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp_e * vp_e * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(
              cornerCoords, localFields, localWork,
              [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                  float(&flux)[3][3]) {
                float c11, c12, c13, c33, c44, c66;
                if constexpr (IS_MODEL_ON_NODES) {
                  int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  float const vp = mesh_local.getModelVpOnNodes(gIndex);
                  float const vs = mesh_local.getModelVsOnNodes(gIndex);
                  float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                  float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
                  float const epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
                  float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
                  float const rho_vp2 = rho * vp * vp;
                  float const rho_vs2 = rho * vs * vs;
                  c33 = rho_vp2;
                  c44 = rho_vs2;
                  c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
                  c66 = rho_vs2 * (1.0f + 2.0f * gamma);
                  float const vp2_vs2 = vp * vp - vs * vs;
                  c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;
                  c12 = c11 - 2.0f * c66;
                } else {
                  c11 = c11_e;
                  c12 = c12_e;
                  c13 = c13_e;
                  c33 = c33_e;
                  c44 = c44_e;
                  c66 = c66_e;
                }

                flux::elasticFluxVti(J_inv, c11, c12, c13, c33, c44, c66, grad_u_ref, flux);
              });

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Vti_Flat(data);
  } else {
    computeElementContributions_Vti_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti_Team(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  // avoid extended __host__ __device__ lambda cannot first-capture variable in
  // constexpr-if context
  (void)mesh_local;
  (void)list_on;
  (void)list_local;
  (void)local_workVectorsGlobal;

  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;
    // Warp-aligned preferred size, clamped to what the backend can provide (a
    // host backend caps it at its thread count; asking for more aborts).
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;

    TeamPolicyType policy(n_iter, team_size);
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Vti Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float c11_e = 0, c12_e = 0, c13_e = 0, c33_e = 0, c44_e = 0, c66_e = 0;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            float const delta_e = mesh_local.getModelDeltaOnElement(elementNumber);
            float const epsilon_e = mesh_local.getModelEpsilonOnElement(elementNumber);
            float const gamma_e = mesh_local.getModelGammaOnElement(elementNumber);
            float const rho_vp2 = rho_e * vp_e * vp_e;
            float const rho_vs2 = rho_e * vs_e * vs_e;
            c33_e = rho_vp2;
            c44_e = rho_vs2;
            c11_e = rho_vp2 * (1.0f + 2.0f * epsilon_e);
            c66_e = rho_vs2 * (1.0f + 2.0f * gamma_e);
            float const vp2_vs2 = vp_e * vp_e - vs_e * vs_e;
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp_e * vp_e * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          auto const vti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            float c11, c12, c13, c33, c44, c66;
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float const vp = mesh_local.getModelVpOnNodes(gIndex);
              float const vs = mesh_local.getModelVsOnNodes(gIndex);
              float const rho = mesh_local.getModelRhoOnNodes(gIndex);
              float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
              float const epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
              float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
              float const rho_vp2 = rho * vp * vp;
              float const rho_vs2 = rho * vs * vs;
              c33 = rho_vp2;
              c44 = rho_vs2;
              c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
              c66 = rho_vs2 * (1.0f + 2.0f * gamma);
              float const vp2_vs2 = vp * vp - vs * vs;
              c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;
              c12 = c11 - 2.0f * c66;
            } else {
              c11 = c11_e;
              c12 = c12_e;
              c13 = c13_e;
              c33 = c33_e;
              c44 = c44_e;
              c66 = c66_e;
            }
            flux::elasticFluxVti(J_inv, c11, c12, c13, c33, c44, c66, grad_u_ref, flux);
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), vti_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), vti_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

// Stores the compact TTI description of every node (see computeTtiCompact()) into cttiNodes_.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::precomputeTtiTensorsOnNodes() {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic || !IS_MODEL_ON_NODES) {
    return;
  } else {
    if (cttiNodesReady_) {
      return;
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    int const numNodes = m_mesh.getNumberOfNodes();
    cttiNodes_ = allocateArray2D<CttiView>(numNodes, flux::kTtiCompactSize, "cttiNodes");

    auto ctti = cttiNodes_;
    auto mesh_pc = m_mesh;
    Kokkos::parallel_for(
        "Tti Precompute Compact On Nodes", Kokkos::RangePolicy<ExecSpace>(0, numNodes), KOKKOS_LAMBDA(const int g) {
          float p[flux::kTtiCompactSize];
          computeTtiCompact(mesh_pc.getModelVpOnNodes(g), mesh_pc.getModelVsOnNodes(g), mesh_pc.getModelRhoOnNodes(g),
                            mesh_pc.getModelDeltaOnNodes(g), mesh_pc.getModelEpsilonOnNodes(g),
                            mesh_pc.getModelGammaOnNodes(g), mesh_pc.getModelPhiOnNodes(g),
                            mesh_pc.getModelThetaOnNodes(g), p);
          for (int k = 0; k < flux::kTtiCompactSize; ++k) ctti(g, k) = p[k];
        });
    Kokkos::fence();

    cttiNodesReady_ = true;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_Flat(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    precomputeTtiTensorsOnNodes();

    auto mesh_local = m_mesh;
    [[maybe_unused]] auto ctti_local = cttiNodes_;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    Kokkos::parallel_for(
        "Solver Element Contribution Tti Flat",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
        KOKKOS_LAMBDA(const int _loop_idx) {
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          int const dim = mesh_local.getOrder() + 1;
          float localFields[kNumFields][kPointsPerElement] = {{0}};
          float localWork[kNumFields][kPointsPerElement] = {{0}};

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
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

          // Captured here, not at first use: nvcc rejects an extended lambda that first-captures
          // a variable inside an if-constexpr branch.
          [[maybe_unused]] auto const ctti = ctti_local;

          // Model on elements only: one tensor for the whole element.
          [[maybe_unused]] float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          auto ttiFlux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                             float(&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
              float p[flux::kTtiCompactSize];
              for (int k = 0; k < flux::kTtiCompactSize; ++k) p[k] = ctti(gIndex, k);
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            } else {
              flux::elasticFluxTti(J_inv, CTTI, grad_u_ref, flux);
            }
          };

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork, ttiFlux);

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
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
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Tti_Flat(data);
  } else {
    computeElementContributions_Tti_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_Team(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    if constexpr (!HasConstantJacobian<MESH_TYPE>::value && IS_MODEL_ON_NODES &&
                  detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
      prepareZDeformedGeometry();
      if (zDeformedEnabled_) {
        precomputeTtiTensorsOnNodes();
        computeElementContributions_Tti_TeamZ(data);
        return;
      }
    }
    precomputeTtiTensorsOnNodes();

    auto mesh_local = m_mesh;
    [[maybe_unused]] auto ctti_local = cttiNodes_;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;
    // Warp-aligned preferred size, clamped to what the backend can provide (a
    // host backend caps it at its thread count; asking for more aborts).
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;

    TeamPolicyType policy(n_iter, team_size);
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Tti Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          // Captured here, not at first use: nvcc rejects an extended lambda that first-captures
          // a variable inside an if-constexpr branch.
          [[maybe_unused]] auto const ctti = ctti_local;

          // Model on elements only: one tensor for the whole element.
          [[maybe_unused]] float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          auto const tti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float p[flux::kTtiCompactSize];
              for (int k = 0; k < flux::kTtiCompactSize; ++k) p[k] = ctti(gIndex, k);
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            } else {
              flux::elasticFluxTti(J_inv, CTTI, grad_u_ref, flux);
            }
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), tti_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), tti_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_TeamZ(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic || !IS_MODEL_ON_NODES ||
                !detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    throw std::runtime_error("computeElementContributions_Tti_TeamZ: needs an elastic on-nodes z-deformed setup.");
  } else {
    auto ctti_local = cttiNodes_;
    // Two 16-byte loads per node instead of eight scalar ones: each scalar load of a warp
    // touches the same cache lines again.
    static_assert(flux::kTtiCompactSize == 8, "two 16-byte loads per node");
    if (ctti_local.stride(0) != flux::kTtiCompactSize ||
        reinterpret_cast<std::uintptr_t>(ctti_local.data()) % alignof(detail::Float4) != 0)
      throw std::runtime_error("computeElementContributions_Tti_TeamZ: cttiNodes_ not packed for 16-byte loads.");
    detail::Float4 const* ctti_quads = reinterpret_cast<detail::Float4 const*>(ctti_local.data());
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = m_mesh.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto zgeom = zDeformedGeom_;
    int const* elemNodes = zDeformedNodeTable();
    float const* basisTabGlobal = zDeformedBasisTab_.data();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    constexpr int dim = ORDER + 1;
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    // Caps registers at 64 per thread; at most 16 teams, the lowest resident-block limit per SM.
    constexpr int kMinTeamsPerSM = std::clamp(1024 / kPreferredTeamSize, 1, 16);

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<kPreferredTeamSize, kMinTeamsPerSM>>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    // Several consecutive elements per team, fewer on small meshes so that the league still fills the device.
    constexpr int kMaxElemsPerTeam = 8;
    constexpr int kMinWaves = 4;
    int const resident_teams = std::max<int>(1, ExecSpace::concurrency() / kPreferredTeamSize);
    // J00, J11 and the eight vertex z; the last entry of the geometry table is unused in elastic.
    constexpr int kGeomEntries = kZGeomStride - 1;
    bool const pipelined = team_size >= std::max(kPointsPerElement, kGeomEntries);
    int const elems_per_team = pipelined ? std::clamp(n_iter / (kMinWaves * resident_teams), 1, kMaxElemsPerTeam) : 1;
    TeamPolicyType policy((n_iter + elems_per_team - 1) / elems_per_team, team_size);
    // One float block: fields (then forces), fluxes, basis table and element geometry. The geometry is
    // read by the whole team as three 16-byte blocks, so it starts on a 16-byte boundary.
    constexpr int kFluxOffset = kNumFields * kPointsPerElement;
    constexpr int kTabOffset = kFluxOffset + 9 * kPointsPerElement;
    constexpr int kGeomOffset = ((kTabOffset + INTEGRAL_TYPE::kBasisTableSize + 3) / 4) * 4;
    constexpr int kGeomPadded = ((kGeomEntries + 3) / 4) * 4;
    static_assert(kGeomPadded == 12, "the geometry is read as three 16-byte blocks");
    constexpr int kScratchFloats = kGeomOffset + kGeomPadded;
    constexpr size_t kScratchBytes = kScratchFloats * sizeof(float);
    size_t const bytes = kScratchBytes + alignof(detail::Float4) + ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes));

    Kokkos::parallel_for(
        "Solver Element Contribution Tti TeamZ", detail::lightWeight(policy), KOKKOS_LAMBDA(const TeamMember& team) {
          float* scratch =
              static_cast<float*>(team.team_scratch(0).get_shmem_aligned(kScratchBytes, alignof(detail::Float4)));
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Displacements in, forces out, same storage.
          float* localFields = scratch;
          float* fluxScratch = scratch + kFluxOffset;
          float* basisTab = scratch + kTabOffset;
          float* geomSh = scratch + kGeomOffset;
          detail::Float4 const* geom4 = reinterpret_cast<detail::Float4 const*>(geomSh);

          int const first_rank = team.league_rank() * elems_per_team;
          int const n_here = Kokkos::min(elems_per_team, n_iter - first_rank);
          int const point = team.team_rank();
          bool const owns_point = pipelined && point < kPointsPerElement;
          bool const carries_geom = pipelined && point < kGeomEntries;

          // Same two-stage pipeline as computeElementContributions_Acoustic_TeamZ().
          int ahead_node = 0;
          int ahead_e = 0;
          int next_node = 0;
          float next_fields[kNumFields] = {};
          float next_geom = 0.0f;
          auto const load_index = [&](int rank) {
            // Indexed by the league position, so the gather does not wait for the element list.
            if (owns_point) ahead_node = elemNodes[rank * kPointsPerElement + point];
            if (carries_geom) ahead_e = list_on ? list_local[rank] : rank;
          };
          auto const load_values = [&]() {
            if (owns_point) {
              next_node = ahead_node;
              for (int f = 0; f < kNumFields; ++f) next_fields[f] = data.getCurrentField(f)(next_node);
            }
            if (carries_geom) next_geom = zgeom(point * nElems + ahead_e);
          };
          if (pipelined) {
            load_index(first_rank);
            load_values();
            if (n_here > 1) load_index(first_rank + 1);
          }

          for (int k = 0; k < n_here; ++k) {
            int const rank = first_rank + k;

            if (k == 0) {
              // Filled before the gather barrier, so it needs no barrier of its own.
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                                   [&](const int idx) { basisTab[idx] = basisTabGlobal[idx]; });
            }

            if (pipelined) {
              if (owns_point) {
                nodeIdx(point) = next_node;
                for (int f = 0; f < kNumFields; ++f) localFields[f * kPointsPerElement + point] = next_fields[f];
              }
              if (carries_geom) geomSh[point] = next_geom;
            } else {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
                int const globalIdx = elemNodes[rank * kPointsPerElement + localIdx];
                nodeIdx(localIdx) = globalIdx;
                for (int f = 0; f < kNumFields; ++f)
                  localFields[f * kPointsPerElement + localIdx] = data.getCurrentField(f)(globalIdx);
              });
              int const elementNumber = list_on ? list_local[rank] : rank;
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kGeomEntries),
                                   [&](const int j) { geomSh[j] = zgeom(j * nElems + elementNumber); });
            }
            team.team_barrier();

            if (k + 1 < n_here) load_values();
            if (k + 2 < n_here) load_index(rank + 2);

            // Entries J00, J11 and the eight vertex z, in that order.
            detail::Float4 const g0 = geom4[0];
            detail::Float4 const g1 = geom4[1];
            detail::Float4 const g2 = geom4[2];
            float const Z[8] = {g0.v[2], g0.v[3], g1.v[0], g1.v[1], g1.v[2], g1.v[3], g2.v[0], g2.v[1]};

            auto const tti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                      float(&flux)[3][3]) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float p[flux::kTtiCompactSize];
              detail::Float4 const lo = ctti_quads[2 * static_cast<size_t>(gIndex)];
              detail::Float4 const hi = ctti_quads[2 * static_cast<size_t>(gIndex) + 1];
              for (int j = 0; j < 4; ++j) {
                p[j] = lo.v[j];
                p[j + 4] = hi.v[j];
              }
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            };

            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeamZDeformed(team, g0.v[0], g0.v[1], Z, localFields,
                                                                       localFields, fluxScratch, tti_flux, basisTab);
            team.team_barrier();

            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
              int const globalIdx = nodeIdx(localIdx);
              for (int f = 0; f < kNumFields; ++f) {
                ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields[f * kPointsPerElement + localIdx]);
              }
            });
          }
        });
  }
}


}  // namespace solver

#endif  // FUNTIDES_SOLVER_IMPL_PHYSICS_ELASTIC_INCLUDE_SEM_SOLVER_ELASTIC_KERNELS_H_
