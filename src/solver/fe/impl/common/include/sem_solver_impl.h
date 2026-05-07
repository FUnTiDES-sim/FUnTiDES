#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#include <data_type.h>

#include <array>
#include <cstdlib>

#include "Integrals.h"
#include "mesh_type_traits.h"
#include "sem_solver.h"

namespace solver {
namespace fe {

//============================================================================
// computeFEInit - Initialize finite element structures
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& sponge_size, const bool surface_sponge,
    const float taper_delta) {
  if (auto* typed_mesh = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh = *typed_mesh;
  } else {
    throw std::runtime_error("Incompatible mesh type in solver");
  }

  sponge_size_[0] = sponge_size[0];
  sponge_size_[1] = sponge_size[1];
  sponge_size_[2] = sponge_size[2];
  surface_sponge_ = surface_sponge;
  taper_delta_ = taper_delta;

  allocateFEarrays();
  initFEarrays();

  if constexpr (MeshTypeTraits<MESH_TYPE>::value == utils::enums::meshType::kStruct) {
    float const hx = m_mesh.getHx();
    float const hy = m_mesh.getHy();
    float const hz = m_mesh.getHz();
    m_struct_metrics_.detJ  = hx * hy * hz / 8.0f;
    m_struct_metrics_.inv_jx = 2.0f / hx;
    m_struct_metrics_.inv_jy = 2.0f / hy;
    m_struct_metrics_.inv_jz = 2.0f / hz;
    m_has_struct_metrics_ = true;
  }

  // Compute Local Mass Matrix
  computeGlobalMassMatrix();
  // Compute Local Damping Matrix
  computeDampingMatrix();

  if (attenuationEnabled_ && nSls_ > 0) {
    float minQVal = std::numeric_limits<float>::max();
    for (int e = 0; e < m_mesh.getNumberOfElements(); ++e) {
      if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
        float q = IS_MODEL_ON_NODES ? m_mesh.getModelQpOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                                    : m_mesh.getModelQpOnElement(e);
        minQVal = std::min(minQVal, q);
      } else {
        float qp = IS_MODEL_ON_NODES ? m_mesh.getModelQpOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                                     : m_mesh.getModelQpOnElement(e);
        float qs = IS_MODEL_ON_NODES ? m_mesh.getModelQsOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                                     : m_mesh.getModelQsOnElement(e);
        minQVal = std::min(minQVal, std::min(qp, qs));
      }
    }

    for (int l = 0; l < nSls_; ++l) {
      if (slsAnelasticityCoefficients_.extent(0) > 0 && slsAnelasticityCoefficients_[l] < 0.0f) {
        slsAnelasticityCoefficients_[l] = 2.0f * minQVal / (std::max(1.0001f, minQVal) - 1.0f);
      }
    }
  }
}

//============================================================================
// Compute Forces (Phase 1)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeForces(const float& dt,
                                                                                           const int& timeSample,
                                                                                           Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, myData);
  FENCE
  computeElementContributions(myData);
  FENCE
  if (attenuationEnabled_ && nSls_ > 0) {
    computeAttenuationContributions(myData);
    FENCE
  }
}

//============================================================================
// Update Solution (Phase 2)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolution(const float& dt,
                                                                                            Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  updateFields(dt, myData);
  FENCE
}

//============================================================================
// resetGlobalVectors
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::resetGlobalVectors(int numNodes) {
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  std::array<std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>, kNumFields>
      local_attenuationWorkVectorsGlobal;

  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    if (has_attenuation) {
      local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];
    }
  }

  Kokkos::parallel_for(
      "Solver Reset GVector", numNodes, KOKKOS_LAMBDA(const int i) {
        for (int f = 0; f < kNumFields; ++f) {
          local_workVectorsGlobal[f][i] = 0;
          if (has_attenuation) {
            local_attenuationWorkVectorsGlobal[f][i] = 0;
          }
        }
      });
}

//============================================================================
// applyRHSTerm
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyRHSTerm(int timeSample, float dt,
                                                                                          const DataType& data) {
  int nb_rhs_element = data.getRhsElement().extent(0);
  auto mesh_local = m_mesh;  // Capture mesh for lambda

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_LAMBDA(const int i) {
        for (int z = 0; z < ORDER + 1; z++) {
          for (int y = 0; y < ORDER + 1; y++) {
            for (int x = 0; x < ORDER + 1; x++) {
              int localNodeId = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
              int nodeRHS = mesh_local.globalNodeIndex(data.getRhsElement()[i], x, y, z);

              for (int f = 0; f < kNumRhs; ++f) {
                float source = data.getRhsTerm(f)(i, timeSample) * data.getRhsWeights()(i, localNodeId);
                local_workVectorsGlobal[f](nodeRHS) -= source;
              }
            }
          }
        }
      });
}

//============================================================================
// computeElementContributions - DISPATCH
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions(
    const DataType& data) {
  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    if (anisotropyType_ == model::AnisotropyType::kIso) {
      computeElementContributions_Iso(data);
    } else if (anisotropyType_ == model::AnisotropyType::kVTI) {
      computeElementContributions_Vti(data);
    } else  // kTTI
    {
      computeElementContributions_Tti(data);
    }
  } else  // Acoustic - DISPATCH
  {
    computeElementContributions_Acoustic(data);
  }
}

//============================================================================
// computeElementContributions_Acoustic - ACOUSTIC
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  auto const struct_metrics_local = m_struct_metrics_;

  Kokkos::parallel_for(
      "Solver Element Contribution Acoustic", Kokkos::RangePolicy(0, n_iter), KOKKOS_LAMBDA(const int _loop_idx) {
        (void)struct_metrics_local;
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

        real_t inv_density = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density = 1.0f / mesh_local.getModelRhoOnElement(elementNumber);
        }

        if constexpr (MeshTypeTraits<MESH_TYPE>::value == utils::enums::meshType::kStruct) {
          INTEGRAL_TYPE::computeStiffnessTermSumFactCartesian(
              localFields[0], localWork[0], struct_metrics_local,
              [&](const int qa, const int qb, const int qc) -> real_t {
                if constexpr (IS_MODEL_ON_NODES) {
                  int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  return 1.0f / mesh_local.getModelRhoOnNodes(gIndex);
                } else {
                  return inv_density;
                }
              });
        } else {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
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

          }

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

//============================================================================
// computeAttenuationContributions - Assemble attenuation stiffness terms
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributions(
    const DataType& data) {
  if (!attenuationEnabled_ || nSls_ <= 0) return;

  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic)
    computeAttenuationContributionsAcoustic(data);
  else
    computeAttenuationContributionsElastic(data);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsAcoustic(
    const DataType& data) {
  auto mesh_local = m_mesh;
  Kokkos::Array<VECTOR_REAL_VIEW, kNumFields> local_attenuationWorkVectorsGlobal;
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsElastic(
    const DataType& data) {
  if (anisotropyType_ != model::AnisotropyType::kIso) return;

  auto mesh_local = m_mesh;

  Kokkos::Array<VECTOR_REAL_VIEW, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Elastic",
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

//============================================================================
// computeElementContributions_Iso - ISOTROPIC OPTIMIZED
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  auto const struct_metrics_local = m_struct_metrics_;

  Kokkos::parallel_for(
      "Solver Element Contribution Iso",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;
        (void)struct_metrics_local;

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
          float mu_e = 0.0f, lambda_e = 0.0f, lam2mu_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
            lam2mu_e = lambda_e + 2.0f * mu_e;
          }

          // Pre-gather physics per node to avoid scattered global reads inside the kernel.
          float mu_local[kPointsPerElement] = {0};
          float lam2mu_local[kPointsPerElement] = {0};
          if constexpr (IS_MODEL_ON_NODES) {
            for (int i = 0; i < dim; ++i)
              for (int j = 0; j < dim; ++j)
                for (int k = 0; k < dim; ++k) {
                  int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                  int const lIdx = i + j * dim + k * dim * dim;
                  float const vp = mesh_local.getModelVpOnNodes(gIdx);
                  float const vs = mesh_local.getModelVsOnNodes(gIdx);
                  float const rho = mesh_local.getModelRhoOnNodes(gIdx);
                  float const mu = rho * vs * vs;
                  float const lam = rho * (vp * vp - 2.0f * vs * vs);
                  mu_local[lIdx] = mu;
                  lam2mu_local[lIdx] = lam + 2.0f * mu;
                }
          }

          if constexpr (solver::fe::MeshTypeTraits<MESH_TYPE>::value == utils::enums::meshType::kStruct) {
            float const a = struct_metrics_local.inv_jx;
            float const b = struct_metrics_local.inv_jy;
            float const c = struct_metrics_local.inv_jz;
            float const aa = a * a, bb = b * b, cc = c * c;
            float const ab = a * b, ac = a * c, bc = b * c;
            auto const iso_func_cart = [&](int qa, int qb, int qc, float const (&)[3][3],
                                           float const (&g)[3][3], float (&flux)[3][3]) {
              float mu, lambda, lam2mu;
              if constexpr (IS_MODEL_ON_NODES) {
                int const lIdx = qa + qb * dim + qc * dim * dim;
                mu = mu_local[lIdx];
                lam2mu = lam2mu_local[lIdx];
                lambda = lam2mu - 2.0f * mu;
              } else {
                mu = mu_e;
                lambda = lambda_e;
                lam2mu = lam2mu_e;
              }
              flux[0][0] = lam2mu * aa * g[0][0] + lambda * (ab * g[1][1] + ac * g[2][2]);
              flux[0][1] = mu * (aa * g[0][1] + ab * g[1][0]);
              flux[0][2] = mu * (aa * g[0][2] + ac * g[2][0]);
              flux[1][0] = mu * (ab * g[0][1] + bb * g[1][0]);
              flux[1][1] = lambda * (ab * g[0][0] + bc * g[2][2]) + lam2mu * bb * g[1][1];
              flux[1][2] = mu * (bb * g[1][2] + bc * g[2][1]);
              flux[2][0] = mu * (ac * g[0][2] + cc * g[2][0]);
              flux[2][1] = mu * (bc * g[1][2] + cc * g[2][1]);
              flux[2][2] = lambda * (ac * g[0][0] + bc * g[1][1]) + lam2mu * cc * g[2][2];
            };
            INTEGRAL_TYPE::computeElasticStiffnessSumFactCartesian(struct_metrics_local, localFields, localWork,
                                                                   iso_func_cart);
          } else {
            auto const iso_func = [&](int qa, int qb, int qc, float const (&J_inv)[3][3],
                                      float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
              float mu, lambda, lam2mu;
              if constexpr (IS_MODEL_ON_NODES) {
                int const lIdx = qa + qb * dim + qc * dim * dim;
                mu = mu_local[lIdx];
                lam2mu = lam2mu_local[lIdx];
                lambda = lam2mu - 2.0f * mu;
              } else {
                mu = mu_e;
                lambda = lambda_e;
                lam2mu = lam2mu_e;
              }

              for (int p = 0; p < 3; ++p) {
                float const Jp0 = J_inv[p][0];
                float const Jp1 = J_inv[p][1];
                float const Jp2 = J_inv[p][2];
                flux[p][0] = 0.0f;
                flux[p][1] = 0.0f;
                flux[p][2] = 0.0f;
                for (int r = 0; r < 3; ++r) {
                  float const Jr0 = J_inv[r][0];
                  float const Jr1 = J_inv[r][1];
                  float const Jr2 = J_inv[r][2];
                  float const v0 = lam2mu * Jp0 * Jr0 + mu * (Jp1 * Jr1 + Jp2 * Jr2);
                  float const v1 = mu * Jp0 * Jr0 + lam2mu * Jp1 * Jr1 + mu * Jp2 * Jr2;
                  float const v2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) + lam2mu * Jp2 * Jr2;
                  float const v3 = lambda * Jp0 * Jr1 + mu * Jp1 * Jr0;
                  float const v4 = lambda * Jp0 * Jr2 + mu * Jp2 * Jr0;
                  float const v5 = lambda * Jp1 * Jr2 + mu * Jp2 * Jr1;
                  float const g0 = grad_u_ref[r][0];
                  float const g1 = grad_u_ref[r][1];
                  float const g2 = grad_u_ref[r][2];
                  flux[p][0] += v0 * g0 + v3 * g1 + v4 * g2;
                  flux[p][1] += v3 * g0 + v1 * g1 + v5 * g2;
                  flux[p][2] += v4 * g0 + v5 * g1 + v2 * g2;
                }
              }
            };
            float cornerCoords[8][3];
            {
              auto const eIdx = mesh_local.elementIndex(elementNumber);
              int I = 0;
              for (int kv = 0; kv < 2; ++kv)
                for (int jv = 0; jv < 2; ++jv)
                  for (int iv = 0; iv < 2; ++iv)
                    mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
            }
            INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork, iso_func);
          }

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

//============================================================================
// computeElementContributions_VTI - VTI
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  auto const struct_metrics_local = m_struct_metrics_;

  Kokkos::parallel_for(
      "Solver Element Contribution Vti",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;
        (void)struct_metrics_local;
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
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * rho_vp2 * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          // Pre-gather Cij per node to avoid scattered global reads inside the kernel.
          // Layout: vti_c[0]=c11, [1]=c12, [2]=c13, [3]=c33, [4]=c44, [5]=c66
          float vti_c[6][kPointsPerElement] = {};
          if constexpr (IS_MODEL_ON_NODES) {
            for (int i = 0; i < dim; ++i)
              for (int j = 0; j < dim; ++j)
                for (int k = 0; k < dim; ++k) {
                  int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                  int const lIdx = i + j * dim + k * dim * dim;
                  float const vp = mesh_local.getModelVpOnNodes(gIdx);
                  float const vs = mesh_local.getModelVsOnNodes(gIdx);
                  float const rho = mesh_local.getModelRhoOnNodes(gIdx);
                  float const delta = mesh_local.getModelDeltaOnNodes(gIdx);
                  float const epsilon = mesh_local.getModelEpsilonOnNodes(gIdx);
                  float const gamma = mesh_local.getModelGammaOnNodes(gIdx);
                  float const rho_vp2 = rho * vp * vp;
                  float const rho_vs2 = rho * vs * vs;
                  float const vp2_vs2 = vp * vp - vs * vs;
                  vti_c[0][lIdx] = rho_vp2 * (1.0f + 2.0f * epsilon);
                  vti_c[3][lIdx] = rho_vp2;
                  vti_c[4][lIdx] = rho_vs2;
                  vti_c[5][lIdx] = rho_vs2 * (1.0f + 2.0f * gamma);
                  vti_c[2][lIdx] =
                      rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * rho_vp2 * delta * vp2_vs2) - rho_vs2;
                  vti_c[1][lIdx] = vti_c[0][lIdx] - 2.0f * vti_c[5][lIdx];
                }
          }

          auto const vti_func = [&](int qa, int qb, int qc, float const (&J_inv)[3][3],
                                    float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
            float c11, c12, c13, c33, c44, c66;
            if constexpr (IS_MODEL_ON_NODES) {
              int const lIdx = qa + qb * dim + qc * dim * dim;
              c11 = vti_c[0][lIdx];
              c12 = vti_c[1][lIdx];
              c13 = vti_c[2][lIdx];
              c33 = vti_c[3][lIdx];
              c44 = vti_c[4][lIdx];
              c66 = vti_c[5][lIdx];
            } else {
              c11 = c11_e;
              c12 = c12_e;
              c13 = c13_e;
              c33 = c33_e;
              c44 = c44_e;
              c66 = c66_e;
            }

            for (int p = 0; p < 3; ++p) {
              float const Jp0 = J_inv[p][0];
              float const Jp1 = J_inv[p][1];
              float const Jp2 = J_inv[p][2];
              flux[p][0] = 0.0f;
              flux[p][1] = 0.0f;
              flux[p][2] = 0.0f;
              for (int r = 0; r < 3; ++r) {
                float const Jr0 = J_inv[r][0];
                float const Jr1 = J_inv[r][1];
                float const Jr2 = J_inv[r][2];
                float const p0r0 = Jp0 * Jr0, p0r1 = Jp0 * Jr1, p0r2 = Jp0 * Jr2;
                float const p1r0 = Jp1 * Jr0, p1r1 = Jp1 * Jr1, p1r2 = Jp1 * Jr2;
                float const p2r0 = Jp2 * Jr0, p2r1 = Jp2 * Jr1, p2r2 = Jp2 * Jr2;
                float const v0 = c11 * p0r0 + c66 * p1r1 + c44 * p2r2;
                float const v1 = c66 * p0r0 + c11 * p1r1 + c44 * p2r2;
                float const v2 = c44 * p0r0 + c44 * p1r1 + c33 * p2r2;
                float const v3 = c66 * p0r1 + c12 * p1r0;
                float const v4 = c44 * p0r2 + c13 * p2r0;
                float const v5 = c44 * p1r2 + c13 * p2r1;
                float const g0 = grad_u_ref[r][0];
                float const g1 = grad_u_ref[r][1];
                float const g2 = grad_u_ref[r][2];
                flux[p][0] += v0 * g0 + v3 * g1 + v4 * g2;
                flux[p][1] += v3 * g0 + v1 * g1 + v5 * g2;
                flux[p][2] += v4 * g0 + v5 * g1 + v2 * g2;
              }
            }
          };

          if constexpr (solver::fe::MeshTypeTraits<MESH_TYPE>::value == utils::enums::meshType::kStruct) {
            float const a = struct_metrics_local.inv_jx;
            float const b = struct_metrics_local.inv_jy;
            float const c = struct_metrics_local.inv_jz;
            float const aa = a * a, bb = b * b, cc = c * c;
            float const ab = a * b, ac = a * c, bc = b * c;
            auto const vti_func_cart = [&](int qa, int qb, int qc, float const (&)[3][3],
                                           float const (&g)[3][3], float (&flux)[3][3]) {
              float c11, c12, c13, c33, c44, c66;
              if constexpr (IS_MODEL_ON_NODES) {
                int const lIdx = qa + qb * dim + qc * dim * dim;
                c11 = vti_c[0][lIdx];
                c12 = vti_c[1][lIdx];
                c13 = vti_c[2][lIdx];
                c33 = vti_c[3][lIdx];
                c44 = vti_c[4][lIdx];
                c66 = vti_c[5][lIdx];
              } else {
                c11 = c11_e;
                c12 = c12_e;
                c13 = c13_e;
                c33 = c33_e;
                c44 = c44_e;
                c66 = c66_e;
              }
              flux[0][0] = c11 * aa * g[0][0] + c12 * ab * g[1][1] + c13 * ac * g[2][2];
              flux[0][1] = c66 * (ab * g[1][0] + aa * g[0][1]);
              flux[0][2] = c44 * (ac * g[2][0] + aa * g[0][2]);
              flux[1][0] = c66 * (bb * g[1][0] + ab * g[0][1]);
              flux[1][1] = c12 * ab * g[0][0] + c11 * bb * g[1][1] + c13 * bc * g[2][2];
              flux[1][2] = c44 * (bc * g[2][1] + bb * g[1][2]);
              flux[2][0] = c44 * (cc * g[2][0] + ac * g[0][2]);
              flux[2][1] = c44 * (cc * g[2][1] + bc * g[1][2]);
              flux[2][2] = c13 * ac * g[0][0] + c13 * bc * g[1][1] + c33 * cc * g[2][2];
            };
            INTEGRAL_TYPE::computeElasticStiffnessSumFactCartesian(struct_metrics_local, localFields, localWork,
                                                                   vti_func_cart);
          } else {
            float cornerCoords[8][3];
            {
              auto const eIdx = mesh_local.elementIndex(elementNumber);
              int I = 0;
              for (int kv = 0; kv < 2; ++kv)
                for (int jv = 0; jv < 2; ++jv)
                  for (int iv = 0; iv < 2; ++iv)
                    mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
            }
            INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork, vti_func);
          }

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

//============================================================================
// computeElementContributions_TTI - TTI
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    auto mesh_local = m_mesh;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    auto const struct_metrics_local = m_struct_metrics_;

    Kokkos::parallel_for(
        "Solver Element Contribution Tti",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
        KOKKOS_LAMBDA(const int _loop_idx) {
          (void)struct_metrics_local;
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

          float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          // Pre-gather TTI params per node to avoid scattered global reads inside the kernel.
          // Layout: tti_p[0]=vp, [1]=vs, [2]=rho, [3]=delta, [4]=epsilon, [5]=gamma, [6]=phi, [7]=theta
          float tti_p[8][kPointsPerElement] = {};
          if constexpr (IS_MODEL_ON_NODES) {
            for (int i = 0; i < dim; ++i)
              for (int j = 0; j < dim; ++j)
                for (int k = 0; k < dim; ++k) {
                  int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                  int const lIdx = i + j * dim + k * dim * dim;
                  tti_p[0][lIdx] = mesh_local.getModelVpOnNodes(gIdx);
                  tti_p[1][lIdx] = mesh_local.getModelVsOnNodes(gIdx);
                  tti_p[2][lIdx] = mesh_local.getModelRhoOnNodes(gIdx);
                  tti_p[3][lIdx] = mesh_local.getModelDeltaOnNodes(gIdx);
                  tti_p[4][lIdx] = mesh_local.getModelEpsilonOnNodes(gIdx);
                  tti_p[5][lIdx] = mesh_local.getModelGammaOnNodes(gIdx);
                  tti_p[6][lIdx] = mesh_local.getModelPhiOnNodes(gIdx);
                  tti_p[7][lIdx] = mesh_local.getModelThetaOnNodes(gIdx);
                }
          }

          auto const tti_func = [&](int qa, int qb, int qc, float const (&J_inv)[3][3],
                                    float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const lIdx = qa + qb * dim + qc * dim * dim;
              for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j)
                  CTTI[i][j] = 0.0f;
              computeCMatrix(tti_p[0][lIdx], tti_p[1][lIdx], tti_p[2][lIdx], tti_p[3][lIdx],
                             tti_p[4][lIdx], tti_p[5][lIdx], tti_p[6][lIdx], tti_p[7][lIdx], CTTI);
            }

            float const C00 = CTTI[0][0], C01 = CTTI[0][1], C02 = CTTI[0][2];
            float const C03 = CTTI[0][3], C04 = CTTI[0][4], C05 = CTTI[0][5];
            float const C11 = CTTI[1][1], C12 = CTTI[1][2], C13 = CTTI[1][3];
            float const C14 = CTTI[1][4], C15 = CTTI[1][5];
            float const C22 = CTTI[2][2], C23 = CTTI[2][3], C24 = CTTI[2][4], C25 = CTTI[2][5];
            float const C33 = CTTI[3][3], C34 = CTTI[3][4], C35 = CTTI[3][5];
            float const C44 = CTTI[4][4], C45 = CTTI[4][5];
            float const C55 = CTTI[5][5];

            for (int p = 0; p < 3; ++p) {
              float const Jp0 = J_inv[p][0];
              float const Jp1 = J_inv[p][1];
              float const Jp2 = J_inv[p][2];
              flux[p][0] = 0.0f;
              flux[p][1] = 0.0f;
              flux[p][2] = 0.0f;
              for (int r = 0; r < 3; ++r) {
                float const Jr0 = J_inv[r][0];
                float const Jr1 = J_inv[r][1];
                float const Jr2 = J_inv[r][2];
                float const p0r0 = Jp0 * Jr0, p0r1 = Jp0 * Jr1, p0r2 = Jp0 * Jr2;
                float const p1r0 = Jp1 * Jr0, p1r1 = Jp1 * Jr1, p1r2 = Jp1 * Jr2;
                float const p2r0 = Jp2 * Jr0, p2r1 = Jp2 * Jr1, p2r2 = Jp2 * Jr2;
                float const v0 = C00 * p0r0 + C05 * p0r1 + C04 * p0r2 + C05 * p1r0 + C55 * p1r1 + C45 * p1r2 +
                                 C04 * p2r0 + C45 * p2r1 + C44 * p2r2;
                float const v1 = C55 * p0r0 + C15 * p0r1 + C35 * p0r2 + C15 * p1r0 + C11 * p1r1 + C13 * p1r2 +
                                 C35 * p2r0 + C13 * p2r1 + C33 * p2r2;
                float const v2 = C44 * p0r0 + C34 * p0r1 + C24 * p0r2 + C34 * p1r0 + C33 * p1r1 + C23 * p1r2 +
                                 C24 * p2r0 + C23 * p2r1 + C22 * p2r2;
                float const v3 = C05 * p0r0 + C01 * p0r1 + C03 * p0r2 + C55 * p1r0 + C15 * p1r1 + C35 * p1r2 +
                                 C45 * p2r0 + C14 * p2r1 + C34 * p2r2;
                float const v4 = C04 * p0r0 + C03 * p0r1 + C02 * p0r2 + C45 * p1r0 + C35 * p1r1 + C25 * p1r2 +
                                 C44 * p2r0 + C34 * p2r1 + C24 * p2r2;
                float const v5 = C45 * p0r0 + C35 * p0r1 + C25 * p0r2 + C14 * p1r0 + C13 * p1r1 + C12 * p1r2 +
                                 C34 * p2r0 + C33 * p2r1 + C23 * p2r2;
                float const g0 = grad_u_ref[r][0];
                float const g1 = grad_u_ref[r][1];
                float const g2 = grad_u_ref[r][2];
                flux[p][0] += v0 * g0 + v3 * g1 + v4 * g2;
                flux[p][1] += v3 * g0 + v1 * g1 + v5 * g2;
                flux[p][2] += v4 * g0 + v5 * g1 + v2 * g2;
              }
            }
          };

          if constexpr (solver::fe::MeshTypeTraits<MESH_TYPE>::value == utils::enums::meshType::kStruct) {
            float const a = struct_metrics_local.inv_jx;
            float const b = struct_metrics_local.inv_jy;
            float const c = struct_metrics_local.inv_jz;
            auto const tti_func_cart = [&](int qa, int qb, int qc, float const (&)[3][3],
                                           float const (&g)[3][3], float (&flux)[3][3]) {
              if constexpr (IS_MODEL_ON_NODES) {
                int const lIdx = qa + qb * dim + qc * dim * dim;
                for (int i = 0; i < 6; ++i)
                  for (int j = 0; j < 6; ++j)
                    CTTI[i][j] = 0.0f;
                computeCMatrix(tti_p[0][lIdx], tti_p[1][lIdx], tti_p[2][lIdx], tti_p[3][lIdx],
                               tti_p[4][lIdx], tti_p[5][lIdx], tti_p[6][lIdx], tti_p[7][lIdx], CTTI);
              }
              float const C00 = CTTI[0][0], C01 = CTTI[0][1], C02 = CTTI[0][2];
              float const C03 = CTTI[0][3], C04 = CTTI[0][4], C05 = CTTI[0][5];
              float const C11 = CTTI[1][1], C12 = CTTI[1][2], C13 = CTTI[1][3];
              float const C14 = CTTI[1][4], C15 = CTTI[1][5];
              float const C22 = CTTI[2][2], C23 = CTTI[2][3], C24 = CTTI[2][4], C25 = CTTI[2][5];
              float const C33 = CTTI[3][3], C34 = CTTI[3][4], C35 = CTTI[3][5];
              float const C44 = CTTI[4][4], C45 = CTTI[4][5];
              float const C55 = CTTI[5][5];
              float const e0 = a * g[0][0];
              float const e1 = b * g[1][1];
              float const e2 = c * g[2][2];
              float const e3 = c * g[2][1] + b * g[1][2];
              float const e4 = c * g[2][0] + a * g[0][2];
              float const e5 = b * g[1][0] + a * g[0][1];
              float const s0 = C00 * e0 + C01 * e1 + C02 * e2 + C03 * e3 + C04 * e4 + C05 * e5;
              float const s1 = C01 * e0 + C11 * e1 + C12 * e2 + C13 * e3 + C14 * e4 + C15 * e5;
              float const s2 = C02 * e0 + C12 * e1 + C22 * e2 + C23 * e3 + C24 * e4 + C25 * e5;
              float const s3 = C03 * e0 + C13 * e1 + C23 * e2 + C33 * e3 + C34 * e4 + C35 * e5;
              float const s4 = C04 * e0 + C14 * e1 + C24 * e2 + C34 * e3 + C44 * e4 + C45 * e5;
              float const s5 = C05 * e0 + C15 * e1 + C25 * e2 + C35 * e3 + C45 * e4 + C55 * e5;
              flux[0][0] = a * s0;
              flux[0][1] = a * s5;
              flux[0][2] = a * s4;
              flux[1][0] = b * s5;
              flux[1][1] = b * s1;
              flux[1][2] = b * s3;
              flux[2][0] = c * s4;
              flux[2][1] = c * s3;
              flux[2][2] = c * s2;
            };
            INTEGRAL_TYPE::computeElasticStiffnessSumFactCartesian(struct_metrics_local, localFields, localWork,
                                                                   tti_func_cart);
          } else {
            float cornerCoords[8][3];
            {
              auto const eIdx = mesh_local.elementIndex(elementNumber);
              int I = 0;
              for (int kv = 0; kv < 2; ++kv)
                for (int jv = 0; jv < 2; ++jv)
                  for (int iv = 0; iv < 2; ++iv)
                    mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
            }
            INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork, tti_func);
          }

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

//============================================================================
// updateFields - Time integration update
//============================================================================
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFields(float dt,
                                                                                          const DataType& data) {
  // Extract scalar constants to local variables
  float const dt_local = dt;
  float const dt2_local = dt * dt;
  int const n_sls = nSls_;
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);

  // Extract single Views and objects to local variables
  auto mesh_local = m_mesh;
  auto mass_matrix = massMatrixGlobal_;
  auto taper_coeff = spongeTaperCoeff_;

  auto sls_w = slsReferenceAngularFrequencies_;
  auto sls_beta = slsAnelasticityCoefficients_;

  std::array<std::remove_reference_t<decltype(data.getCurrentField(0))>, kNumFields> current_field;
  std::array<std::remove_reference_t<decltype(data.getPreviousField(0))>, kNumFields> prev_field;
  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> damping_matrix;
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> work_vector;
  std::array<std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>, kNumFields> atten_work_vec;
  std::array<std::remove_reference_t<decltype(attenuationMemoryVariables_[0])>, kNumFields> atten_mem_vars;

  for (int f = 0; f < kNumFields; ++f) {
    current_field[f] = data.getCurrentField(f);
    prev_field[f] = data.getPreviousField(f);
    damping_matrix[f] = dampingMatrixGlobal_[f];
    work_vector[f] = workVectorsGlobal_[f];
    if (has_attenuation) {
      atten_work_vec[f] = attenuationWorkVectorsGlobal_[f];
      atten_mem_vars[f] = attenuationMemoryVariables_[f];
    }
  }

  bool const list_on = m_node_list_mode_;
  auto list_local = m_node_list_;

  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
    int const n_iter = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();
    Kokkos::parallel_for(
        "Solver Update Field Acoustic", n_iter, KOKKOS_LAMBDA(const int _node_idx) {
          if (_node_idx >= n_iter) return;
          int const I = list_on ? list_local[_node_idx] : _node_idx;
          if (mass_matrix[I] <= 0.0f) return;

          if (mesh_local.isFreeSurface(I)) {
            current_field[0](I) = 0.0f;
            prev_field[0](I) = 0.0f;
          } else {
            float next_val = (2.0f * mass_matrix(I) * current_field[0](I) -
                              (mass_matrix(I) - 0.5f * dt_local * damping_matrix[0](I)) * prev_field[0](I) -
                              dt2_local * work_vector[0](I));

            if (has_attenuation) {
              for (int l = 0; l < n_sls; ++l) {
                float const w = sls_w[l];
                float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                float const gamma_p = 0.5f + 0.5f * gamma;
                float const beta_p = 0.5f * beta;

                next_val += dt2_local * (gamma_p * atten_mem_vars[0](I, l) + beta_p * atten_work_vec[0](I));

                atten_mem_vars[0](I, l) = gamma * atten_mem_vars[0](I, l) + beta * atten_work_vec[0](I);
              }
            }

            prev_field[0](I) = next_val / (mass_matrix(I) + 0.5f * dt_local * damping_matrix[0](I));
            prev_field[0](I) *= taper_coeff(I);
            current_field[0](I) *= taper_coeff(I);
          }
        });
  } else  // ELASTIC
  {
    int const n_iter_el = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();

    Kokkos::parallel_for(
        "Solver Update Field Elastic", n_iter_el, KOKKOS_LAMBDA(const int _node_idx) {
          if (_node_idx >= n_iter_el) return;
          int const I = list_on ? list_local[_node_idx] : _node_idx;
          if (mass_matrix[I] <= 0.0f) return;
          if (mesh_local.isFreeSurface(I)) {
            for (int f = 0; f < kNumFields; ++f) {
              float next_val = (2.0f * mass_matrix(I) * current_field[f](I) - mass_matrix(I) * prev_field[f](I) -
                                dt2_local * work_vector[f](I));

              if (has_attenuation) {
                for (int l = 0; l < n_sls; ++l) {
                  float const w = sls_w[l];
                  float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) + beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) + beta * atten_work_vec[f](I);
                }
              }

              prev_field[f](I) = next_val / mass_matrix(I);
              prev_field[f](I) *= taper_coeff(I);
              current_field[f](I) *= taper_coeff(I);
            }
          } else {
            for (int f = 0; f < kNumFields; ++f) {
              float next_val = (2.0f * mass_matrix(I) * current_field[f](I) -
                                (mass_matrix(I) - 0.5f * dt_local * damping_matrix[f](I)) * prev_field[f](I) -
                                dt2_local * work_vector[f](I));

              if (has_attenuation) {
                for (int l = 0; l < n_sls; ++l) {
                  float const w = sls_w[l];
                  float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) + beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) + beta * atten_work_vec[f](I);
                }
              }

              prev_field[f](I) = next_val / (mass_matrix(I) + 0.5f * dt_local * damping_matrix[f](I));
              prev_field[f](I) *= taper_coeff(I);
              current_field[f](I) *= taper_coeff(I);
            }
          }
        });
  }
}

//============================================================================
// computeGlobalMassMatrix - Assemble mass matrix
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrix() {
  auto mesh_local = m_mesh;
  auto local_massMatrixGlobal = massMatrixGlobal_;

  Kokkos::parallel_for(
      "Solver Compute GMatrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        float massMatrixLocal[kPointsPerElement] = {0};
        int const dim = mesh_local.getOrder() + 1;

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        INTEGRAL_TYPE::computeMassTerm(cornerCoords, [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

        real_t model_factor = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            model_factor =
                1.0f / (mesh_local.getModelVpOnElement(elementNumber) * mesh_local.getModelVpOnElement(elementNumber) *
                        mesh_local.getModelRhoOnElement(elementNumber));
          } else {
            model_factor = mesh_local.getModelRhoOnElement(elementNumber);
          }
        }

        for (int i = 0; i < mesh_local.getNumberOfPointsPerElement(); ++i) {
          int x = i % dim;
          int z = (i / dim) % dim;
          int y = i / (dim * dim);
          int const gIndex = mesh_local.globalNodeIndex(elementNumber, x, y, z);

          if constexpr (IS_MODEL_ON_NODES) {
            if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
              model_factor = 1.0f / (mesh_local.getModelVpOnNodes(gIndex) * mesh_local.getModelVpOnNodes(gIndex) *
                                     mesh_local.getModelRhoOnNodes(gIndex));
            } else {
              model_factor = mesh_local.getModelRhoOnNodes(gIndex);
            }
          }

          massMatrixLocal[i] *= model_factor;
          ATOMICADD(local_massMatrixGlobal[gIndex], massMatrixLocal[i]);
        }
      });
}

//============================================================================
// computeGlobalDampingMatrix - Assemble damping matrix
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrix() {
  auto mesh_local = m_mesh;

  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> local_dampingMatrixGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_dampingMatrixGlobal[f] = dampingMatrixGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Compute Damping Matrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        (void)local_dampingMatrixGlobal;
        for (int i = 0; i < 6; ++i) {
          // Get global face ID for this element face
          int f = mesh_local.getGlobalFace(elementNumber, static_cast<model::CubicFace>(i));

          // Skip internal faces (only process boundary faces)
          if (!mesh_local.isBoundaryFace(f)) continue;

          // Get corner coordinates of the face for integration
          float coords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) {
              coords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
            }
          }

          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            // Acoustic damping
            real_t model_rho = 0.0f;
            real_t model_vp = 0.0f;
            real_t alpha = 0.0f;

            if constexpr (!IS_MODEL_ON_NODES) {
              model_rho = mesh_local.getModelRhoOnElement(elementNumber);
              model_vp = mesh_local.getModelVpOnElement(elementNumber);
              alpha = 1.0 / (model_rho * model_vp);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Skip free surface nodes (no damping on free surface)
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                model_rho = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                model_vp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                alpha = 1.0 / (model_rho * model_vp);
              }

              real_t localIncrement = alpha * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrement);
            }
          } else  // Elastic
          {
            // Elastic damping
            float normal[3];
            mesh_local.faceNormal(elementNumber, static_cast<model::CubicFace>(i), normal);
            real_t nx = normal[0], ny = normal[1], nz = normal[2];

            real_t density, velocityVp, velocityVs;

            if constexpr (!IS_MODEL_ON_NODES) {
              density = mesh_local.getModelRhoOnElement(elementNumber);
              velocityVp = mesh_local.getModelVpOnElement(elementNumber);
              velocityVs = mesh_local.getModelVsOnElement(elementNumber);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Skip free surface nodes (no damping on free surface)
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                density = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                velocityVp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                velocityVs = mesh_local.getModelVsOnNodes(globalNodeIndex);
              }

              real_t aux = density * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              real_t localIncrementx = aux * (velocityVp * fabs(nx) + velocityVs * sqrt(ny * ny + nz * nz));
              real_t localIncrementy = aux * (velocityVp * fabs(ny) + velocityVs * sqrt(nx * nx + nz * nz));
              real_t localIncrementz = aux * (velocityVp * fabs(nz) + velocityVs * sqrt(nx * nx + ny * ny));

              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrementx);
              ATOMICADD(local_dampingMatrixGlobal[1][globalNodeIndex], localIncrementy);
              ATOMICADD(local_dampingMatrixGlobal[2][globalNodeIndex], localIncrementz);
            }
          }
        }
      });
}

//============================================================================
// allocateFEarrays - Allocate memory for FE arrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  massMatrixGlobal_ = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), "massMatrixGlobal");

  static constexpr const char* dampingNames[3] = {"dampingX", "dampingY", "dampingZ"};
  for (int f = 0; f < kNumFields; ++f) {
    dampingMatrixGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), dampingNames[f]);
  }

  // Allocate work vectors for each field
  static constexpr const char* workVectorNames[3] = {"workVec0", "workVec1", "workVec2"};
  for (int f = 0; f < kNumFields; ++f) {
    workVectorsGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), workVectorNames[f]);
  }

  if (attenuationEnabled_ && nSls_ > 0) {
    static constexpr const char* attWorkNames[3] = {"attWorkVec0", "attWorkVec1", "attWorkVec2"};
    static constexpr const char* attMemNames[3] = {"attMemory0", "attMemory1", "attMemory2"};
    for (int f = 0; f < kNumFields; ++f) {
      attenuationWorkVectorsGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), attWorkNames[f]);
      attenuationMemoryVariables_[f] =
          allocateArray2D<ARRAY_REAL_VIEW>(m_mesh.getNumberOfNodes(), nSls_, attMemNames[f]);
    }
  }

  spongeTaperCoeff_ = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), "spongeTaperCoeff");
}

//============================================================================
// initFEarrays - Initialize FE arrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initFEarrays() {
  initSpongeValues();

  if (attenuationEnabled_ && nSls_ > 0) {
    for (int n = 0; n < m_mesh.getNumberOfNodes(); ++n) {
      for (int f = 0; f < kNumFields; ++f) {
        attenuationWorkVectorsGlobal_[f](n) = 0.0f;
        for (int l = 0; l < nSls_; ++l) {
          attenuationMemoryVariables_[f](n, l) = 0.0f;
        }
      }
    }
    FENCE
  }
}

//============================================================================
// initSpongeValues - Initialize absorbing boundary coefficients
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initSpongeValues() {
  const double sigma_max = 0.15;

  for (int n = 0; n < m_mesh.getNumberOfNodes(); n++) {
    const double x = m_mesh.nodeCoord(n, 0);
    const double y = m_mesh.nodeCoord(n, 1);
    const double z = m_mesh.nodeCoord(n, 2);

    const double distToFrontierX = (surface_sponge_) ? m_mesh.domainSize(0) - x : min(m_mesh.domainSize(0) - x, x);
    const double distToFrontierY = min(m_mesh.domainSize(1) - y, y);
    const double distToFrontierZ = min(m_mesh.domainSize(2) - z, z);

    double minDistToFrontier = max(m_mesh.domainSize(0), max(m_mesh.domainSize(1), m_mesh.domainSize(2)));

    bool is_sponge = false;
    if (distToFrontierX < sponge_size_[0]) {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierX);
    }
    if (distToFrontierY < sponge_size_[1]) {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierY);
    }
    if (distToFrontierZ < sponge_size_[2]) {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierZ);
    }

    if (is_sponge) {
      double d = minDistToFrontier;
      double delta = taper_delta_;
      double sigma = sigma_max * std::exp(-((d / delta) * (d / delta)));
      spongeTaperCoeff_(n) = 1.0 / (1.0 + sigma);
    } else {
      spongeTaperCoeff_(n) = 1.0;
    }
  }

  FENCE
}

//============================================================================
// outputSolutionValues - Output field values for diagnostics
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const VECTOR_REAL_VIEW& fieldGlobal, const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << fieldGlobal(m_mesh.globalNodeIndex(e, 0, 0, 0)) << endl;
}

//============================================================================
// computeCMatrix - Compute TTI elasticity tensor (for nodes only)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
template <physicType P, typename>
PROXY_HOST_DEVICE void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeCMatrix(
    float const vp, float const vs, float const rho, float const delta, float const epsilon, float const gamma,
    float const phi, float const theta, float (&CTTI)[6][6]) {
  const float rho_vp2 = rho * vp * vp;
  const float rho_vs2 = rho * vs * vs;
  const float two_eps = 2.0f * epsilon;
  const float two_gam = 2.0f * gamma;

  float CVTI[6][6] = {0.0f};
  CVTI[0][0] = rho_vp2 * (1.0f + two_eps);
  CVTI[1][1] = CVTI[0][0];
  CVTI[2][2] = rho_vp2;
  CVTI[3][3] = rho_vs2;
  CVTI[4][4] = CVTI[3][3];
  CVTI[5][5] = rho_vs2 * (1.0f + two_gam);
  CVTI[0][1] = CVTI[0][0] - 2.0f * CVTI[5][5];
  CVTI[1][0] = CVTI[0][1];

  const float vp2_vs2 = vp * vp - vs * vs;
  const float sqrt_arg = vp2_vs2 * vp2_vs2 + 2.0f * rho_vp2 * delta * vp2_vs2;
  CVTI[0][2] = rho * sqrtf(sqrt_arg) - rho_vs2;
  CVTI[1][2] = CVTI[0][2];
  CVTI[2][0] = CVTI[0][2];
  CVTI[2][1] = CVTI[0][2];

  constexpr float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  const float theta_rad = theta * DEG_TO_RAD;
  const float phi_rad = phi * DEG_TO_RAD;

  const float ctheta = cosf(theta_rad);
  const float stheta = sinf(theta_rad);
  const float cphi = cosf(phi_rad);
  const float sphi = sinf(phi_rad);

  const float ct_cp = ctheta * cphi;
  const float ct_sp = ctheta * sphi;
  const float st_cp = stheta * cphi;
  const float st_sp = stheta * sphi;

  float R[3][3];
  R[0][0] = ct_cp;
  R[0][1] = ct_sp;
  R[0][2] = -stheta;
  R[1][0] = -sphi;
  R[1][1] = cphi;
  R[1][2] = 0.0f;
  R[2][0] = st_cp;
  R[2][1] = st_sp;
  R[2][2] = ctheta;

  const float R00_2 = R[0][0] * R[0][0];
  const float R01_2 = R[0][1] * R[0][1];
  const float R02_2 = R[0][2] * R[0][2];
  const float R10_2 = R[1][0] * R[1][0];
  const float R11_2 = R[1][1] * R[1][1];
  const float R12_2 = R[1][2] * R[1][2];
  const float R20_2 = R[2][0] * R[2][0];
  const float R21_2 = R[2][1] * R[2][1];
  const float R22_2 = R[2][2] * R[2][2];

  float M[6][6] = {0.0f};
  M[0][0] = R00_2;
  M[0][1] = R01_2;
  M[0][2] = R02_2;
  M[1][0] = R10_2;
  M[1][1] = R11_2;
  M[1][2] = R12_2;
  M[2][0] = R20_2;
  M[2][1] = R21_2;
  M[2][2] = R22_2;
  M[0][3] = R[0][1] * R[0][2];
  M[0][4] = R[0][0] * R[0][2];
  M[0][5] = R[0][0] * R[0][1];
  M[1][3] = R[1][1] * R[1][2];
  M[1][4] = R[1][0] * R[1][2];
  M[1][5] = R[1][0] * R[1][1];
  M[2][3] = R[2][1] * R[2][2];
  M[2][4] = R[2][0] * R[2][2];
  M[2][5] = R[2][0] * R[2][1];
  M[3][0] = 2.0f * R[1][0] * R[2][0];
  M[3][1] = 2.0f * R[1][1] * R[2][1];
  M[3][2] = 2.0f * R[1][2] * R[2][2];
  M[3][3] = R[1][1] * R[2][2] + R[1][2] * R[2][1];
  M[3][4] = R[1][0] * R[2][2] + R[1][2] * R[2][0];
  M[3][5] = R[1][0] * R[2][1] + R[1][1] * R[2][0];
  M[4][0] = 2.0f * R[0][0] * R[2][0];
  M[4][1] = 2.0f * R[0][1] * R[2][1];
  M[4][2] = 2.0f * R[0][2] * R[2][2];
  M[4][3] = R[0][1] * R[2][2] + R[0][2] * R[2][1];
  M[4][4] = R[0][0] * R[2][2] + R[0][2] * R[2][0];
  M[4][5] = R[0][0] * R[2][1] + R[0][1] * R[2][0];
  M[5][0] = 2.0f * R[0][0] * R[1][0];
  M[5][1] = 2.0f * R[0][1] * R[1][1];
  M[5][2] = 2.0f * R[0][2] * R[1][2];
  M[5][3] = R[0][1] * R[1][2] + R[0][2] * R[1][1];
  M[5][4] = R[0][0] * R[1][2] + R[0][2] * R[1][0];
  M[5][5] = R[0][0] * R[1][1] + R[0][1] * R[1][0];

  float temp[6][6];

  // M * CVTI
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++) {
        sum += M[i][k] * CVTI[k][j];
      }
      temp[i][j] = sum;
    }
  }

  // temp * M^T
  for (int i = 0; i < 6; i++) {
    for (int j = i; j < 6; j++) {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++) {
        sum += temp[i][k] * M[j][k];
      }
      CTTI[i][j] = sum;
      if (i != j) CTTI[j][i] = sum;
    }
  }
}

//============================================================================
// computeGlobalMassMatrixMasked - domain-masked mass matrix assembly
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrixMasked(
    const VECTOR_INT_VIEW& elem_mask, int active_value) {
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeGlobalMassMatrix();
  m_mask_enabled_ = false;
}

//============================================================================
// computeDampingMatrixMasked - domain-masked damping matrix assembly
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrixMasked(
    const VECTOR_INT_VIEW& elem_mask, int active_value) {
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeDampingMatrix();
  m_mask_enabled_ = false;
}

//============================================================================
// computeElementContributionsMasked - domain-masked stiffness assembly
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsMasked(
    const DataType& data, const VECTOR_INT_VIEW& elem_mask, int active_value) {
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeElementContributions(data);
  m_mask_enabled_ = false;
}

//============================================================================
// computeElementContributionsFromList - stiffness assembly from compact list
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsFromList(
    const DataType& data, const VECTOR_INT_VIEW& elem_list, int n_elems) {
  m_elem_list_ = elem_list;
  m_n_elem_list_ = n_elems;
  m_list_mode_ = true;
  computeElementContributions(data);
  m_list_mode_ = false;
}

//============================================================================
// updateFieldsFromList - Verlet update restricted to a compact node list
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromList(
    float dt, const DataType& data, const VECTOR_INT_VIEW& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFields(dt, data);
  m_node_list_mode_ = false;
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
