#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#include <data_type.h>

#include <array>
#include <cstdlib>

#include "fe/Integrals.hpp"
#include "model_discretization_interface.h"
#include "sem_solver.h"

namespace solver
{
namespace fe
{

//============================================================================
// computeFEInit - Initialize finite element structures
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeFEInit(model::ModelApi<float, int>& mesh_in,
                                       const std::array<float, 3>& sponge_size,
                                       const bool surface_sponge,
                                       const float taper_delta)
{
  if (auto* typed_mesh = dynamic_cast<MESH_TYPE*>(&mesh_in))
  {
    m_mesh = *typed_mesh;
  }
  else
  {
    throw std::runtime_error("Incompatible mesh type in solver");
  }

  sponge_size_[0] = sponge_size[0];
  sponge_size_[1] = sponge_size[1];
  sponge_size_[2] = sponge_size[2];
  surface_sponge_ = surface_sponge;
  taper_delta_ = taper_delta;

  allocateFEarrays();
  initFEarrays();

  // Compute Local Mass Matrix
  computeGlobalMassMatrix();
  // Compute Local Damping Matrix
  computeDampingMatrix();

  if (attenuationEnabled_ && nSls_ > 0)
  {
    float minQVal = std::numeric_limits<float>::max();
    for (int e = 0; e < m_mesh.getNumberOfElements(); ++e)
    {
      if constexpr (PHYSICS == enums::physicType::kAcoustic)
      {
        float q =
            IS_MODEL_ON_NODES
                ? m_mesh.getModelQpOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                : m_mesh.getModelQpOnElement(e);
        minQVal = std::min(minQVal, q);
      }
      else
      {
        float qp =
            IS_MODEL_ON_NODES
                ? m_mesh.getModelQpOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                : m_mesh.getModelQpOnElement(e);
        float qs =
            IS_MODEL_ON_NODES
                ? m_mesh.getModelQsOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                : m_mesh.getModelQsOnElement(e);
        minQVal = std::min(minQVal, std::min(qp, qs));
      }
    }

    for (int l = 0; l < nSls_; ++l)
    {
      if (slsAnelasticityCoefficients_.extent(0) > 0 &&
          slsAnelasticityCoefficients_[l] < 0.0f)
      {
        slsAnelasticityCoefficients_[l] =
            2.0f * minQVal / (std::max(1.0001f, minQVal) - 1.0f);
      }
    }
  }
}

//============================================================================
// Compute Forces (Phase 1)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeForces(const float& dt, const int& timeSample,
                                       Solver::DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, myData);
  FENCE
  computeElementContributions(myData);
  FENCE
  if (attenuationEnabled_ && nSls_ > 0)
  {
    computeAttenuationContributions(myData);
    FENCE
  }
}

//============================================================================
// Update Solution (Phase 2)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::updateSolution(const float& dt,
                                        Solver::DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);
  updateFields(dt, myData);
  FENCE
}

//============================================================================
// resetGlobalVectors
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::resetGlobalVectors(int numNodes)
{
  Kokkos::parallel_for(
      "Solver Reset GVector", numNodes, KOKKOS_CLASS_LAMBDA(const int i) {
        for (int f = 0; f < kNumFields; ++f)
        {
          workVectorsGlobal_[f][i] = 0;
          if (attenuationEnabled_ && nSls_ > 0)
          {
            attenuationWorkVectorsGlobal_[f][i] = 0;
          }
        }
      });
}

//============================================================================
// applyRHSTerm
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::applyRHSTerm(int timeSample, float dt,
                                      const DataType& data)
{
  int nb_rhs_element = data.getRhsElement().extent(0);

  auto mesh_local = m_mesh;  // Capture mesh for lambda

  Kokkos::parallel_for(
      "Solver Apply RHSTerm", nb_rhs_element, KOKKOS_CLASS_LAMBDA(const int i) {
        for (int z = 0; z < ORDER + 1; z++)
        {
          for (int y = 0; y < ORDER + 1; y++)
          {
            for (int x = 0; x < ORDER + 1; x++)
            {
              int localNodeId =
                  x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
              int nodeRHS =
                  mesh_local.globalNodeIndex(data.getRhsElement()[i], x, y, z);

              for (int f = 0; f < kNumRhs; ++f)
              {
                float source = data.getRhsTerm(f)(i, timeSample) *
                               data.getRhsWeights()(i, localNodeId);
                workVectorsGlobal_[f](nodeRHS) -= source;
              }
            }
          }
        }
      });
}

//============================================================================
// computeElementContributions - DISPATCH
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeElementContributions(const DataType& data)
{
  if constexpr (PHYSICS == enums::physicType::kElastic)
  {
    if (anisotropyType_ == model::AnisotropyType::kIso)
    {
      computeElementContributions_Iso(data);
    }
    else if (anisotropyType_ == model::AnisotropyType::kVTI)
    {
      computeElementContributions_Vti(data);
    }
    else  // kTTI
    {
      computeElementContributions_Tti(data);
    }
  }
  else  // Acoustic - DISPATCH
  {
    auto mesh_local = m_mesh;

    Kokkos::parallel_for(
        "Solver Element Contribution",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,
                                                 LaunchMinBlocksPerSM>>(
            0, mesh_local.getNumberOfElements()),
        KOKKOS_CLASS_LAMBDA(const int elementNumber) {
          if (elementNumber >= mesh_local.getNumberOfElements()) return;

          int const dim = mesh_local.getOrder() + 1;
          float localFields[kNumFields][kPointsPerElement] = {{0}};
          float localWork[kNumFields][kPointsPerElement] = {{0}};

          for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
            {
              for (int k = 0; k < dim; ++k)
              {
                int const globalIdx =
                    mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f)
                {
                  localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
                }
              }
            }
          }

          typename INTEGRAL_TYPE::TransformType transformData;
          model_discretization_interface::gatherTransformData(
              elementNumber, mesh_local, transformData);

          real_t inv_density = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES)
          {
            inv_density = 1.0f / mesh_local.getModelRhoOnElement(elementNumber);
          }

          INTEGRAL_TYPE::computeStiffnessTerm(
              transformData,
              [&](const int qa, const int qb, const int qc) {
                if constexpr (IS_MODEL_ON_NODES)
                {
                  int const gIndex =
                      mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  inv_density = 1.0f / mesh_local.getModelRhoOnNodes(gIndex);
                }
              },
              [&](const int i, const int j, const real_t val) {
                float localIncrement = inv_density * val * localFields[0][j];
                localWork[0][i] += localIncrement;
              });

          for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
            {
              for (int k = 0; k < dim; ++k)
              {
                int const globalIdx =
                    mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f)
                {
                  ATOMICADD(workVectorsGlobal_[f][globalIdx],
                            localWork[f][localIdx]);
                }
              }
            }
          }
        });
  }
}

//============================================================================
// computeAttenuationContributions - Assemble attenuation stiffness terms
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeAttenuationContributions(const DataType& data)
{
  if (!attenuationEnabled_ || nSls_ <= 0) return;

  auto mesh_local = m_mesh;

  Kokkos::parallel_for(
      "Solver Attenuation Contributions",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh_local.getNumberOfElements()) return;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWorkA[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            for (int k = 0; k < dim; ++k)
            {
              int const globalIdx =
                  mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
              {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        typename INTEGRAL_TYPE::TransformType transformData;
        model_discretization_interface::gatherTransformData(
            elementNumber, mesh_local, transformData);

        if constexpr (PHYSICS == enums::physicType::kAcoustic)
        {
          real_t inv_density_q = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES)
          {
            inv_density_q =
                1.0f / (mesh_local.getModelRhoOnElement(elementNumber) *
                        mesh_local.getModelQpOnElement(elementNumber));
          }

          INTEGRAL_TYPE::computeStiffnessTerm(
              transformData,
              [&](const int qa, const int qb, const int qc) {
                if constexpr (IS_MODEL_ON_NODES)
                {
                  int const gIndex =
                      mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  inv_density_q =
                      1.0f / (mesh_local.getModelRhoOnNodes(gIndex) *
                              mesh_local.getModelQpOnNodes(gIndex));
                }
              },
              [&](const int i, const int j, const real_t val) {
                localWorkA[0][i] += inv_density_q * val * localFields[0][j];
              });
        }
        else
        {
          if (anisotropyType_ != model::AnisotropyType::kIso) return;

          struct CJPacked
          {
            float a0, a1, a2, a3;
            float b0, b1;
          };
          CJPacked CJflat[3 * 3];

          INTEGRAL_TYPE::computeStiffNessTermwithJac(
              transformData,
              [&](int qa, int qb, int qc, float const(&J)[3][3]) {
                float vp, vs, rho, qp, qs;
                if constexpr (IS_MODEL_ON_NODES)
                {
                  int const gIndex =
                      mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  vp = mesh_local.getModelVpOnNodes(gIndex);
                  vs = mesh_local.getModelVsOnNodes(gIndex);
                  rho = mesh_local.getModelRhoOnNodes(gIndex);
                  qp = mesh_local.getModelQpOnNodes(gIndex);
                  qs = mesh_local.getModelQsOnNodes(gIndex);
                }
                else
                {
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

                for (int p = 0; p < 3; ++p)
                {
                  float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
                  for (int r = 0; r < 3; ++r)
                  {
                    float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                    int const idx = p * 3 + r;
                    CJflat[idx].a0 = lambda_plus_2mu * Jp0 * Jr0 +
                                     mu * (Jp1 * Jr1 + Jp2 * Jr2);
                    CJflat[idx].a1 = mu * Jp0 * Jr0 +
                                     lambda_plus_2mu * Jp1 * Jr1 +
                                     mu * Jp2 * Jr2;
                    CJflat[idx].a2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) +
                                     lambda_plus_2mu * Jp2 * Jr2;
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
                localWorkA[0][i] +=
                    val * (CJflat[idx].a0 * uxj + CJflat[idx].a3 * uyj +
                           CJflat[idx].b0 * uzj);
                localWorkA[1][i] +=
                    val * (CJflat[idx].a3 * uxj + CJflat[idx].a1 * uyj +
                           CJflat[idx].b1 * uzj);
                localWorkA[2][i] +=
                    val * (CJflat[idx].b0 * uxj + CJflat[idx].b1 * uyj +
                           CJflat[idx].a2 * uzj);
              });
        }

        for (int i = 0; i < dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            for (int k = 0; k < dim; ++k)
            {
              int const globalIdx =
                  mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
              {
                ATOMICADD(attenuationWorkVectorsGlobal_[f][globalIdx],
                          localWorkA[f][localIdx]);
              }
            }
          }
        }
      });
}

//============================================================================
// computeElementContributions_Iso - ISOTROPIC
//============================================================================

//============================================================================
// computeElementContributions_Iso - ISOTROPIC OPTIMIZED
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeElementContributions_Iso(const DataType& data)
{
  auto mesh_local = m_mesh;

  Kokkos::parallel_for(
      "Solver Element Contribution Iso",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh_local.getNumberOfElements()) return;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            for (int k = 0; k < dim; ++k)
            {
              int const globalIdx =
                  mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f)
              {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        typename INTEGRAL_TYPE::TransformType transformData;
        model_discretization_interface::gatherTransformData(
            elementNumber, mesh_local, transformData);

#if defined(__CUDACC__) || defined(__HIPCC__)

        struct CJPacked
        {
          float4 a;
          float2 b;
        };
#else
        struct CJPacked
        {
          alignas(16) float a0, a1, a2, a3;
          float b0, b1;
          float pad[2];
        };
#endif
        CJPacked CJflat[3 * 3];

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            transformData,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              // Get material properties
              float vp, vs, rho;
              if constexpr (IS_MODEL_ON_NODES)
              {
                int const gIndex =
                    mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                vp = mesh_local.getModelVpOnNodes(gIndex);
                vs = mesh_local.getModelVsOnNodes(gIndex);
                rho = mesh_local.getModelRhoOnNodes(gIndex);
              }
              else
              {
                vp = mesh_local.getModelVpOnElement(elementNumber);
                vs = mesh_local.getModelVsOnElement(elementNumber);
                rho = mesh_local.getModelRhoOnElement(elementNumber);
              }

              // Lamé parameters
              float const mu = rho * vs * vs;
              float const lambda = rho * (vp * vp - 2.0f * vs * vs);
              float const lambda_plus_2mu = lambda + 2.0f * mu;

              for (int p = 0; p < 3; ++p)
              {
                float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];

                for (int r = 0; r < 3; ++r)
                {
                  float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                  int const idx = p * 3 + r;

                  float const v0 = lambda_plus_2mu * Jp0 * Jr0 +
                                   mu * (Jp1 * Jr1 + Jp2 * Jr2);
                  float const v1 = mu * Jp0 * Jr0 +
                                   lambda_plus_2mu * Jp1 * Jr1 + mu * Jp2 * Jr2;
                  float const v2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) +
                                   lambda_plus_2mu * Jp2 * Jr2;

                  float const v3 = lambda * Jp0 * Jr1 + mu * Jp1 * Jr0;
                  float const v4 = lambda * Jp0 * Jr2 + mu * Jp2 * Jr0;
                  float const v5 = lambda * Jp1 * Jr2 + mu * Jp2 * Jr1;

#if defined(__CUDACC__) || defined(__HIPCC__)

                  CJflat[idx].a = make_float4(v0, v1, v2, v3);
                  CJflat[idx].b = make_float2(v4, v5);
#else
                  CJflat[idx].a0 = v0;
                  CJflat[idx].a1 = v1;
                  CJflat[idx].a2 = v2;
                  CJflat[idx].a3 = v3;
                  CJflat[idx].b0 = v4;
                  CJflat[idx].b1 = v5;
#endif
                }
              }
            },
            [&](int i, int j, float val, const int p, const int r) {
              int const idx = p * 3 + r;
#if defined(__CUDACC__) || defined(__HIPCC__)

              float3 const u_local = make_float3(
                  localFields[0][j], localFields[1][j], localFields[2][j]);
              float4 const a = CJflat[idx].a;
              float2 const b = CJflat[idx].b;
              localWork[0][i] +=
                  fmaf(val * a.x, u_local.x,
                       fmaf(val * a.w, u_local.y, val * b.x * u_local.z));
              localWork[1][i] +=
                  fmaf(val * a.w, u_local.x,
                       fmaf(val * a.y, u_local.y, val * b.y * u_local.z));
              localWork[2][i] +=
                  fmaf(val * b.x, u_local.x,
                       fmaf(val * b.y, u_local.y, val * a.z * u_local.z));
#else
              float const uxj = localFields[0][j];
              float const uyj = localFields[1][j];
              float const uzj = localFields[2][j];
              float const rxx = CJflat[idx].a0;
              float const ryy = CJflat[idx].a1;
              float const rzz = CJflat[idx].a2;
              float const rxy = CJflat[idx].a3;
              float const rxz = CJflat[idx].b0;
              float const ryz = CJflat[idx].b1;
              localWork[0][i] +=
                  fmaf(val * rxx, uxj, fmaf(val * rxy, uyj, val * rxz * uzj));
              localWork[1][i] +=
                  fmaf(val * rxy, uxj, fmaf(val * ryy, uyj, val * ryz * uzj));
              localWork[2][i] +=
                  fmaf(val * rxz, uxj, fmaf(val * ryz, uyj, val * rzz * uzj));
#endif
            });

        for (int i = 0; i < dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            for (int k = 0; k < dim; ++k)
            {
              int const globalIdx =
                  mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f)
              {
                ATOMICADD(workVectorsGlobal_[f][globalIdx],
                          localWork[f][localIdx]);
              }
            }
          }
        }
      });
}

//============================================================================
// computeElementContributions_VTI - VTI
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeElementContributions_Vti(const DataType& data)
{
  auto mesh_local = m_mesh;

  Kokkos::parallel_for(
      "Solver Element Contribution Vti",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh_local.getNumberOfElements()) return;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            for (int k = 0; k < dim; ++k)
            {
              int const globalIdx =
                  mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f)
              {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        typename INTEGRAL_TYPE::TransformType transformData;
        model_discretization_interface::gatherTransformData(
            elementNumber, mesh_local, transformData);

#if defined(__CUDACC__) || defined(__HIPCC__)

        struct CJPacked
        {
          float4 a;
          float2 b;
        };
#else
        struct CJPacked
        {
          alignas(16) float a0, a1, a2, a3;
          float b0, b1;
          float pad[2];
        };
#endif
        CJPacked CJflat[3 * 3];

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            transformData,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              // Get material properties + Thomsen parameters
              float vp, vs, rho, delta, epsilon, gamma;

              if constexpr (IS_MODEL_ON_NODES)
              {
                int const gIndex =
                    mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                vp = mesh_local.getModelVpOnNodes(gIndex);
                vs = mesh_local.getModelVsOnNodes(gIndex);
                rho = mesh_local.getModelRhoOnNodes(gIndex);
                delta = mesh_local.getModelDeltaOnNodes(gIndex);
                epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
                gamma = mesh_local.getModelGammaOnNodes(gIndex);
              }
              else
              {
                vp = mesh_local.getModelVpOnElement(elementNumber);
                vs = mesh_local.getModelVsOnElement(elementNumber);
                rho = mesh_local.getModelRhoOnElement(elementNumber);
                delta = mesh_local.getModelDeltaOnElement(elementNumber);
                epsilon = mesh_local.getModelEpsilonOnElement(elementNumber);
                gamma = mesh_local.getModelGammaOnElement(elementNumber);
              }

              // Compute 5 independent VTI coefficients
              float const rho_vp2 = rho * vp * vp;
              float const rho_vs2 = rho * vs * vs;
              float const c33 = rho_vp2;
              float const c44 = rho_vs2;
              float const c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
              float const c66 = rho_vs2 * (1.0f + 2.0f * gamma);

              float const vp2_vs2 = vp * vp - vs * vs;
              float const sqrt_arg =
                  vp2_vs2 * vp2_vs2 + 2.0f * rho_vp2 * delta * vp2_vs2;
              float const c13 = rho * sqrtf(sqrt_arg) - rho_vs2;
              float const c12 = c11 - 2.0f * c66;

              for (int p = 0; p < 3; ++p)
              {
                float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];

                for (int r = 0; r < 3; ++r)
                {
                  float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                  float const p0r0 = Jp0 * Jr0, p0r1 = Jp0 * Jr1,
                              p0r2 = Jp0 * Jr2;
                  float const p1r0 = Jp1 * Jr0, p1r1 = Jp1 * Jr1,
                              p1r2 = Jp1 * Jr2;
                  float const p2r0 = Jp2 * Jr0, p2r1 = Jp2 * Jr1,
                              p2r2 = Jp2 * Jr2;
                  int const idx = p * 3 + r;

                  float const v0 = c11 * p0r0 + c66 * p1r1 + c44 * p2r2;
                  float const v1 = c66 * p0r0 + c11 * p1r1 + c44 * p2r2;
                  float const v2 = c44 * p0r0 + c44 * p1r1 + c33 * p2r2;
                  float const v3 = c66 * p0r1 + c12 * p1r0;
                  float const v4 = c44 * p0r2 + c13 * p2r0;
                  float const v5 = c44 * p1r2 + c13 * p2r1;

#if defined(__CUDACC__) || defined(__HIPCC__)

                  CJflat[idx].a = make_float4(v0, v1, v2, v3);
                  CJflat[idx].b = make_float2(v4, v5);
#else
                  CJflat[idx].a0 = v0;
                  CJflat[idx].a1 = v1;
                  CJflat[idx].a2 = v2;
                  CJflat[idx].a3 = v3;
                  CJflat[idx].b0 = v4;
                  CJflat[idx].b1 = v5;
#endif
                }
              }
            },
            [&](int i, int j, float val, const int p, const int r) {
              int const idx = p * 3 + r;
#if defined(__CUDACC__) || defined(__HIPCC__)

              float3 const u_local = make_float3(
                  localFields[0][j], localFields[1][j], localFields[2][j]);
              float4 const a = CJflat[idx].a;
              float2 const b = CJflat[idx].b;
              localWork[0][i] +=
                  fmaf(val * a.x, u_local.x,
                       fmaf(val * a.w, u_local.y, val * b.x * u_local.z));
              localWork[1][i] +=
                  fmaf(val * a.w, u_local.x,
                       fmaf(val * a.y, u_local.y, val * b.y * u_local.z));
              localWork[2][i] +=
                  fmaf(val * b.x, u_local.x,
                       fmaf(val * b.y, u_local.y, val * a.z * u_local.z));
#else
              float const uxj = localFields[0][j];
              float const uyj = localFields[1][j];
              float const uzj = localFields[2][j];
              float const rxx = CJflat[idx].a0;
              float const ryy = CJflat[idx].a1;
              float const rzz = CJflat[idx].a2;
              float const rxy = CJflat[idx].a3;
              float const rxz = CJflat[idx].b0;
              float const ryz = CJflat[idx].b1;
              localWork[0][i] +=
                  fmaf(val * rxx, uxj, fmaf(val * rxy, uyj, val * rxz * uzj));
              localWork[1][i] +=
                  fmaf(val * rxy, uxj, fmaf(val * ryy, uyj, val * ryz * uzj));
              localWork[2][i] +=
                  fmaf(val * rxz, uxj, fmaf(val * ryz, uyj, val * rzz * uzj));
#endif
            });

        for (int i = 0; i < dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            for (int k = 0; k < dim; ++k)
            {
              int const globalIdx =
                  mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f)
              {
                ATOMICADD(workVectorsGlobal_[f][globalIdx],
                          localWork[f][localIdx]);
              }
            }
          }
        }
      });
}
//============================================================================
// computeElementContributions_TTI - TTI (CODE EXISTANT)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeElementContributions_Tti(const DataType& data)
{
  if constexpr (PHYSICS != enums::physicType::kElastic)
  {
  }
  else
  {
    auto mesh_local = m_mesh;

    Kokkos::parallel_for(
        "Solver Element Contribution Tti",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,
                                                 LaunchMinBlocksPerSM>>(
            0, mesh_local.getNumberOfElements()),
        KOKKOS_CLASS_LAMBDA(const int elementNumber) {
          if (elementNumber >= mesh_local.getNumberOfElements()) return;

          int const dim = mesh_local.getOrder() + 1;
          float localFields[kNumFields][kPointsPerElement] = {{0}};
          float localWork[kNumFields][kPointsPerElement] = {{0}};

          for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
            {
              for (int k = 0; k < dim; ++k)
              {
                int const globalIdx =
                    mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f)
                {
                  localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
                }
              }
            }
          }

          typename INTEGRAL_TYPE::TransformType transformData;
          model_discretization_interface::gatherTransformData(
              elementNumber, mesh_local, transformData);

          float CTTI[6][6];

          // For elements, preload tensor
          if constexpr (!IS_MODEL_ON_NODES)
          {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

#if defined(__CUDACC__) || defined(__HIPCC__)

          struct CJPacked
          {
            float4 a;
            float2 b;
          };
#else
          struct CJPacked
          {
            alignas(16) float a0, a1, a2, a3;
            float b0, b1;
            float pad[2];
          };
#endif
          CJPacked CJflat[3 * 3];

          INTEGRAL_TYPE::computeStiffNessTermwithJac(
              transformData,
              [&](int qa, int qb, int qc, float const(&J)[3][3]) {
                // For nodes, compute on-the-fly with rotation
                if constexpr (IS_MODEL_ON_NODES)
                {
                  int const gIndex =
                      mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  float const vp = mesh_local.getModelVpOnNodes(gIndex);
                  float const vs = mesh_local.getModelVsOnNodes(gIndex);
                  float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                  float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
                  float const epsilon =
                      mesh_local.getModelEpsilonOnNodes(gIndex);
                  float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
                  float const phi = mesh_local.getModelPhiOnNodes(gIndex);
                  float const theta = mesh_local.getModelThetaOnNodes(gIndex);
                  computeCMatrix(vp, vs, rho, delta, epsilon, gamma, phi, theta,
                                 CTTI);
                }

                // Apply tensor C (identique aux autres)
                float const C00 = CTTI[0][0], C01 = CTTI[0][1],
                            C02 = CTTI[0][2];
                float const C03 = CTTI[0][3], C04 = CTTI[0][4],
                            C05 = CTTI[0][5];
                float const C11 = CTTI[1][1], C12 = CTTI[1][2],
                            C13 = CTTI[1][3];
                float const C14 = CTTI[1][4], C15 = CTTI[1][5];
                float const C22 = CTTI[2][2], C23 = CTTI[2][3],
                            C24 = CTTI[2][4], C25 = CTTI[2][5];
                float const C33 = CTTI[3][3], C34 = CTTI[3][4],
                            C35 = CTTI[3][5];
                float const C44 = CTTI[4][4], C45 = CTTI[4][5];
                float const C55 = CTTI[5][5];

                for (int p = 0; p < 3; ++p)
                {
                  const float Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
                  for (int r = 0; r < 3; ++r)
                  {
                    const float Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                    const float p0r0 = Jp0 * Jr0, p0r1 = Jp0 * Jr1,
                                p0r2 = Jp0 * Jr2;
                    const float p1r0 = Jp1 * Jr0, p1r1 = Jp1 * Jr1,
                                p1r2 = Jp1 * Jr2;
                    const float p2r0 = Jp2 * Jr0, p2r1 = Jp2 * Jr1,
                                p2r2 = Jp2 * Jr2;
                    const int idx = p * 3 + r;

                    float v0 = C00 * p0r0 + C05 * p0r1 + C04 * p0r2 +
                               C05 * p1r0 + C55 * p1r1 + C45 * p1r2 +
                               C04 * p2r0 + C45 * p2r1 + C44 * p2r2;
                    float v1 = C55 * p0r0 + C15 * p0r1 + C35 * p0r2 +
                               C15 * p1r0 + C11 * p1r1 + C13 * p1r2 +
                               C35 * p2r0 + C13 * p2r1 + C33 * p2r2;
                    float v2 = C44 * p0r0 + C34 * p0r1 + C24 * p0r2 +
                               C34 * p1r0 + C33 * p1r1 + C23 * p1r2 +
                               C24 * p2r0 + C23 * p2r1 + C22 * p2r2;
                    float v3 = C05 * p0r0 + C01 * p0r1 + C03 * p0r2 +
                               C55 * p1r0 + C15 * p1r1 + C35 * p1r2 +
                               C45 * p2r0 + C14 * p2r1 + C34 * p2r2;
                    float v4 = C04 * p0r0 + C03 * p0r1 + C02 * p0r2 +
                               C45 * p1r0 + C35 * p1r1 + C25 * p1r2 +
                               C44 * p2r0 + C34 * p2r1 + C24 * p2r2;
                    float v5 = C45 * p0r0 + C35 * p0r1 + C25 * p0r2 +
                               C14 * p1r0 + C13 * p1r1 + C12 * p1r2 +
                               C34 * p2r0 + C33 * p2r1 + C23 * p2r2;

#if defined(__CUDACC__) || defined(__HIPCC__)

                    CJflat[idx].a = make_float4(v0, v1, v2, v3);
                    CJflat[idx].b = make_float2(v4, v5);
#else
                    CJflat[idx].a0 = v0;
                    CJflat[idx].a1 = v1;
                    CJflat[idx].a2 = v2;
                    CJflat[idx].a3 = v3;
                    CJflat[idx].b0 = v4;
                    CJflat[idx].b1 = v5;
#endif
                  }
                }
              },
              [&](int i, int j, float val, const int p, const int r) {
                const int idx = p * 3 + r;
#if defined(__CUDACC__) || defined(__HIPCC__)

                const float3 u_local = make_float3(
                    localFields[0][j], localFields[1][j], localFields[2][j]);
                const float4 a = CJflat[idx].a;
                const float2 b = CJflat[idx].b;
                localWork[0][i] +=
                    fmaf(val * a.x, u_local.x,
                         fmaf(val * a.w, u_local.y, val * b.x * u_local.z));
                localWork[1][i] +=
                    fmaf(val * a.w, u_local.x,
                         fmaf(val * a.y, u_local.y, val * b.y * u_local.z));
                localWork[2][i] +=
                    fmaf(val * b.x, u_local.x,
                         fmaf(val * b.y, u_local.y, val * a.z * u_local.z));
#else
                const float uxj = localFields[0][j];
                const float uyj = localFields[1][j];
                const float uzj = localFields[2][j];
                const float rxx = CJflat[idx].a0;
                const float ryy = CJflat[idx].a1;
                const float rzz = CJflat[idx].a2;
                const float rxy = CJflat[idx].a3;
                const float rxz = CJflat[idx].b0;
                const float ryz = CJflat[idx].b1;
                localWork[0][i] +=
                    fmaf(val * rxx, uxj, fmaf(val * rxy, uyj, val * rxz * uzj));
                localWork[1][i] +=
                    fmaf(val * rxy, uxj, fmaf(val * ryy, uyj, val * ryz * uzj));
                localWork[2][i] +=
                    fmaf(val * rxz, uxj, fmaf(val * ryz, uyj, val * rzz * uzj));
#endif
              });

          for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
            {
              for (int k = 0; k < dim; ++k)
              {
                int const globalIdx =
                    mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f)
                {
                  ATOMICADD(workVectorsGlobal_[f][globalIdx],
                            localWork[f][localIdx]);
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
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::updateFields(float dt, const DataType& data)
{
  // Extract scalar constants to local variables
  float const dt_local = dt;
  float const dt2_local = dt * dt;
  int const n_sls = nSls_;
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);

  // Extract single Views and objects to local variables
  auto mesh_local = m_mesh;
  auto mass_matrix = massMatrixGlobal_;
  auto taper_coeff = spongeTaperCoeff_;

  // (Assuming these two are Kokkos::Views and NOT std::vectors.
  // If they are std::vectors, you must copy them to device Views first).
  auto sls_w = slsReferenceAngularFrequencies_;
  auto sls_beta = slsAnelasticityCoefficients_;

  // Extract std::arrays of Views to local arrays so the lambda can capture
  // them by value. We use std::remove_reference_t to ensure we are creating
  // arrays of Views, not arrays of references!
  std::array<std::remove_reference_t<decltype(data.getCurrentField(0))>,
             kNumFields>
      current_field;
  std::array<std::remove_reference_t<decltype(data.getPreviousField(0))>,
             kNumFields>
      prev_field;
  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>,
             kNumFields>
      damping_matrix;
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>,
             kNumFields>
      work_vector;
  std::array<
      std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>,
      kNumFields>
      atten_work_vec;
  std::array<std::remove_reference_t<decltype(attenuationMemoryVariables_[0])>,
             kNumFields>
      atten_mem_vars;

  for (int f = 0; f < kNumFields; ++f)
  {
    current_field[f] = data.getCurrentField(f);
    prev_field[f] = data.getPreviousField(f);
    damping_matrix[f] = dampingMatrixGlobal_[f];
    work_vector[f] = workVectorsGlobal_[f];
    if (has_attenuation)
    {
      atten_work_vec[f] = attenuationWorkVectorsGlobal_[f];
      atten_mem_vars[f] = attenuationMemoryVariables_[f];
    }
  }

  if constexpr (PHYSICS == enums::physicType::kAcoustic)
  {
    Kokkos::parallel_for(
        "Solver Update Field Acoustic", mesh_local.getNumberOfNodes(),
        // Use standard KOKKOS_LAMBDA to capture local variables by value
        KOKKOS_LAMBDA(const int I) {
          if (mesh_local.isFreeSurface(I))
          {
            current_field[0](I) = 0.0f;
            prev_field[0](I) = 0.0f;
          }
          else
          {
            float next_val =
                (2.0f * mass_matrix(I) * current_field[0](I) -
                 (mass_matrix(I) - 0.5f * dt_local * damping_matrix[0](I)) *
                     prev_field[0](I) -
                 dt2_local * work_vector[0](I));

            if (has_attenuation)
            {
              for (int l = 0; l < n_sls; ++l)
              {
                float const w = sls_w[l];
                float const gamma =
                    (2.0f - w * dt_local) / (2.0f + w * dt_local);
                float const beta =
                    sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                float const gamma_p = 0.5f + 0.5f * gamma;
                float const beta_p = 0.5f * beta;

                next_val += dt2_local * (gamma_p * atten_mem_vars[0](I, l) +
                                         beta_p * atten_work_vec[0](I));

                atten_mem_vars[0](I, l) = gamma * atten_mem_vars[0](I, l) +
                                          beta * atten_work_vec[0](I);
              }
            }

            prev_field[0](I) =
                next_val /
                (mass_matrix(I) + 0.5f * dt_local * damping_matrix[0](I));
            prev_field[0](I) *= taper_coeff(I);
            current_field[0](I) *= taper_coeff(I);
          }
        });
  }
  else  // ELASTIC
  {
    Kokkos::parallel_for(
        "Solver Update Field Elastic", mesh_local.getNumberOfNodes(),
        KOKKOS_LAMBDA(const int I) {
          if (mesh_local.isFreeSurface(I))
          {
            for (int f = 0; f < kNumFields; ++f)
            {
              float next_val = (2.0f * mass_matrix(I) * current_field[f](I) -
                                mass_matrix(I) * prev_field[f](I) -
                                dt2_local * work_vector[f](I));

              if (has_attenuation)
              {
                for (int l = 0; l < n_sls; ++l)
                {
                  float const w = sls_w[l];
                  float const gamma =
                      (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta =
                      sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) +
                                           beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) +
                                            beta * atten_work_vec[f](I);
                }
              }

              prev_field[f](I) = next_val / mass_matrix(I);
              prev_field[f](I) *= taper_coeff(I);
              current_field[f](I) *= taper_coeff(I);
            }
          }
          else
          {
            for (int f = 0; f < kNumFields; ++f)
            {
              float next_val =
                  (2.0f * mass_matrix(I) * current_field[f](I) -
                   (mass_matrix(I) - 0.5f * dt_local * damping_matrix[f](I)) *
                       prev_field[f](I) -
                   dt2_local * work_vector[f](I));

              if (has_attenuation)
              {
                for (int l = 0; l < n_sls; ++l)
                {
                  float const w = sls_w[l];
                  float const gamma =
                      (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta =
                      sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) +
                                           beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) +
                                            beta * atten_work_vec[f](I);
                }
              }

              prev_field[f](I) =
                  next_val /
                  (mass_matrix(I) + 0.5f * dt_local * damping_matrix[f](I));
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeGlobalMassMatrix()
{
  auto mesh_local = m_mesh;

  Kokkos::parallel_for(
      "Solver Compute GMatrix",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh_local.getNumberOfElements()) return;

        float massMatrixLocal[kPointsPerElement] = {0};
        int const dim = mesh_local.getOrder() + 1;

        typename INTEGRAL_TYPE::TransformType transformData;
        model_discretization_interface::gatherTransformData(
            elementNumber, mesh_local, transformData);

        INTEGRAL_TYPE::computeMassTerm(
            transformData,
            [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

        real_t model_factor = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES)
        {
          if constexpr (PHYSICS == enums::physicType::kAcoustic)
          {
            model_factor =
                1.0f / (mesh_local.getModelVpOnElement(elementNumber) *
                        mesh_local.getModelVpOnElement(elementNumber) *
                        mesh_local.getModelRhoOnElement(elementNumber));
          }
          else
          {
            model_factor = mesh_local.getModelRhoOnElement(elementNumber);
          }
        }

        for (int i = 0; i < mesh_local.getNumberOfPointsPerElement(); ++i)
        {
          int x = i % dim;
          int z = (i / dim) % dim;
          int y = i / (dim * dim);
          int const gIndex = mesh_local.globalNodeIndex(elementNumber, x, y, z);

          if constexpr (IS_MODEL_ON_NODES)
          {
            if constexpr (PHYSICS == enums::physicType::kAcoustic)
            {
              model_factor = 1.0f / (mesh_local.getModelVpOnNodes(gIndex) *
                                     mesh_local.getModelVpOnNodes(gIndex) *
                                     mesh_local.getModelRhoOnNodes(gIndex));
            }
            else
            {
              model_factor = mesh_local.getModelRhoOnNodes(gIndex);
            }
          }

          massMatrixLocal[i] *= model_factor;
          ATOMICADD(massMatrixGlobal_[gIndex], massMatrixLocal[i]);
        }
      });
}

//============================================================================
// computeGlobalDampingMatrix - Assemble damping matrix
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeDampingMatrix()
{
  auto mesh_local = m_mesh;

  Kokkos::parallel_for(
      "Solver Compute Damping Matrix",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh_local.getNumberOfElements()) return;

        for (int i = 0; i < 6; ++i)
        {
          // Get global face ID for this element face
          int f = mesh_local.getGlobalFace(elementNumber,
                                           static_cast<model::CubicFace>(i));

          // Skip internal faces (only process boundary faces)
          if (!mesh_local.isBoundaryFace(f)) continue;

          // Get corner coordinates of the face for integration
          float coords[4][3];
          for (int j = 0; j < 4; ++j)
          {
            int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(
                f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d)
            {
              coords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
            }
          }

          if constexpr (PHYSICS == enums::physicType::kAcoustic)
          {
            // Acoustic damping
            real_t model_rho = 0.0f;
            real_t model_vp = 0.0f;
            real_t alpha = 0.0f;

            if constexpr (!IS_MODEL_ON_NODES)
            {
              model_rho = mesh_local.getModelRhoOnElement(elementNumber);
              model_vp = mesh_local.getModelVpOnElement(elementNumber);
              alpha = 1.0 / (model_rho * model_vp);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q)
            {
              int const globalNodeIndex =
                  mesh_local.getGlobalNodeFromFace(f, q);

              // Skip free surface nodes (no damping on free surface)
              if (mesh_local.isFreeSurface(globalNodeIndex))
              {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES)
              {
                model_rho = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                model_vp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                alpha = 1.0 / (model_rho * model_vp);
              }

              real_t localIncrement =
                  alpha * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              ATOMICADD(dampingMatrixGlobal_[0][globalNodeIndex],
                        localIncrement);
            }
          }
          else  // Elastic
          {
            // Elastic damping
            float normal[3];
            mesh_local.faceNormal(elementNumber,
                                  static_cast<model::CubicFace>(i), normal);
            real_t nx = normal[0], ny = normal[1], nz = normal[2];

            real_t density, velocityVp, velocityVs;

            if constexpr (!IS_MODEL_ON_NODES)
            {
              density = mesh_local.getModelRhoOnElement(elementNumber);
              velocityVp = mesh_local.getModelVpOnElement(elementNumber);
              velocityVs = mesh_local.getModelVsOnElement(elementNumber);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q)
            {
              int const globalNodeIndex =
                  mesh_local.getGlobalNodeFromFace(f, q);

              // Skip free surface nodes (no damping on free surface)
              if (mesh_local.isFreeSurface(globalNodeIndex))
              {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES)
              {
                density = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                velocityVp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                velocityVs = mesh_local.getModelVsOnNodes(globalNodeIndex);
              }

              real_t aux =
                  density * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              real_t localIncrementx =
                  aux * (velocityVp * fabs(nx) +
                         velocityVs * sqrt(ny * ny + nz * nz));
              real_t localIncrementy =
                  aux * (velocityVp * fabs(ny) +
                         velocityVs * sqrt(nx * nx + nz * nz));
              real_t localIncrementz =
                  aux * (velocityVp * fabs(nz) +
                         velocityVs * sqrt(nx * nx + ny * ny));

              ATOMICADD(dampingMatrixGlobal_[0][globalNodeIndex],
                        localIncrementx);
              ATOMICADD(dampingMatrixGlobal_[1][globalNodeIndex],
                        localIncrementy);
              ATOMICADD(dampingMatrixGlobal_[2][globalNodeIndex],
                        localIncrementz);
            }
          }
        }
      });
}

//============================================================================
// allocateFEarrays - Allocate memory for FE arrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::allocateFEarrays()
{
  massMatrixGlobal_ = allocateVector<VECTOR_REAL_VIEW>(
      m_mesh.getNumberOfNodes(), "massMatrixGlobal");

  static constexpr const char* dampingNames[3] = {"dampingX", "dampingY",
                                                  "dampingZ"};
  for (int f = 0; f < kNumFields; ++f)
  {
    dampingMatrixGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(
        m_mesh.getNumberOfNodes(), dampingNames[f]);
  }

  // Allocate work vectors for each field
  static constexpr const char* workVectorNames[3] = {"workVec0", "workVec1",
                                                     "workVec2"};
  for (int f = 0; f < kNumFields; ++f)
  {
    workVectorsGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(
        m_mesh.getNumberOfNodes(), workVectorNames[f]);
  }

  if (attenuationEnabled_ && nSls_ > 0)
  {
    static constexpr const char* attWorkNames[3] = {
        "attWorkVec0", "attWorkVec1", "attWorkVec2"};
    static constexpr const char* attMemNames[3] = {"attMemory0", "attMemory1",
                                                   "attMemory2"};
    for (int f = 0; f < kNumFields; ++f)
    {
      attenuationWorkVectorsGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(
          m_mesh.getNumberOfNodes(), attWorkNames[f]);
      attenuationMemoryVariables_[f] = allocateArray2D<ARRAY_REAL_VIEW>(
          m_mesh.getNumberOfNodes(), nSls_, attMemNames[f]);
    }
  }

  spongeTaperCoeff_ = allocateVector<VECTOR_REAL_VIEW>(
      m_mesh.getNumberOfNodes(), "spongeTaperCoeff");
}

//============================================================================
// initFEarrays - Initialize FE arrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::initFEarrays()
{
  initSpongeValues();

  if (attenuationEnabled_ && nSls_ > 0)
  {
    for (int n = 0; n < m_mesh.getNumberOfNodes(); ++n)
    {
      for (int f = 0; f < kNumFields; ++f)
      {
        attenuationWorkVectorsGlobal_[f](n) = 0.0f;
        for (int l = 0; l < nSls_; ++l)
        {
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::initSpongeValues()
{
  const double sigma_max = 0.15;

  for (int n = 0; n < m_mesh.getNumberOfNodes(); n++)
  {
    const double x = m_mesh.nodeCoord(n, 0);
    const double y = m_mesh.nodeCoord(n, 1);
    const double z = m_mesh.nodeCoord(n, 2);

    const double distToFrontierX = (surface_sponge_)
                                       ? m_mesh.domainSize(0) - x
                                       : min(m_mesh.domainSize(0) - x, x);
    const double distToFrontierY = min(m_mesh.domainSize(1) - y, y);
    const double distToFrontierZ = min(m_mesh.domainSize(2) - z, z);

    double minDistToFrontier = max(
        m_mesh.domainSize(0), max(m_mesh.domainSize(1), m_mesh.domainSize(2)));

    bool is_sponge = false;
    if (distToFrontierX < sponge_size_[0])
    {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierX);
    }
    if (distToFrontierY < sponge_size_[1])
    {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierY);
    }
    if (distToFrontierZ < sponge_size_[2])
    {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierZ);
    }

    if (is_sponge)
    {
      double d = minDistToFrontier;
      double delta = taper_delta_;
      double sigma = sigma_max * std::exp(-((d / delta) * (d / delta)));
      spongeTaperCoeff_(n) = 1.0 / (1.0 + sigma);
    }
    else
    {
      spongeTaperCoeff_(n) = 1.0;
    }
  }

  FENCE
}

//============================================================================
// outputSolutionValues - Output field values for diagnostics
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::outputSolutionValues(const int& t, int& e,
                                              const VECTOR_REAL_VIEW&
                                                  fieldGlobal,
                                              const char* fieldName)
{
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location "
       << e << " after computeOneStep = "
       << fieldGlobal(m_mesh.globalNodeIndex(e, 0, 0, 0)) << endl;
}

//============================================================================
// computeCMatrix - Compute TTI elasticity tensor (for nodes only)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
template <physicType P, typename>
PROXY_HOST_DEVICE void
SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
          PHYSICS>::computeCMatrix(float const vp, float const vs,
                                   float const rho, float const delta,
                                   float const epsilon, float const gamma,
                                   float const phi, float const theta,
                                   float (&CTTI)[6][6]) const
{
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
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++)
      {
        sum += M[i][k] * CVTI[k][j];
      }
      temp[i][j] = sum;
    }
  }

  // temp * M^T
  for (int i = 0; i < 6; i++)
  {
    for (int j = i; j < 6; j++)
    {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++)
      {
        sum += temp[i][k] * M[j][k];
      }
      CTTI[i][j] = sum;
      if (i != j) CTTI[j][i] = sum;
    }
  }
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
