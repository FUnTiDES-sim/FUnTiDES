//************************************************************************
//   FUnTiDES - Finite/Unstructured waveform propagation
//
//  sem_solver_impl.h: Unified SEM solver implementation
//
//  Template implementation for the unified SEMsolver class supporting
//  both acoustic and elastic wave propagation.
//************************************************************************

#ifndef SRC_SOLVER_FE_IMPL_INCLUDE_SEMSOLVERIMPL_H_
#define SRC_SOLVER_FE_IMPL_INCLUDE_SEMSOLVERIMPL_H_

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
  computeGlobalMassMatrix();
}

//============================================================================
// computeOneStep - Main time stepping routine
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeOneStep(const float& dt, const int& timeSample,
                                        SolverBase::DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE

  applyRHSTerm(timeSample, dt, myData);
  FENCE

  computeElementContributions(myData);
  FENCE

  updateFields(dt, myData);
  FENCE
}

//============================================================================
// resetGlobalVectors - Zero out work vectors
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::resetGlobalVectors(int numNodes)
{
  LOOPHEAD(numNodes, i)
  {
    for (int f = 0; f < kNumFields; ++f)
    {
      workVectorsGlobal_[f][i] = 0;
    }
  }
  LOOPEND
}

//============================================================================
// applyRHSTerm - Apply source forcing terms
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::applyRHSTerm(int timeSample, float dt,
                                      const DataType& data)
{
  int nb_rhs_element = data.getRhsElement().extent(0);

  LOOPHEAD(nb_rhs_element, i)
  {
    for (int z = 0; z < ORDER + 1; z++)
    {
      for (int y = 0; y < ORDER + 1; y++)
      {
        for (int x = 0; x < ORDER + 1; x++)
        {
          int localNodeId = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
          int nodeRHS =
              m_mesh.globalNodeIndex(data.getRhsElement()[i], x, y, z);

          // Apply RHS for each field component
          for (int f = 0; f < kNumRhs; ++f)
          {
            float source = data.getRhsTerm(f)(i, timeSample) *
                           data.getRhsWeights()(i, localNodeId);
            workVectorsGlobal_[f](nodeRHS) -= source;
          }
        }
      }
    }
  }
  LOOPEND
}

//============================================================================
// computeElementContributions - Assemble element contributions
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeElementContributions(const DataType& data)
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  int const dim = m_mesh.getOrder() + 1;

  // Local field arrays - size depends on physics
  float localFields[kNumFields][kPointsPerElement] = {{0}};
  float localWork[kNumFields][kPointsPerElement] = {{0}};

  // Gather global to local
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        int const globalIdx = m_mesh.globalNodeIndex(elementNumber, i, j, k);
        int const localIdx = i + j * dim + k * dim * dim;

        for (int f = 0; f < kNumFields; ++f)
        {
          localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
        }
      }
    }
  }

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, m_mesh,
                                                      transformData);

  // Physics-specific element contribution computation
  if constexpr (PHYSICS == enums::physicType::kAcoustic)
  {
    // Acoustic: scalar stiffness term with 1/rho coefficient
    real_t inv_density = 0.0f;
    if constexpr (!IS_MODEL_ON_NODES)
    {
      inv_density = 1.0f / m_mesh.getModelRhoOnElement(elementNumber);
    }

    INTEGRAL_TYPE::computeStiffnessTerm(
        transformData,
        [&](const int qa, const int qb, const int qc) {
          if constexpr (IS_MODEL_ON_NODES)
          {
            int const gIndex =
                m_mesh.globalNodeIndex(elementNumber, qa, qb, qc);
            inv_density = 1.0f / m_mesh.getModelRhoOnNodes(gIndex);
          }
        },
        [&](const int i, const int j, const real_t val) {
          float localIncrement = inv_density * val * localFields[0][j];
          localWork[0][i] += localIncrement;
        });
  }
  else if constexpr (PHYSICS == enums::physicType::kElastic)
  {
    // Elastic: tensor stiffness with elasticity matrix
    float CTTI[6][6];
    if constexpr (!IS_MODEL_ON_NODES)
    {
      m_mesh.getCTensorOnElement(elementNumber, CTTI);
    }

    // Packed structure for efficient memory access
#ifdef __CUDACC__
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
          if constexpr (IS_MODEL_ON_NODES)
          {
            int const gIndex =
                m_mesh.globalNodeIndex(elementNumber, qa, qb, qc);
            float const vp = m_mesh.getModelVpOnNodes(gIndex);
            float const vs = m_mesh.getModelVsOnNodes(gIndex);
            float const rho = m_mesh.getModelRhoOnNodes(gIndex);
            float const delta = m_mesh.getModelDeltaOnNodes(gIndex);
            float const epsilon = m_mesh.getModelEpsilonOnNodes(gIndex);
            float const gamma = m_mesh.getModelGammaOnNodes(gIndex);
            float const phi = m_mesh.getModelPhiOnNodes(gIndex);
            float const theta = m_mesh.getModelThetaOnNodes(gIndex);
            computeCMatrix(vp, vs, rho, delta, epsilon, gamma, phi, theta,
                           CTTI);
          }

          // Load CTTI into cache
          float const C00 = CTTI[0][0], C01 = CTTI[0][1], C02 = CTTI[0][2];
          float const C03 = CTTI[0][3], C04 = CTTI[0][4], C05 = CTTI[0][5];
          float const C11 = CTTI[1][1], C12 = CTTI[1][2], C13 = CTTI[1][3];
          float const C14 = CTTI[1][4], C15 = CTTI[1][5];
          float const C22 = CTTI[2][2], C23 = CTTI[2][3], C24 = CTTI[2][4],
                      C25 = CTTI[2][5];
          float const C33 = CTTI[3][3], C34 = CTTI[3][4], C35 = CTTI[3][5];
          float const C44 = CTTI[4][4], C45 = CTTI[4][5];
          float const C55 = CTTI[5][5];

          for (int p = 0; p < 3; ++p)
          {
            const float Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
            for (int r = 0; r < 3; ++r)
            {
              const float Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];

              // Precompute products
              const float p0r0 = Jp0 * Jr0, p0r1 = Jp0 * Jr1, p0r2 = Jp0 * Jr2;
              const float p1r0 = Jp1 * Jr0, p1r1 = Jp1 * Jr1, p1r2 = Jp1 * Jr2;
              const float p2r0 = Jp2 * Jr0, p2r1 = Jp2 * Jr1, p2r2 = Jp2 * Jr2;

              const int idx = p * 3 + r;

              float v0 = C00 * p0r0 + C05 * p0r1 + C04 * p0r2 + C05 * p1r0 +
                         C55 * p1r1 + C45 * p1r2 + C04 * p2r0 + C45 * p2r1 +
                         C44 * p2r2;  // Rxx
              float v1 = C55 * p0r0 + C15 * p0r1 + C35 * p0r2 + C15 * p1r0 +
                         C11 * p1r1 + C13 * p1r2 + C35 * p2r0 + C13 * p2r1 +
                         C33 * p2r2;  // Ryy
              float v2 = C44 * p0r0 + C34 * p0r1 + C24 * p0r2 + C34 * p1r0 +
                         C33 * p1r1 + C23 * p1r2 + C24 * p2r0 + C23 * p2r1 +
                         C22 * p2r2;  // Rzz
              float v3 = C05 * p0r0 + C01 * p0r1 + C03 * p0r2 + C55 * p1r0 +
                         C15 * p1r1 + C35 * p1r2 + C45 * p2r0 + C14 * p2r1 +
                         C34 * p2r2;  // Rxy
              float v4 = C04 * p0r0 + C03 * p0r1 + C02 * p0r2 + C45 * p1r0 +
                         C35 * p1r1 + C25 * p1r2 + C44 * p2r0 + C34 * p2r1 +
                         C24 * p2r2;  // Rxz
              float v5 = C45 * p0r0 + C35 * p0r1 + C25 * p0r2 + C14 * p1r0 +
                         C13 * p1r1 + C12 * p1r2 + C34 * p2r0 + C33 * p2r1 +
                         C23 * p2r2;  // Ryz

#ifdef __CUDACC__
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

#ifdef __CUDACC__
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
  }

  // Scatter local to global
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        int const globalIdx = m_mesh.globalNodeIndex(elementNumber, i, j, k);
        int const localIdx = i + j * dim + k * dim * dim;

        for (int f = 0; f < kNumFields; ++f)
        {
          ATOMICADD(workVectorsGlobal_[f][globalIdx], localWork[f][localIdx]);
        }
      }
    }
  }

  MAINLOOPEND
}

//============================================================================
// updateFields - Time integration update
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::updateFields(float dt, const DataType& data)
{
  float const dt2 = dt * dt;

  LOOPHEAD(m_mesh.getNumberOfNodes(), I)
  {
    for (int f = 0; f < kNumFields; ++f)
    {
      data.getPreviousField(f)(I) =
          2 * data.getCurrentField(f)(I) - data.getPreviousField(f)(I) -
          dt2 * workVectorsGlobal_[f][I] / massMatrixGlobal_[I];
      data.getPreviousField(f)(I) *= spongeTaperCoeff_(I);
      data.getCurrentField(f)(I) *= spongeTaperCoeff_(I);
    }
  }
  LOOPEND
}

//============================================================================
// computeGlobalMassMatrix - Assemble mass matrix
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::computeGlobalMassMatrix()
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  float massMatrixLocal[kPointsPerElement] = {0};

  int const dim = m_mesh.getOrder() + 1;

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, m_mesh,
                                                      transformData);

  INTEGRAL_TYPE::computeMassTerm(
      transformData,
      [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

  // Physics-specific mass matrix scaling
  real_t model_factor = 0.0f;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    if constexpr (PHYSICS == enums::physicType::kAcoustic)
    {
      // Acoustic: 1/(vp^2 * rho)
      model_factor = 1.0f / (m_mesh.getModelVpOnElement(elementNumber) *
                             m_mesh.getModelVpOnElement(elementNumber) *
                             m_mesh.getModelRhoOnElement(elementNumber));
    }
    else
    {
      // Elastic: rho (density)
      model_factor = m_mesh.getModelRhoOnElement(elementNumber);
    }
  }

  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);

    if constexpr (IS_MODEL_ON_NODES)
    {
      if constexpr (PHYSICS == enums::physicType::kAcoustic)
      {
        model_factor = 1.0f / (m_mesh.getModelVpOnNodes(gIndex) *
                               m_mesh.getModelVpOnNodes(gIndex) *
                               m_mesh.getModelRhoOnNodes(gIndex));
      }
      else
      {
        model_factor = m_mesh.getModelRhoOnNodes(gIndex);
      }
    }

    massMatrixLocal[i] *= model_factor;
    ATOMICADD(massMatrixGlobal_[gIndex], massMatrixLocal[i]);
  }

  MAINLOOPEND
}

//============================================================================
// allocateFEarrays - Allocate memory for FE arrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES,
               PHYSICS>::allocateFEarrays()
{
  // Shared arrays
  massMatrixGlobal_ = allocateVector<VECTOR_REAL_VIEW>(
      m_mesh.getNumberOfNodes(), "massMatrixGlobal");

  // Allocate work vectors for each field
  static constexpr const char* workVectorNames[3] = {"workVec0", "workVec1",
                                                     "workVec2"};
  for (int f = 0; f < kNumFields; ++f)
  {
    workVectorsGlobal_[f] = allocateVector<VECTOR_REAL_VIEW>(
        m_mesh.getNumberOfNodes(), workVectorNames[f]);
  }

  // Sponge allocation
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
// computeCMatrix - Compute elasticity tensor (elastic only)
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

#endif  // SRC_SOLVER_FE_IMPL_INCLUDE_SEM_SOLVER_IMPL_H_
