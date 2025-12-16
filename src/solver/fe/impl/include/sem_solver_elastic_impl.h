//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolverElastic.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolverElastic class servers as a base class for the SEM solver
//
//************************************************************************

#include <data_type.h>

#include <array>
#include <cstdlib>

#include "fe/Integrals.hpp"
#include "model_discretization_interface.h"
#include "sem_solver_elastic.h"

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeFEInit(model::ModelApi<float, int> &mesh_in,
                  const std::array<float, 3> &sponge_size,
                  const bool surface_sponge, const float taper_delta)
{
  if (auto *typed_mesh = dynamic_cast<MESH_TYPE *>(&mesh_in))
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeOneStep(const float &dt, const int &timeSample,
                   SolverBase::DataStruct &data)
{
  // Cast to the specific DataStruct type
  auto &myData = dynamic_cast<SEMsolverDataElastic &>(data);

  int const &i1 = myData.m_i1;
  int const &i2 = myData.m_i2;
  ARRAY_REAL_VIEW const &rhsTermx = myData.m_rhsTermx;
  ARRAY_REAL_VIEW const &rhsTermy = myData.m_rhsTermy;
  ARRAY_REAL_VIEW const &rhsTermz = myData.m_rhsTermz;
  ARRAY_REAL_VIEW const &uxnGlobal = myData.m_uxnGlobal;
  ARRAY_REAL_VIEW const &uynGlobal = myData.m_uynGlobal;
  ARRAY_REAL_VIEW const &uznGlobal = myData.m_uznGlobal;
  VECTOR_INT_VIEW const &rhsElement = myData.m_rhsElement;
  ARRAY_REAL_VIEW const &rhsWeightsX = myData.m_rhsWeightsX;
  ARRAY_REAL_VIEW const &rhsWeightsY = myData.m_rhsWeightsY;
  ARRAY_REAL_VIEW const &rhsWeightsZ = myData.m_rhsWeightsZ;

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, i2, rhsTermx, rhsTermy, rhsTermz, rhsElement,
               rhsWeightsX, rhsWeightsY, rhsWeightsZ);
  FENCE
  computeElementContributions(i2, uxnGlobal, uynGlobal, uznGlobal);
  FENCE
  updateDisplacementField(dt, i1, i2, uxnGlobal, uynGlobal, uznGlobal);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                      IS_MODEL_ON_NODES>::resetGlobalVectors(int numNodes)
{
  LOOPHEAD(numNodes, i)
  {
    uxGlobal[i] = 0;
    uyGlobal[i] = 0;
    uzGlobal[i] = 0;
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    applyRHSTerm(int timeSample, float dt, int i2,
                 const ARRAY_REAL_VIEW &rhsTermx,
                 const ARRAY_REAL_VIEW &rhsTermy,
                 const ARRAY_REAL_VIEW &rhsTermz,
                 const VECTOR_INT_VIEW &rhsElement,
                 const ARRAY_REAL_VIEW &rhsWeightsX,
                 const ARRAY_REAL_VIEW &rhsWeightsY,
                 const ARRAY_REAL_VIEW &rhsWeightsZ)
{
  float const dt2 = dt * dt;
  int nb_rhs_element = rhsElement.extent(0);
  LOOPHEAD(nb_rhs_element, i)
  {
    for (int z = 0; z < ORDER + 1; z++)
    {
      for (int y = 0; y < ORDER + 1; y++)
      {
        for (int x = 0; x < ORDER + 1; x++)
        {
          int localNodeId = x + y * (ORDER + 1) + z * (ORDER + 1) * (ORDER + 1);
          int nodeRHS = m_mesh.globalNodeIndex(rhsElement[i], x, y, z);
          float sourcex = rhsTermx(i, timeSample) * rhsWeightsX(i, localNodeId);
          float sourcey = rhsTermy(i, timeSample) * rhsWeightsY(i, localNodeId);
          float sourcez = rhsTermz(i, timeSample) * rhsWeightsZ(i, localNodeId);
          uxGlobal(nodeRHS) -= sourcex;
          uyGlobal(nodeRHS) -= sourcey;
          uzGlobal(nodeRHS) -= sourcez;
        }
      }
    }
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeElementContributions(int i2, const ARRAY_REAL_VIEW &uxnGlobal,
                                const ARRAY_REAL_VIEW &uynGlobal,
                                const ARRAY_REAL_VIEW &uznGlobal)
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  float uxnLocal[nPointsElement] = {0};
  float uynLocal[nPointsElement] = {0};
  float uznLocal[nPointsElement] = {0};
  float ux[nPointsElement] = {0};
  float uy[nPointsElement] = {0};
  float uz[nPointsElement] = {0};

  int dim = m_mesh.getOrder() + 1;
  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const globalIdx = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    uxnLocal[i] = uxnGlobal(globalIdx, i2);
    uynLocal[i] = uynGlobal(globalIdx, i2);
    uznLocal[i] = uznGlobal(globalIdx, i2);
  }

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, m_mesh,
                                                      transformData);

  // Get material properties
  float vp, vs, rho, mu, lambda;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    vp = m_mesh.getModelVpOnElement(elementNumber);
    vs = m_mesh.getModelVsOnElement(elementNumber);
    rho = m_mesh.getModelRhoOnElement(elementNumber);
    mu = rho * vs * vs;
    lambda = rho * (vp * vp - 2.0f * vs * vs);
  }

  // Storage for R_ij tensors at each quadrature point
  float Rxx[3][3], Ryy[3][3], Rzz[3][3];
  float Rxy[3][3], Rxz[3][3], Ryz[3][3];

  INTEGRAL_TYPE::computeStiffNessTermwithJac(
      transformData,
      [&](int qa, int qb, int qc, float const(&J)[3][3]) {
        // Update material properties for model on nodes
        if constexpr (IS_MODEL_ON_NODES)
        {
          int const gIndex = m_mesh.globalNodeIndex(elementNumber, qa, qb, qc);
          vp = m_mesh.getModelVpOnNodes(gIndex);
          vs = m_mesh.getModelVsOnNodes(gIndex);
          rho = m_mesh.getModelRhoOnNodes(gIndex);
          mu = rho * vs * vs;
          lambda = rho * (vp * vp - 2.0f * vs * vs);
        }

        // Compute R_ij for all direction pairs (p,r)
        for (int p = 0; p < 3; ++p)
        {
          for (int r = 0; r < 3; ++r)
          {
            // J[p] = gradient of basis function p in physical space
            // J[p][0] = dN_p/dx, J[p][1] = dN_p/dy, J[p][2] = dN_p/dz
            float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
            float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];

            // Correct isotropic elasticity tensor components (Voigt)
            Rxx[p][r] = (lambda + 2.0f * mu) * Jp0 * Jr0 + 
                        lambda * (Jp1 * Jr1 + Jp2 * Jr2);
            Ryy[p][r] = (lambda + 2.0f * mu) * Jp1 * Jr1 + 
                        lambda * (Jp0 * Jr0 + Jp2 * Jr2);
            Rzz[p][r] = (lambda + 2.0f * mu) * Jp2 * Jr2 + 
                        lambda * (Jp0 * Jr0 + Jp1 * Jr1);
            Rxy[p][r] = mu * (Jp0 * Jr1 + Jp1 * Jr0);
            Rxz[p][r] = mu * (Jp0 * Jr2 + Jp2 * Jr0);
            Ryz[p][r] = mu * (Jp1 * Jr2 + Jp2 * Jr1);
          }
        }
      },
      [&](int i, int j, float val, const int p, const int r) {
        // Apply stiffness matrix
        float const Rxx_ij = val * Rxx[p][r];
        float const Ryy_ij = val * Ryy[p][r];
        float const Rzz_ij = val * Rzz[p][r];
        float const Rxy_ij = val * Rxy[p][r];
        float const Rxz_ij = val * Rxz[p][r];
        float const Ryz_ij = val * Ryz[p][r];

        ux[i] += Rxx_ij * uxnLocal[j] + Rxy_ij * uynLocal[j] + Rxz_ij * uznLocal[j];
        uy[i] += Rxy_ij * uxnLocal[j] + Ryy_ij * uynLocal[j] + Ryz_ij * uznLocal[j];
        uz[i] += Rxz_ij * uxnLocal[j] + Ryz_ij * uynLocal[j] + Rzz_ij * uznLocal[j];
      });

  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    ATOMICADD(uxGlobal[gIndex], ux[i]);
    ATOMICADD(uyGlobal[gIndex], uy[i]);
    ATOMICADD(uzGlobal[gIndex], uz[i]);
  }

  MAINLOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    updateDisplacementField(float dt, int i1, int i2,
                            const ARRAY_REAL_VIEW &uxnGlobal,
                            const ARRAY_REAL_VIEW &uynGlobal,
                            const ARRAY_REAL_VIEW &uznGlobal)
{
  float const dt2 = dt * dt;
  LOOPHEAD(m_mesh.getNumberOfNodes(), I)
  {
    uxnGlobal(I, i1) = 2 * uxnGlobal(I, i2) - uxnGlobal(I, i1) -
                       dt2 * uxGlobal[I] / massMatrixGlobal[I];
    uxnGlobal(I, i1) *= spongeTaperCoeff(I);
    uxnGlobal(I, i2) *= spongeTaperCoeff(I);
    uynGlobal(I, i1) = 2 * uynGlobal(I, i2) - uynGlobal(I, i1) -
                       dt2 * uyGlobal[I] / massMatrixGlobal[I];
    uynGlobal(I, i1) *= spongeTaperCoeff(I);
    uynGlobal(I, i2) *= spongeTaperCoeff(I);
    uznGlobal(I, i1) = 2 * uznGlobal(I, i2) - uznGlobal(I, i1) -
                       dt2 * uzGlobal[I] / massMatrixGlobal[I];
    uznGlobal(I, i1) *= spongeTaperCoeff(I);
    uznGlobal(I, i2) *= spongeTaperCoeff(I);
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    outputSolutionValues(const int &indexTimeStep, int &i1,
                         int &myElementSource,
                         const ARRAY_REAL_VIEW &fieldGlobal,
                         const char *fieldName)
{
  cout << "TimeStep=" << indexTimeStep << ";  " << fieldName
       << " @ elementSource location " << myElementSource
       << " after computeOneStep = "
       << fieldGlobal(m_mesh.globalNodeIndex(myElementSource, 0, 0, 0), i1)
       << endl;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                      IS_MODEL_ON_NODES>::initFEarrays()
{
  initSpongeValues();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                      IS_MODEL_ON_NODES>::allocateFEarrays()
{
  int nbQuadraturePoints = (m_mesh.getOrder() + 1) * (m_mesh.getOrder() + 1) *
                           (m_mesh.getOrder() + 1);

  // shared arrays
  massMatrixGlobal = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(),
                                                      "massMatrixGlobal");
  uxGlobal =
      allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), "uxGlobal");
  uyGlobal =
      allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), "uyGlobal");
  uzGlobal =
      allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), "uzGlobal");

  // sponge allocation
  spongeTaperCoeff = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(),
                                                      "spongeTaperCoeff");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                      IS_MODEL_ON_NODES>::computeGlobalMassMatrix()
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  float massMatrixLocal[nPointsElement] = {0};

  int dim = m_mesh.getOrder() + 1;

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, m_mesh,
                                                      transformData);

  real_t density = 0.0f;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    density = m_mesh.getModelRhoOnElement(elementNumber);
  }

  INTEGRAL_TYPE::computeMassTerm(
      transformData,
      [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    if constexpr (IS_MODEL_ON_NODES)
    {
      density = m_mesh.getModelRhoOnNodes(gIndex);
    }
    massMatrixLocal[i] *= density;
    ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
  }

  MAINLOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
PROXY_HOST_DEVICE void SEMsolverElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::computeCMatrix(float const vp, float const vs,
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
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                      IS_MODEL_ON_NODES>::initSpongeValues()
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

    // Compute taper coefficient using the original Gaussian formula
    if (is_sponge)
    {
      // d = distance from absorption boundary
      double d = minDistToFrontier;
      // δ = characteristic width of the Gaussian
      double delta = taper_delta_;
      // σ(d) = σ_max * exp(-(d/δ)²)
      double sigma = sigma_max * std::exp(-((d / delta) * (d / delta)));
      // Convert to taper coefficient
      spongeTaperCoeff(n) = 1.0 / (1.0 + sigma);
    }
    else
    {
      // No damping in physical domain
      spongeTaperCoeff(n) = 1.0;
    }
  }

  FENCE
}
