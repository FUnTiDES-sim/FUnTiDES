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
#include "sem_solver_elastic.h"

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::computeFEInit(
    model::ModelApi<float, int> &mesh_in,
    const std::array<float, 3> &sponge_size, const bool surface_sponge,
    const float taper_delta)
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::computeOneStep(
    const float &dt, const int &timeSample, SolverBase::DataStruct &data)
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
  ARRAY_REAL_VIEW const &rhsWeights = myData.m_rhsWeights;

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, i2, rhsTermx,rhsTermy,rhsTermz, rhsElement, rhsWeights);
  FENCE
  computeElementContributions(i2, uxnGlobal,uynGlobal,uznGlobal);
  FENCE
  updateDisplacementField(dt, i1, i2, uxnGlobal,uynGlobal,uznGlobal);
  FENCE
}


template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::resetGlobalVectors(
    int numNodes)
{
  LOOPHEAD(numNodes, i)
  {
    uxGlobal[i] = 0;
    uyGlobal[i] = 0;
    uzGlobal[i] = 0;
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::applyRHSTerm(
    int timeSample, float dt, int i2, const ARRAY_REAL_VIEW &rhsTermx,
    const ARRAY_REAL_VIEW &rhsTermy,const ARRAY_REAL_VIEW &rhsTermz,
    const VECTOR_INT_VIEW &rhsElement, 
    const ARRAY_REAL_VIEW &rhsWeights)
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
          float sourcex = rhsTermx(i, timeSample) * rhsWeights(i, localNodeId);
          float sourcey = rhsTermy(i, timeSample) * rhsWeights(i, localNodeId);
          float sourcez = rhsTermz(i, timeSample) * rhsWeights(i, localNodeId);
          uxGlobal(nodeRHS) -= sourcex;
          uyGlobal(nodeRHS) -= sourcey;
          uzGlobal(nodeRHS) -= sourcez;
        }
      }
    }
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::computeElementContributions(
    int i2, const ARRAY_REAL_VIEW &uxnGlobal, const ARRAY_REAL_VIEW &uynGlobal, const ARRAY_REAL_VIEW &uznGlobal)
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
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

  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh.getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh.globalNodeIndex(elementNumber, i, j, k);
        cornerCoords[I][0] = m_mesh.nodeCoord(nodeIdx, 0);
        cornerCoords[I][1] = m_mesh.nodeCoord(nodeIdx, 1);
        cornerCoords[I][2] = m_mesh.nodeCoord(nodeIdx, 2);
        I++;
      }
    }
  }

  float (CTTI)[6][6];

  if constexpr (!IS_MODEL_ON_NODES)
  {
    float const vp = m_mesh.getModelVpOnElement(elementNumber);
    float const vs = m_mesh.getModelVsOnElement(elementNumber);
    float const rho = m_mesh.getModelRhoOnElement(elementNumber);
    float const delta = m_mesh.getModelDeltaOnElement(elementNumber);
    float const epsilon = m_mesh.getModelEpsilonOnElement(elementNumber);
    float const gamma = m_mesh.getModelGammaOnElement(elementNumber);
    float const phi = m_mesh.getModelPhiOnElement(elementNumber);
    float const theta = m_mesh.getModelThetaOnElement(elementNumber);
    computeCMatrix(vp,vs,rho,delta,epsilon,gamma,phi,theta,CTTI);

  }

  INTEGRAL_TYPE::computeStiffNessTermwithJac(cornerCoords, [&] (int qa, int qb, int qc)
  { 
    if constexpr (IS_MODEL_ON_NODES) {
      int const gIndex = m_mesh.globalNodeIndex(elementNumber, qa, qb, qc);
      float const vp = m_mesh.getModelVpOnNodes(gIndex);
      float const vs = m_mesh.getModelVsOnNodes(gIndex);
      float const rho = m_mesh.getModelRhoOnNodes(gIndex);
      float const delta = m_mesh.getModelDeltaOnNodes(gIndex);
      float const epsilon = m_mesh.getModelEpsilonOnNodes(gIndex);
      float const gamma = m_mesh.getModelGammaOnNodes(gIndex);
      float const phi = m_mesh.getModelPhiOnNodes(gIndex);
      float const theta = m_mesh.getModelThetaOnNodes(gIndex);
      computeCMatrix(vp,vs,rho,delta,epsilon,gamma,phi,theta,CTTI);
    } 

  },

  [&] (int i, int j, float val, float const (&J)[3][3], const int p, const int r)
  {
    float Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
    float Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
    
    float C00 = CTTI[0][0], C01 = CTTI[0][1], C02 = CTTI[0][2];
    float C03 = CTTI[0][3], C04 = CTTI[0][4], C05 = CTTI[0][5];
    float C11 = CTTI[1][1], C12 = CTTI[1][2], C13 = CTTI[1][3];
    float C14 = CTTI[1][4], C15 = CTTI[1][5];
    float C22 = CTTI[2][2], C23 = CTTI[2][3], C24 = CTTI[2][4], C25 = CTTI[2][5];
    float C33 = CTTI[3][3], C34 = CTTI[3][4], C35 = CTTI[3][5];
    float C44 = CTTI[4][4], C45 = CTTI[4][5];
    float C55 = CTTI[5][5];
    
    float Jp0Jr0 = Jp0*Jr0, Jp0Jr1 = Jp0*Jr1, Jp0Jr2 = Jp0*Jr2;
    float Jp1Jr0 = Jp1*Jr0, Jp1Jr1 = Jp1*Jr1, Jp1Jr2 = Jp1*Jr2;
    float Jp2Jr0 = Jp2*Jr0, Jp2Jr1 = Jp2*Jr1, Jp2Jr2 = Jp2*Jr2;
    
    float const Rxx = val*(C00*Jp0Jr0 + C05*Jp0Jr1 + C04*Jp0Jr2 +
                           C05*Jp1Jr0 + C55*Jp1Jr1 + C45*Jp1Jr2 +
                           C04*Jp2Jr0 + C45*Jp2Jr1 + C44*Jp2Jr2);
    
    float const Ryy = val*(C55*Jp0Jr0 + C15*Jp0Jr1 + C35*Jp0Jr2 +
                           C15*Jp1Jr0 + C11*Jp1Jr1 + C13*Jp1Jr2 +
                           C35*Jp2Jr0 + C13*Jp2Jr1 + C33*Jp2Jr2);
    
    float const Rzz = val*(C44*Jp0Jr0 + C34*Jp0Jr1 + C24*Jp0Jr2 +
                           C34*Jp1Jr0 + C33*Jp1Jr1 + C23*Jp1Jr2 +
                           C24*Jp2Jr0 + C23*Jp2Jr1 + C22*Jp2Jr2);
    
    float const Rxy = val*(C05*Jp0Jr0 + C01*Jp0Jr1 + C03*Jp0Jr2 +
                           C55*Jp1Jr0 + C15*Jp1Jr1 + C35*Jp1Jr2 +
                           C45*Jp2Jr0 + C14*Jp2Jr1 + C34*Jp2Jr2);
    
    float const Rxz = val*(C04*Jp0Jr0 + C03*Jp0Jr1 + C02*Jp0Jr2 +
                           C45*Jp1Jr0 + C35*Jp1Jr1 + C25*Jp1Jr2 +
                           C44*Jp2Jr0 + C34*Jp2Jr1 + C24*Jp2Jr2);
    
    float const Ryz = val*(C45*Jp0Jr0 + C35*Jp0Jr1 + C25*Jp0Jr2 +
                           C14*Jp1Jr0 + C13*Jp1Jr1 + C12*Jp1Jr2 +
                           C34*Jp2Jr0 + C33*Jp2Jr1 + C23*Jp2Jr2);
    
    float const localIncrementx = Rxx * uxnLocal[j] + Rxy * uynLocal[j] + Rxz * uznLocal[j];
    float const localIncrementy = Rxy * uxnLocal[j] + Ryy * uynLocal[j] + Ryz * uznLocal[j];
    float const localIncrementz = Rxz * uxnLocal[j] + Ryz * uynLocal[j] + Rzz * uznLocal[j];
    
    
    ux[i] += localIncrementx;
    uy[i] += localIncrementy;
    uz[i] += localIncrementz;

  } ); 


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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::updateDisplacementField(
    float dt, int i1, int i2, const ARRAY_REAL_VIEW &uxnGlobal, const ARRAY_REAL_VIEW &uynGlobal, const ARRAY_REAL_VIEW &uznGlobal)
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
    outputSolutionValues(const int &indexTimeStep, int &i1, int &myElementSource,
                         const ARRAY_REAL_VIEW &fieldGlobal,
                         const char* fieldName)  
{
  cout << "TimeStep=" << indexTimeStep
       << ";  " << fieldName << " @ elementSource location " << myElementSource
       << " after computeOneStep = "
       << fieldGlobal(m_mesh.globalNodeIndex(myElementSource, 0, 0, 0), i1)
       << endl;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::initFEarrays()
{
  INTEGRAL_TYPE::init(m_precomputedIntegralData);
  initSpongeValues();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::allocateFEarrays()
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
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::computeGlobalMassMatrix()
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  float massMatrixLocal[nPointsElement] = {0};

  int dim = m_mesh.getOrder() + 1;

  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh.getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh.globalNodeIndex(elementNumber, i, j, k);
        cornerCoords[I][0] = m_mesh.nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh.nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh.nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  real_t density = 0.0f;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    density = m_mesh.getModelRhoOnElement(elementNumber);
  }

  INTEGRAL_TYPE::computeMassTerm(
      cornerCoords,
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



template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
PROXY_HOST_DEVICE
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::computeCMatrix( float const vp,
     float const vs, float const rho, float const delta, float const epsilon, float const gamma, 
     float const phi, float const  theta, float (&CTTI)[6][6]) const
{

  float CVTI[6][6] = {0.0f};
  CVTI[0][0] = rho * vp * vp * (1.0f + 2.0f*epsilon);
  CVTI[1][1] = CVTI[0][0];
  CVTI[2][2] = rho * vp * vp;
  CVTI[3][3] = rho * vs * vs;
  CVTI[4][4] = CVTI[3][3];
  CVTI[5][5] = rho * vs * vs * (1.0f + 2.0f*gamma);
  
  CVTI[0][1] = CVTI[0][0] - 2.0f * CVTI[5][5];
  CVTI[1][0] = CVTI[0][1];
  
  CVTI[0][2] = rho * sqrtf( (vp*vp - vs*vs)*(vp*vp - vs*vs) + 2.0f*vp*vp*delta*(vp*vp - vs*vs) ) - rho * vs*vs;
  CVTI[1][2] = CVTI[0][2];
  CVTI[2][0] = CVTI[0][2];
  CVTI[2][1] = CVTI[0][2];
  
  float R[3][3];
  float ctheta = cosf(theta), stheta = sinf(theta);
  float cphi   = cosf(phi),   sphi   = sinf(phi);
  
  R[0][0] =  ctheta*cphi;  R[0][1] =  ctheta*sphi;  R[0][2] = -stheta;
  R[1][0] = -sphi;         R[1][1] =  cphi;        R[1][2] = 0.0f;
  R[2][0] =  stheta*cphi;  R[2][1] =  stheta*sphi; R[2][2] = ctheta;
  
  float M[6][6] = {0.0f};
  
  M[0][0] = R[0][0]*R[0][0]; M[0][1] = R[0][1]*R[0][1]; M[0][2] = R[0][2]*R[0][2];
  M[1][0] = R[1][0]*R[1][0]; M[1][1] = R[1][1]*R[1][1]; M[1][2] = R[1][2]*R[1][2];
  M[2][0] = R[2][0]*R[2][0]; M[2][1] = R[2][1]*R[2][1]; M[2][2] = R[2][2]*R[2][2];
  
  M[0][3] = R[0][1]*R[0][2]; M[0][4] = R[0][0]*R[0][2]; M[0][5] = R[0][0]*R[0][1];
  M[1][3] = R[1][1]*R[1][2]; M[1][4] = R[1][0]*R[1][2]; M[1][5] = R[1][0]*R[1][1];
  M[2][3] = R[2][1]*R[2][2]; M[2][4] = R[2][0]*R[2][2]; M[2][5] = R[2][0]*R[2][1];
  
  M[3][0] = 2*R[1][0]*R[2][0]; M[3][1] = 2*R[1][1]*R[2][1]; M[3][2] = 2*R[1][2]*R[2][2];
  M[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
  M[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
  M[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];
  
  M[4][0] = 2*R[0][0]*R[2][0]; M[4][1] = 2*R[0][1]*R[2][1]; M[4][2] = 2*R[0][2]*R[2][2];
  M[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
  M[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
  M[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];
  
  M[5][0] = 2*R[0][0]*R[1][0]; M[5][1] = 2*R[0][1]*R[1][1]; M[5][2] = 2*R[0][2]*R[1][2];
  M[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
  M[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
  M[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];
  
  // 4. Multiplication CTTI = M * CVTI * M^T (déroulé)
  float temp[6][6] = {0.0f};
  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
      temp[i][j] = 0.0f;
  
  for(int i=0; i<6; i++)
      for(int k=0; k<6; k++)
          for(int j=0; j<6; j++)
              temp[i][j] += M[i][k] * CVTI[k][j];
  
  for(int i=0; i<6; i++)
      for(int j=i; j<6; j++) // Triangle supérieur
      {
          CTTI[i][j] = 0.0f;
          for(int k=0; k<6; k++)
              CTTI[i][j] += temp[i][k] * M[j][k];
          if(i!=j) CTTI[j][i] = CTTI[i][j]; // Symétrie
      }
  
  
      


}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,IS_MODEL_ON_NODES>::initSpongeValues()
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
