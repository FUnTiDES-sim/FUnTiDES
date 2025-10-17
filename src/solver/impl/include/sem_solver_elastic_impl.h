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
        cornerCoords[I][2] = m_mesh.nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh.nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  float CTTInodes[nPointsElement][6][6];
  float CTTIelement[6][6];
  //Compute elasticity matrix
  if constexpr (IS_MODEL_ON_NODES)
  {
    for (int i = 0; i < nPointsElement; i++)
    {
      int x = i % dim;
      int z = (i / dim) % dim;
      int y = i / (dim * dim);
      int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
      float const vp = m_mesh.getModelVpOnNodes(gIndex);
      float const vs = m_mesh.getModelVsOnNodes(gIndex);
      float const rho = m_mesh.getModelRhoOnNodes(gIndex);
      float const delta = m_mesh.getModelDeltaOnNodes(gIndex);
      float const epsilon = m_mesh.getModelEpsilonOnNodes(gIndex);
      float const gamma = m_mesh.getModelGammaOnNodes(gIndex);
      float const phi = m_mesh.getModelPhiOnNodes(gIndex);
      float const theta = m_mesh.getModelThetaOnNodes(gIndex);
      computeCMatrix(vp,vs,rho,delta,epsilon,gamma,phi,theta,CTTInodes[i]);
    }
  }

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
    computeCMatrix(vp,vs,rho,delta,epsilon,gamma,phi,theta,CTTIelement);

  }

  float CTTItemp[6][6]; 
  float (*CTTI)[6][6];

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
      computeCMatrix(vp,vs,rho,delta,epsilon,gamma,phi,theta,CTTItemp);
      CTTI = &CTTItemp;
    } else {
      CTTI = &CTTIelement;
    }

  },

  [&] (int i, int j, float val, float const (&J)[3][3], const int p, const int r)
  {
          float const Rxx_ij = val*((*CTTI)[0][0]*J[p][0]*J[r][0]+(*CTTI)[5][0]*J[p][0]*J[r][1]+(*CTTI)[4][0]*J[p][0]*J[r][2]+
                                 (*CTTI)[0][5]*J[p][1]*J[r][0]+(*CTTI)[5][5]*J[p][1]*J[r][1]+(*CTTI)[4][5]*J[p][1]*J[r][2]+
                                 (*CTTI)[0][4]*J[p][2]*J[r][0]+(*CTTI)[5][4]*J[p][2]*J[r][1]+(*CTTI)[4][4]*J[p][2]*J[r][2]);

      float const Ryy_ij = val*((*CTTI)[5][5]*J[p][0]*J[r][0]+(*CTTI)[1][5]*J[p][0]*J[r][1]+(*CTTI)[3][5]*J[p][0]*J[r][2]+
                                 (*CTTI)[5][2]*J[p][1]*J[r][0]+(*CTTI)[2][2]*J[p][1]*J[r][1]+(*CTTI)[3][2]*J[p][1]*J[r][2]+
                                 (*CTTI)[5][3]*J[p][2]*J[r][0]+(*CTTI)[2][3]*J[p][2]*J[r][1]+(*CTTI)[3][3]*J[p][2]*J[r][2]);

      float const Rzz_ij = val*((*CTTI)[4][4]*J[p][0]*J[r][0]+(*CTTI)[3][4]*J[p][0]*J[r][1]+(*CTTI)[2][4]*J[p][0]*J[r][2]+
                                 (*CTTI)[4][3]*J[p][1]*J[r][0]+(*CTTI)[3][3]*J[p][1]*J[r][1]+(*CTTI)[2][3]*J[p][1]*J[r][2]+
                                 (*CTTI)[4][2]*J[p][2]*J[r][0]+(*CTTI)[3][2]*J[p][2]*J[r][1]+(*CTTI)[2][2]*J[p][2]*J[r][2]);

      float const Ryx_ij = val*((*CTTI)[0][5]*J[p][0]*J[r][0]+(*CTTI)[5][5]*J[p][0]*J[r][1]+(*CTTI)[4][5]*J[p][0]*J[r][2]+
                                 (*CTTI)[0][1]*J[p][1]*J[r][0]+(*CTTI)[4][1]*J[p][1]*J[r][1]+(*CTTI)[4][1]*J[p][1]*J[r][2]+
                                 (*CTTI)[0][3]*J[p][2]*J[r][0]+(*CTTI)[5][3]*J[p][2]*J[r][1]+(*CTTI)[4][3]*J[p][2]*J[r][2]);

      float const Rxy_ij = val*((*CTTI)[5][0]*J[p][0]*J[r][0]+(*CTTI)[1][0]*J[p][0]*J[r][1]+(*CTTI)[3][0]*J[p][0]*J[r][2]+
                                 (*CTTI)[5][5]*J[p][1]*J[r][0]+(*CTTI)[1][5]*J[p][1]*J[r][1]+(*CTTI)[3][5]*J[p][1]*J[r][2]+
                                 (*CTTI)[5][4]*J[p][2]*J[r][0]+(*CTTI)[1][4]*J[p][2]*J[r][1]+(*CTTI)[3][4]*J[p][2]*J[r][2]);

      float const Rzx_ij = val*((*CTTI)[0][4]*J[p][0]*J[r][0]+(*CTTI)[5][4]*J[p][0]*J[r][1]+(*CTTI)[4][4]*J[p][0]*J[r][2]+
                                 (*CTTI)[0][3]*J[p][1]*J[r][0]+(*CTTI)[5][3]*J[p][1]*J[r][1]+(*CTTI)[4][3]*J[p][1]*J[r][2]+
                                 (*CTTI)[0][2]*J[p][2]*J[r][0]+(*CTTI)[4][2]*J[p][2]*J[r][1]+(*CTTI)[4][2]*J[p][2]*J[r][2]);

      float const Rxz_ij = val*((*CTTI)[4][0]*J[p][0]*J[r][0]+(*CTTI)[3][0]*J[p][0]*J[r][1]+(*CTTI)[2][0]*J[p][0]*J[r][2]+
                                 (*CTTI)[4][5]*J[p][1]*J[r][0]+(*CTTI)[3][5]*J[p][1]*J[r][1]+(*CTTI)[2][5]*J[p][1]*J[r][2]+
                                 (*CTTI)[4][4]*J[p][2]*J[r][0]+(*CTTI)[3][4]*J[p][2]*J[r][1]+(*CTTI)[2][4]*J[p][2]*J[r][2]);

      float const Rzy_ij = val*((*CTTI)[5][4]*J[p][0]*J[r][0]+(*CTTI)[1][4]*J[p][0]*J[r][1]+(*CTTI)[3][4]*J[p][0]*J[r][2]+
                                 (*CTTI)[5][3]*J[p][1]*J[r][0]+(*CTTI)[1][3]*J[p][1]*J[r][1]+(*CTTI)[3][3]*J[p][1]*J[r][2]+
                                 (*CTTI)[5][2]*J[p][2]*J[r][0]+(*CTTI)[1][2]*J[p][2]*J[r][1]+(*CTTI)[3][2]*J[p][2]*J[r][2]);

      float const Ryz_ij = val*((*CTTI)[4][5]*J[p][0]*J[r][0]+(*CTTI)[3][5]*J[p][0]*J[r][1]+(*CTTI)[2][5]*J[p][0]*J[r][2]+
                                 (*CTTI)[4][1]*J[p][1]*J[r][0]+(*CTTI)[3][1]*J[p][1]*J[r][1]+(*CTTI)[2][1]*J[p][1]*J[r][2]+
                                 (*CTTI)[5][3]*J[p][2]*J[r][0]+(*CTTI)[3][3]*J[p][2]*J[r][1]+(*CTTI)[2][3]*J[p][2]*J[r][2]);

        float const localIncrementx = Rxx_ij * uxnLocal[j] + Rxy_ij * uynLocal[j] + Rxz_ij * uznLocal[j];
        float const localIncrementy = Ryx_ij * uxnLocal[j] + Ryy_ij * uynLocal[j] + Ryz_ij * uznLocal[j];
        float const localIncrementz = Rzx_ij * uxnLocal[j] + Rzy_ij * uynLocal[j] + Rzz_ij * uznLocal[j];
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
  //Build CVTI Matrix
  float CVTI[6][6] = {0.0f};
  CVTI[0][0] = rho * vp * vp * (1.0f + 2.0f*epsilon);
  CVTI[1][1] = CVTI[0][0];
  CVTI[2][0] = rho * sqrtf( powf(vp*vp - vs*vs, 2)
                             + 2.0f*vp*vp*delta*(vp*vp - vs*vs) )
                - rho * vs * vs;
  CVTI[0][2] = CVTI[2][0];
  CVTI[1][2] = CVTI[0][2];
  CVTI[2][1] = CVTI[1][2];
  CVTI[2][2] = rho * vp * vp;
  CVTI[3][3] = rho * vs * vs;
  CVTI[4][4] = CVTI[3][3];
  CVTI[5][5] = rho * vs * vs * (1.0f + 2.0f*gamma);
  CVTI[1][0] = CVTI[0][0] - 2.0f*CVTI[5][5];
  CVTI[0][1] = CVTI[1][0];
  // --- Voigt mapping ---
  int Voigt[3][3] = {{0,5,4},{5,1,3},{4,3,2}};
  // --- Rotation matrix ---
  float R[3][3];
  float ctheta = cosf(theta);
  float stheta = sinf(theta);
  float cphi   = cosf(phi);
  float sphi   = sinf(phi);
  R[0][0] = ctheta*cphi;  R[0][1] = ctheta*sphi;  R[0][2] = -stheta;
  R[1][0] = -sphi;        R[1][1] = cphi;        R[1][2] = 0.0f;
  R[2][0] = stheta*cphi;  R[2][1] = stheta*sphi; R[2][2] = ctheta;
  // ---Build CTTI ---
  for(int l=0; l<3; l++)
      for(int m=0; m<3; m++)
          for(int j=0; j<3; j++)
              for(int i=0; i<3; i++)
              {
                  float coeff = 0.0f;
                  for(int a=0; a<3; a++)
                      for(int b=0; b<3; b++)
                          for(int c=0; c<3; c++)
                              for(int d=0; d<3; d++)
                                  coeff += R[d][i]*R[c][j]*R[b][m]*R[a][l]*CVTI[Voigt[d][c]][Voigt[b][a]];
                  CTTI[Voigt[i][j]][Voigt[m][l]] = coeff;
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
