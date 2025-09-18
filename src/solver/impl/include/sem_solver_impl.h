//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include <data_type.h>

#include <cstdlib>

#include "fe/Integrals.hpp"
#include "sem_solver.h"

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::computeFEInit(
    model::ModelApi<float, int> &mesh_in)
{
  if (auto *typed_mesh = dynamic_cast<MESH_TYPE *>(&mesh_in))
  {
    m_mesh = *typed_mesh;
  }
  else
  {
    throw std::runtime_error("Incompatible mesh type in solver");
  }
  allocateFEarrays();
  initFEarrays();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::computeOneStep(
    const float &dt, const int &timeSample, SolverBase::DataStruct &data)
{
  // Cast to the specific DataStruct type
  auto &myData = dynamic_cast<SEMsolverData &>(data);

  int const &i1 = myData.m_i1;
  int const &i2 = myData.m_i2;
  ARRAY_REAL_VIEW const &rhsTerm = myData.m_rhsTerm;
  ARRAY_REAL_VIEW const &pnGlobal = myData.m_pnGlobal;
  VECTOR_INT_VIEW const &rhsElement = myData.m_rhsElement;
  ARRAY_REAL_VIEW const &rhsWeights = myData.m_rhsWeights;

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, i2, rhsTerm, rhsElement, pnGlobal, rhsWeights);
  FENCE
  computeElementContributions(i2, pnGlobal);
  FENCE
  updatePressureField(dt, i1, i2, pnGlobal);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::resetGlobalVectors(
    int numNodes)
{
  LOOPHEAD(numNodes, i)
  {
    massMatrixGlobal[i] = 0;
    yGlobal[i] = 0;
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::applyRHSTerm(
    int timeSample, float dt, int i2, const ARRAY_REAL_VIEW &rhsTerm,
    const VECTOR_INT_VIEW &rhsElement, const ARRAY_REAL_VIEW &pnGlobal,
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
          float source = rhsTerm(i, timeSample) * rhsWeights(i, localNodeId);
          yGlobal(nodeRHS) -= source;
        }
      }
    }
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::computeElementContributions(
    int i2, const ARRAY_REAL_VIEW &pnGlobal)
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

  float massMatrixLocal[nPointsElement] = {0};
  float pnLocal[nPointsElement] = {0};
  float Y[nPointsElement] = {0};

  int dim = m_mesh.getOrder() + 1;
  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const globalIdx = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    pnLocal[i] = pnGlobal(globalIdx, i2);
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

  INTEGRAL_TYPE::computeMassMatrixAndStiffnessVector(
      elementNumber, m_mesh.getNumberOfPointsPerElement(), cornerCoords,
      m_precomputedIntegralData, massMatrixLocal, pnLocal, Y);

  auto const inv_model2 = 1.0f / (m_mesh.getModelVpOnElement(elementNumber) *
                                  m_mesh.getModelVpOnElement(elementNumber));
  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    massMatrixLocal[i] *= inv_model2;
    ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
    ATOMICADD(yGlobal[gIndex], Y[i]);
  }

  MAINLOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::updatePressureField(
    float dt, int i1, int i2, const ARRAY_REAL_VIEW &pnGlobal)
{
  float const dt2 = dt * dt;
  LOOPHEAD(m_mesh.getNumberOfNodes(), I)
  {
    pnGlobal(I, i1) = 2 * pnGlobal(I, i2) - pnGlobal(I, i1) -
                      dt2 * yGlobal[I] / massMatrixGlobal[I];
    pnGlobal(I, i1) *= spongeTaperCoeff(I);
    pnGlobal(I, i2) *= spongeTaperCoeff(I);
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::outputPnValues(
    const int &indexTimeStep, int &i1, int &myElementSource,
    const ARRAY_REAL_VIEW &pnGlobal)
{
  cout << "TimeStep=" << indexTimeStep
       << ";  pnGlobal @ elementSource location " << myElementSource
       << " after computeOneStep = "
       << pnGlobal(m_mesh.globalNodeIndex(myElementSource, 0, 0, 0), i1)
       << endl;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::initFEarrays()
{
  INTEGRAL_TYPE::init(m_precomputedIntegralData);
  initSpongeValues();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::allocateFEarrays()
{
  int nbQuadraturePoints = (m_mesh.getOrder() + 1) * (m_mesh.getOrder() + 1) *
                           (m_mesh.getOrder() + 1);

  // shared arrays
  massMatrixGlobal = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(),
                                                      "massMatrixGlobal");
  yGlobal =
      allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(), "yGlobal");

  // sponge allocation
  spongeTaperCoeff = allocateVector<VECTOR_REAL_VIEW>(m_mesh.getNumberOfNodes(),
                                                      "spongeTaperCoeff");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE>::initSpongeValues()
{
  // LOOPHEAD(m_mesh.getNumberOfNodes(), n)
  double sigma_max = 0.01;
  for (int n = 0; n < m_mesh.getNumberOfNodes(); n++)
  {
    float x = m_mesh.nodeCoord(n, 0);
    float y = m_mesh.nodeCoord(n, 1);
    float z = m_mesh.nodeCoord(n, 2);
    float distToFrontierX = min(m_mesh.domainSize(0) - x, x);
    float distToFrontierY = min(m_mesh.domainSize(1) - y, y);
    float distToFrontierZ = min(m_mesh.domainSize(2) - z, z);

    // Find closest distance to domain's boundary
    float minDistToFrontier;
    // if (!isSurface) minDistToFrontier = min(distToFrontierX,
    // min(distToFrontierY, distToFrontierZ));
    minDistToFrontier =
        min(distToFrontierX, min(distToFrontierY, distToFrontierZ));
    // else minDistToFrontier = min(distToFrontierY, distToFrontierZ);

    // Compute taper coefficient using the original Gaussian formula
    if (minDistToFrontier < m_spongeSize)
    {
      // d = distance from absorption boundary
      double d = minDistToFrontier;
      // δ = characteristic width of the Gaussian
      double delta = m_spongeSize / 3.0;
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
  // LOOPEND

  FENCE

  // Debuging: Wrinting down all taper coef
  std::ofstream outfile("spongeTaperCoeff.txt");
  for (int i = 0; i < m_mesh.getNumberOfNodes(); i++)
  {
    outfile << spongeTaperCoeff(i) << ' ';
  }

  outfile << endl;
  outfile.close();
}
