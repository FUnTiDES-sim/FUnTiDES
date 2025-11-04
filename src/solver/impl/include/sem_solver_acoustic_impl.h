#include <data_type.h>

#include <array>
#include <cstdlib>

#include "fe/Integrals.hpp"
#include "sem_solver_acoustic.h"

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
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
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeOneStep(const float &dt, const int &timeSample,
                   SolverBase::DataStruct &data)
{
  // Cast to the specific DataStruct type
  auto &myData = dynamic_cast<SEMsolverDataAcoustic &>(data);

  int const &i1 = myData.m_i1;
  int const &i2 = myData.m_i2;
  ARRAY_REAL_VIEW const &rhsTerm = myData.m_rhsTerm;
  ARRAY_REAL_VIEW const &pnGlobal = myData.m_pnGlobal;
  VECTOR_INT_VIEW const &rhsElement = myData.m_rhsElement;
  ARRAY_REAL_VIEW const &rhsWeights = myData.m_rhsWeights;

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, i2, rhsTerm, rhsElement, rhsWeights);
  FENCE
  computeElementContributions(i2, pnGlobal);
  FENCE
  updatePressureField(dt, i1, i2, pnGlobal);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                       IS_MODEL_ON_NODES>::resetGlobalVectors(int numNodes)
{
  LOOPHEAD(numNodes, i) { yGlobal[i] = 0; }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    applyRHSTerm(int timeSample, float dt, int i2,
                 const ARRAY_REAL_VIEW &rhsTerm,
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
          float source = rhsTerm(i, timeSample) * rhsWeights(i, localNodeId);
          yGlobal(nodeRHS) -= source;
        }
      }
    }
  }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeElementContributions(int i2, const ARRAY_REAL_VIEW &pnGlobal)
{
  MAINLOOPHEAD(m_mesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= m_mesh.getNumberOfElements()) return;

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

  real_t inv_density = 0.0f;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    inv_density = 1.0f / m_mesh.getModelRhoOnElement(elementNumber);
  }

  INTEGRAL_TYPE::computeStiffnessTerm(
      cornerCoords, [&](const int qa, const int qb, const int qc, const int i,
                        const int j, const real_t val) {
        if constexpr (IS_MODEL_ON_NODES)
        {
          int const gIndex = m_mesh.globalNodeIndex(elementNumber, qa, qb, qc);
          inv_density = 1.0f / m_mesh.getModelRhoOnNodes(gIndex);
        }
        float localIncrement = inv_density * val * pnLocal[j];
        Y[i] += localIncrement;
      });

  for (int i = 0; i < m_mesh.getNumberOfPointsPerElement(); ++i)
  {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = m_mesh.globalNodeIndex(elementNumber, x, y, z);
    ATOMICADD(yGlobal[gIndex], Y[i]);
  }

  MAINLOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    updatePressureField(float dt, int i1, int i2,
                        const ARRAY_REAL_VIEW &pnGlobal)
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    outputSolutionValues(const int &indexTimeStep, int &i1,
                         int &myElementSource,
                         const ARRAY_REAL_VIEW &fieldGlobal,
                         const char *fieldName)  // ← Nouveau paramètre
{
  cout << "TimeStep=" << indexTimeStep << ";  " << fieldName
       << " @ elementSource location " << myElementSource
       << " after computeOneStep = "
       << fieldGlobal(m_mesh.globalNodeIndex(myElementSource, 0, 0, 0), i1)
       << endl;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                       IS_MODEL_ON_NODES>::initFEarrays()
{
  initSpongeValues();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                       IS_MODEL_ON_NODES>::computeGlobalMassMatrix()
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

  real_t inv_model2 = 0.0f;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    inv_model2 = 1.0f / (m_mesh.getModelVpOnElement(elementNumber) *
                         m_mesh.getModelVpOnElement(elementNumber) *
                         m_mesh.getModelRhoOnElement(elementNumber));
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
      inv_model2 = 1.0f / (m_mesh.getModelVpOnNodes(gIndex) *
                           m_mesh.getModelVpOnNodes(gIndex) *
                           m_mesh.getModelRhoOnNodes(gIndex));
    }
    massMatrixLocal[i] *= inv_model2;
    ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
  }

  MAINLOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                       IS_MODEL_ON_NODES>::allocateFEarrays()
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
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
