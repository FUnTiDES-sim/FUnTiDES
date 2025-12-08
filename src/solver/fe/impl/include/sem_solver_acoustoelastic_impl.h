#include <data_type.h>

#include <array>
#include <cstdlib>

#include "fe/Integrals.hpp"
#include "model_discretization_interface.h"
#include "sem_solver_acoustoelastic.h"

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
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
    throw std::runtime_error("Incompatible mesh type in acoustoelastic solver");
  }

  sponge_size_[0] = sponge_size[0];
  sponge_size_[1] = sponge_size[1];
  sponge_size_[2] = sponge_size[2];
  surface_sponge_ = surface_sponge;
  taper_delta_ = taper_delta;

  int numElements = m_mesh.getNumberOfElements();
  int numNodes = m_mesh.getNumberOfNodes();
  
  m_elementType = allocateVector<VECTOR_INT_VIEW>(numElements, "elementType");
  m_isInterfaceNode = allocateVector<VECTOR_BOOL_VIEW>(numNodes, "isInterfaceNode");

  // Tag elements and nodes
  tagElements();
  FENCE
  tagNodes();
  FENCE

  // Give the mask to the sub-solvers
  m_acousticSolver.setElementMask(&m_elementType, ELEMENT_TYPE_ACOUSTIC);
  m_elasticSolver.setElementMask(&m_elementType, ELEMENT_TYPE_ELASTIC);
  
  m_acousticSolver.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  m_elasticSolver.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  
  std::cout << "AcoustoElastic solver initialized:" << std::endl;
  std::cout << "  - Acoustic elements: " << m_numAcousticElements << std::endl;
  std::cout << "  - Elastic elements: " << m_numElasticElements << std::endl;
  std::cout << "  - Interface nodes: " << m_numInterfaceNodes << std::endl;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    tagElements()
{
  int numElements = m_mesh.getNumberOfElements();
  
  m_numAcousticElements = 0;
  m_numElasticElements = 0;
  
  LOOPHEAD(numElements, elem)
  {
    float mu = 0.0f;
    
    if constexpr (IS_MODEL_ON_NODES)
    {
      float mu_sum = 0.0f;
      int dim = m_mesh.getOrder() + 1;
      for (int i = 0; i < dim; ++i)
      {
        for (int j = 0; j < dim; ++j)
        {
          for (int k = 0; k < dim; ++k)
          {
            int globalIdx = m_mesh.globalNodeIndex(elem, i, j, k);
            float vs = m_mesh.getModelVsOnNodes(globalIdx);
            float rho = m_mesh.getModelRhoOnNodes(globalIdx);
            mu_sum += rho * vs * vs;
          }
        }
      }
      mu = mu_sum / (dim * dim * dim);
    }
    else
    {
      float vs = m_mesh.getModelVsOnElement(elem);
      float rho = m_mesh.getModelRhoOnElement(elem);
      mu = rho * vs * vs;
    }
    
    if (mu < mu_tolerance_)
    {
      m_elementType(elem) = ELEMENT_TYPE_ACOUSTIC; 
      ATOMICADD(m_numAcousticElements, 1);
    }
    else
    {
      m_elementType(elem) = ELEMENT_TYPE_ELASTIC;  
      ATOMICADD(m_numElasticElements, 1);
    }
  }
  LOOPEND
}

template<int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    tagNodes()
{
  int numNodes = m_mesh.getNumberOfNodes();
  int numElements = m_mesh.getNumberOfElements();
  
  // Reset compteur
  m_numInterfaceNodes = 0;
  
  // Initialiser tous les noeuds comme non-interface
  LOOPHEAD(numNodes, node)
  {
    m_isInterfaceNode(node) = false;
  }
  LOOPEND
  FENCE
  
  VECTOR_INT_VIEW acousticCount = allocateVector<VECTOR_INT_VIEW>(numNodes, "acousticCount");
  VECTOR_INT_VIEW elasticCount = allocateVector<VECTOR_INT_VIEW>(numNodes, "elasticCount");
  
  LOOPHEAD(numNodes, node)
  {
    acousticCount(node) = 0;
    elasticCount(node) = 0;
  }
  LOOPEND
  FENCE
  
  LOOPHEAD(numElements, elem)
  {
    int elemType = m_elementType(elem);
    int dim = m_mesh.getOrder() + 1;
    
    for (int i = 0; i < dim; ++i)
    {
      for (int j = 0; j < dim; ++j)
      {
        for (int k = 0; k < dim; ++k)
        {
          int globalIdx = m_mesh.globalNodeIndex(elem, i, j, k);
          
          if (elemType == ELEMENT_TYPE_ACOUSTIC)
          {
            ATOMICADD(acousticCount(globalIdx), 1);
          }
          else if (elemType == ELEMENT_TYPE_ELASTIC)
          {
            ATOMICADD(elasticCount(globalIdx), 1);
          }
        }
      }
    }
  }
  LOOPEND
  FENCE
  
  LOOPHEAD(numNodes, node)
  {
    if (acousticCount(node) > 0 && elasticCount(node) > 0)
    {
      m_isInterfaceNode(node) = true;
      ATOMICADD(m_numInterfaceNodes, 1);
    }
  }
  LOOPEND
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeOneStep(const float &dt, const int &timeSample,
                   SolverBase::DataStruct &data)
{
  // Cast vers le type de données acoustoélastique
  auto &myData = dynamic_cast<SEMsolverDataAcoustoElastic &>(data);

  // Extraire les données acoustiques
  SEMsolverDataAcoustic acousticData(
      myData.m_i1, myData.m_i2,
      myData.m_rhsTerm, myData.m_pnGlobal,
      myData.m_rhsElement, myData.m_rhsWeights);

  // Extraire les données élastiques
  SEMsolverDataElastic elasticData(
      myData.m_i1, myData.m_i2,
      myData.m_rhsTermx, myData.m_rhsTermy, myData.m_rhsTermz,
      myData.m_uxnGlobal, myData.m_uynGlobal, myData.m_uznGlobal,
      myData.m_rhsElement, myData.m_rhsWeights);

  // Appeler le solveur acoustique (traite seulement les éléments acoustiques)
  m_acousticSolver.computeOneStep(dt, timeSample, acousticData);
  FENCE

  // Appeler le solveur élastique (traite seulement les éléments élastiques)
  m_elasticSolver.computeOneStep(dt, timeSample, elasticData);
  FENCE

  // Appliquer le couplage à l'interface
  computeInterfaceCoupling(dt, myData.m_i1, myData.m_i2, data);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeInterfaceCoupling(float dt, int i1, int i2, SolverBase::DataStruct &data)
{
  auto &myData = dynamic_cast<SEMsolverDataAcoustoElastic &>(data);
  
  // 1.  Compute the normals at interface nodes
  // 2. p =-sigma.n
  // 3. ensure continuity of normal displacement
  
  int numNodes = m_mesh.getNumberOfNodes();
  
  LOOPHEAD(numNodes, node)
  {
    if (m_isInterfaceNode(node))
    {
    }
  LOOPEND
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::initFEarrays()
{
  m_acousticSolver.initFEarrays();
  m_elasticSolver.initFEarrays();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::allocateFEarrays()
{
  m_acousticSolver.allocateFEarrays();
  m_elasticSolver.allocateFEarrays();
  
  int numNodes = m_mesh.getNumberOfNodes();
  m_interfaceNormals = allocateArray<ARRAY_REAL_VIEW>(numNodes, 3, "interfaceNormals");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::initSpongeValues()
{
  m_acousticSolver.initSpongeValues();
  m_elasticSolver.initSpongeValues();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::resetGlobalVectors(int numNodes)
{
  m_acousticSolver.resetGlobalVectors(numNodes);
  m_elasticSolver.resetGlobalVectors(numNodes);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::computeGlobalMassMatrix()
{
  m_acousticSolver.computeGlobalMassMatrix();
  m_elasticSolver.computeGlobalMassMatrix();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    outputSolutionValues(const int &indexTimeStep, int &i1,
                         int &myElementSource,
                         const ARRAY_REAL_VIEW &field,
                         const char *fieldName)
{
  if (m_elementType(myElementSource) == ELEMENT_TYPE_ACOUSTIC)
  {
    m_acousticSolver.outputSolutionValues(indexTimeStep, i1, myElementSource, field, fieldName);
  }
  else
  {
    m_elasticSolver.outputSolutionValues(indexTimeStep, i1, myElementSource, field, fieldName);
  }
}






