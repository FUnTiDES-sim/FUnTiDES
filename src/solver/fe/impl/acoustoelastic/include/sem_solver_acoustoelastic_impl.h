#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "data_type.h"
#include "fe/Integrals.hpp"
#include "model_discretization_interface.h"
#include "sem_solver_acoustoelastic.h"
#include "sem_solver_data.h"

namespace solver
{
namespace fe
{

//============================================================================
// computeFEInit
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    computeFEInit(model::ModelApi<float, int>& mesh_in,
                  const std::array<float, 3>& sponge_size,
                  const bool surface_sponge, const float taper_delta)
{
  // --- 1. Store mesh (needed for TagElements, TagNodes, kernels) ---
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in))
  {
    m_mesh_ = *typed;
  }
  else
  {
    throw std::runtime_error(
        "SEMsolverAcoustoElastic: incompatible mesh type in computeFEInit");
  }

  // --- 2. Initialise sub-solvers ---
  // Each sub-solver allocates its arrays, initialises its sponge taper, and
  // computes mass/damping matrices for ALL elements.  Those matrices are
  // intentionally overridden below by the domain-masked computation.
  m_acoustic_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge,
                                   taper_delta);
  m_elastic_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge,
                                  taper_delta);

  // --- 4. Allocate coupled-solver-specific arrays ---
  allocateFEarrays();

  // --- 5. Classify each element as acoustic or elastic ---
  TagElements();
  std::cout << "SEMsolverAcoustoElastic: " << num_acoustic_elements_
            << " acoustic elements, " << num_elastic_elements_
            << " elastic elements." << std::endl;

  // --- 6. Override mass matrices with domain-masked computation ---
  computeGlobalMassMatrix();

  // --- 7. Override damping matrices with domain-masked computation ---
  //computeDampingMatrix();

  // --- 8. Identify interface nodes and compute coupling coefficients ---
  TagNodes();
  std::cout << "SEMsolverAcoustoElastic: " << num_interface_nodes_
            << " interface nodes." << std::endl;

  ComputeInterfaceCouplingCoefficients();
}

//============================================================================
// allocateFEarrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::allocateFEarrays()
{
  int const nElem = m_mesh_.getNumberOfElements();
  int const nNode = m_mesh_.getNumberOfNodes();

  m_element_type_ =
      allocateVector<VECTOR_INT_VIEW>(nElem, "acoustoElasticElementType");
  m_is_interface_node_ =
      allocateVector<VECTOR_INT_VIEW>(nNode, "isInterfaceNode");
  m_interface_node_index_ =
      allocateVector<VECTOR_INT_VIEW>(nNode, "interfaceNodeIndex");

  // Per-global-node coupling coefficients (area-weighted normal).
  // Initialised to zero; accumulated in ComputeInterfaceCouplingCoefficients.
  m_coupling_coeff_x_ =
      allocateVector<VECTOR_REAL_VIEW>(nNode, "couplingCoeffX");
  m_coupling_coeff_y_ =
      allocateVector<VECTOR_REAL_VIEW>(nNode, "couplingCoeffY");
  m_coupling_coeff_z_ =
      allocateVector<VECTOR_REAL_VIEW>(nNode, "couplingCoeffZ");

  // Initialise interface arrays to sentinel / zero values.
  LOOPHEAD(nNode, i)
  {
    m_is_interface_node_[i] = 0;
    m_interface_node_index_[i] = -1;
    m_coupling_coeff_x_[i] = 0.0f;
    m_coupling_coeff_y_[i] = 0.0f;
    m_coupling_coeff_z_[i] = 0.0f;
  }
  LOOPEND
}

//============================================================================
// initFEarrays / initSpongeValues
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::initFEarrays()
{
  // Sub-solvers handle their own sponge via computeFEInit.
  // The coupled solver's sponge taper is initialised in allocateFEarrays.
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::initSpongeValues()
{
  // Sponge taper is owned by each sub-solver; delegate reinitialisation.
  m_acoustic_solver_.initSpongeValues();
  m_elastic_solver_.initSpongeValues();
}

//============================================================================
// resetGlobalVectors
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::resetGlobalVectors(int
                                                                         numNodes)
{
  auto acoustic_w0 = m_acoustic_solver_.getForceVector(0);
  auto elastic_w0 = m_elastic_solver_.getForceVector(0);
  auto elastic_w1 = m_elastic_solver_.getForceVector(1);
  auto elastic_w2 = m_elastic_solver_.getForceVector(2);

  LOOPHEAD(numNodes, i)
  {
    acoustic_w0[i] = 0.0f;
    elastic_w0[i] = 0.0f;
    elastic_w1[i] = 0.0f;
    elastic_w2[i] = 0.0f;
  }
  LOOPEND
}

//============================================================================
// TagElements
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::TagElements()
{
  int const nElem = m_mesh_.getNumberOfElements();
  int n_acoustic = 0;
  int n_elastic = 0;

  // Sequential host-side loop to populate element type tags and counters.
  for (int e = 0; e < nElem; ++e)
  {
    float vs, rho;
    if constexpr (IS_MODEL_ON_NODES)
    {
      // Use the centroid node (0,0,0) as representative for the element.
      int const gIdx = m_mesh_.globalNodeIndex(e, 0, 0, 0);
      vs = m_mesh_.getModelVsOnNodes(gIdx);
      rho = m_mesh_.getModelRhoOnNodes(gIdx);
    }
    else
    {
      vs = m_mesh_.getModelVsOnElement(e);
      rho = m_mesh_.getModelRhoOnElement(e);
    }

    float const mu = rho * vs * vs;
    if (mu < kMuTolerance)
    {
      m_element_type_[e] = kElementTypeAcoustic;
      ++n_acoustic;
    }
    else
    {
      m_element_type_[e] = kElementTypeElastic;
      ++n_elastic;
    }
  }

  num_acoustic_elements_ = n_acoustic;
  num_elastic_elements_ = n_elastic;
}

//============================================================================
// TagNodes
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::TagNodes()
{
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();
  int const dim = ORDER + 1;

  // Temporary per-node counters: how many acoustic / elastic elements touch
  // each node.
  VECTOR_INT_VIEW acoustic_count =
      allocateVector<VECTOR_INT_VIEW>(nNode, "acousticCount");
  VECTOR_INT_VIEW elastic_count =
      allocateVector<VECTOR_INT_VIEW>(nNode, "elasticCount");

  LOOPHEAD(nNode, i)
  {
    acoustic_count[i] = 0;
    elastic_count[i] = 0;
  }
  LOOPEND
  FENCE

  auto elem_type = m_element_type_;
  auto mesh_local = m_mesh_;

  MAINLOOPHEAD(nElem, e)
  if (e >= nElem) return;

  int const etype = elem_type[e];
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
      {
        int const gIdx = mesh_local.globalNodeIndex(e, i, j, k);
        if (etype == kElementTypeAcoustic)
        {
          ATOMICADD(acoustic_count[gIdx], 1);
        }
        else
        {
          ATOMICADD(elastic_count[gIdx], 1);
        }
      }
  MAINLOOPEND
  FENCE

  // First pass: count interface nodes.
  int n_interface = 0;
  for (int n = 0; n < nNode; ++n)
  {
    if (acoustic_count[n] > 0 && elastic_count[n] > 0)
    {
      ++n_interface;
    }
  }
  num_interface_nodes_ = n_interface;

  // Second pass: fill is_interface and index map.
  int idx = 0;
  for (int n = 0; n < nNode; ++n)
  {
    if (acoustic_count[n] > 0 && elastic_count[n] > 0)
    {
      m_is_interface_node_[n] = 1;
      m_interface_node_index_[n] = idx;
      ++idx;
    }
    else
    {
      m_is_interface_node_[n] = 0;
      m_interface_node_index_[n] = -1;
    }
  }
  FENCE
}

//============================================================================
// ComputeInterfaceCouplingCoefficients
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    ComputeInterfaceCouplingCoefficients()
{
  // V1: horizontal bicouche — outward unit normal (solid→fluid) is (0, 0, 1).
  // For general meshes the normal should be derived from faceNormal()
  // (TODO post-V1 — same approach as computeDampingMatrix elastic branch).
  constexpr float kNx = 0.0f, kNy = 0.0f, kNz = 1.0f;
  constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);

  int const nElem = m_mesh_.getNumberOfElements();

  // Host-side sequential loop (mesh face API is not GPU-safe for ModelUnstruct).
  for (int elementNumber = 0; elementNumber < nElem; ++elementNumber)
  {
    if (m_element_type_[elementNumber] != kElementTypeAcoustic) continue;

    for (int fi = 0; fi < 6; ++fi)
    {
      int const f = m_mesh_.getGlobalFace(
          elementNumber, static_cast<model::CubicFace>(fi));

      // Skip faces that do not touch any interface node.
      bool is_iface_face = false;
      for (int q = 0; q < numNodesPerFace && !is_iface_face; ++q)
      {
        if (m_is_interface_node_[m_mesh_.getGlobalNodeFromFace(f, q)])
          is_iface_face = true;
      }
      if (!is_iface_face) continue;

      // Gather 4 corner node coordinates for face quadrature.
      float coords[4][3];
      for (int j = 0; j < 4; ++j)
      {
        int const gn = m_mesh_.getGlobalNodeFromFace(
            f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
        for (int d = 0; d < 3; ++d) coords[j][d] = m_mesh_.nodeCoord(gn, d);
      }

      // Accumulate area-weighted normal at each face node.
      for (int q = 0; q < numNodesPerFace; ++q)
      {
        int const gn = m_mesh_.getGlobalNodeFromFace(f, q);
        float const aux =
            static_cast<float>(INTEGRAL_TYPE::computeDampingTerm(q, coords));
        m_coupling_coeff_x_[gn] += aux * kNx;
        m_coupling_coeff_y_[gn] += aux * kNy;
        m_coupling_coeff_z_[gn] += aux * kNz;
      }
    }
  }
  FENCE
}

//============================================================================
// computeGlobalMassMatrix  (domain-masked override)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::computeGlobalMassMatrix()
{
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();

  // Reset sub-solver mass matrices before re-accumulation.
  auto acoustic_mass = m_acoustic_solver_.getMassMatrix();
  auto elastic_mass = m_elastic_solver_.getMassMatrix();
  LOOPHEAD(nNode, i)
  {
    acoustic_mass[i] = 0.0f;
    elastic_mass[i] = 0.0f;
  }
  LOOPEND
  FENCE

  constexpr int kPts = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  auto mesh_local = m_mesh_;
  auto elem_type = m_element_type_;

  MAINLOOPHEAD(nElem, elementNumber)
  if (elementNumber >= nElem) return;

  float massLocal[kPts] = {0.0f};
  int const dim = mesh_local.getOrder() + 1;

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, mesh_local,
                                                      transformData);

  INTEGRAL_TYPE::computeMassTerm(
      transformData,
      [&](const int j, const real_t val) { massLocal[j] += val; });

  int const etype = elem_type[elementNumber];

  real_t model_factor = 0.0f;
  if constexpr (!IS_MODEL_ON_NODES)
  {
    if (etype == kElementTypeAcoustic)
    {
      float const vp = mesh_local.getModelVpOnElement(elementNumber);
      float const rho = mesh_local.getModelRhoOnElement(elementNumber);
      model_factor = 1.0f / (vp * vp * rho);
    }
    else
    {
      model_factor = mesh_local.getModelRhoOnElement(elementNumber);
    }
  }

  for (int i = 0; i < mesh_local.getNumberOfPointsPerElement(); ++i)
  {
    int const x = i % dim;
    int const z = (i / dim) % dim;
    int const y = i / (dim * dim);
    int const gIdx = mesh_local.globalNodeIndex(elementNumber, x, y, z);

    if constexpr (IS_MODEL_ON_NODES)
    {
      if (etype == kElementTypeAcoustic)
      {
        float const vp = mesh_local.getModelVpOnNodes(gIdx);
        float const rho = mesh_local.getModelRhoOnNodes(gIdx);
        model_factor = 1.0f / (vp * vp * rho);
      }
      else
      {
        model_factor = mesh_local.getModelRhoOnNodes(gIdx);
      }
    }

    float const contrib = massLocal[i] * model_factor;
    if (etype == kElementTypeAcoustic)
    {
      ATOMICADD(acoustic_mass[gIdx], contrib);
    }
    else
    {
      ATOMICADD(elastic_mass[gIdx], contrib);
    }
  }
  MAINLOOPEND
  FENCE
}

//============================================================================
// computeDampingMatrix  (domain-masked override)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::computeDampingMatrix()
{
  int const nNode = m_mesh_.getNumberOfNodes();

  // Reset sub-solver damping matrices before re-accumulation.
  auto acoustic_d0 = m_acoustic_solver_.getDampingMatrix(0);
  auto elastic_d0 = m_elastic_solver_.getDampingMatrix(0);
  auto elastic_d1 = m_elastic_solver_.getDampingMatrix(1);
  auto elastic_d2 = m_elastic_solver_.getDampingMatrix(2);

  LOOPHEAD(nNode, i)
  {
    acoustic_d0[i] = 0.0f;
    elastic_d0[i] = 0.0f;
    elastic_d1[i] = 0.0f;
    elastic_d2[i] = 0.0f;
  }
  LOOPEND
  FENCE

  int const nElem = m_mesh_.getNumberOfElements();
  auto mesh_local = m_mesh_;
  auto elem_type = m_element_type_;

  MAINLOOPHEAD(nElem, elementNumber)
  if (elementNumber >= nElem) return;

  int const etype = elem_type[elementNumber];

  for (int fi = 0; fi < 6; ++fi)
  {
    int const f =
        mesh_local.getGlobalFace(elementNumber, static_cast<model::CubicFace>(fi));
    if (!mesh_local.isBoundaryFace(f)) continue;

    float coords[4][3];
    for (int j = 0; j < 4; ++j)
    {
      int const gn =
          mesh_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
      for (int d = 0; d < 3; ++d) coords[j][d] = mesh_local.nodeCoord(gn, d);
    }

    constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);

    if (etype == kElementTypeAcoustic)
    {
      // Acoustic absorbing boundary: alpha = 1 / (rho * vp)
      real_t model_rho = 0.0f;
      real_t model_vp = 0.0f;
      real_t alpha = 0.0f;

      if constexpr (!IS_MODEL_ON_NODES)
      {
        model_rho = mesh_local.getModelRhoOnElement(elementNumber);
        model_vp = mesh_local.getModelVpOnElement(elementNumber);
        alpha = 1.0f / (model_rho * model_vp);
      }

      for (int q = 0; q < numNodesPerFace; ++q)
      {
        int const gn = mesh_local.getGlobalNodeFromFace(f, q);
        if (mesh_local.isFreeSurface(gn)) continue;

        if constexpr (IS_MODEL_ON_NODES)
        {
          model_rho = mesh_local.getModelRhoOnNodes(gn);
          model_vp = mesh_local.getModelVpOnNodes(gn);
          alpha = 1.0f / (model_rho * model_vp);
        }

        real_t const incr = alpha * INTEGRAL_TYPE::computeDampingTerm(q, coords);
        ATOMICADD(acoustic_d0[gn], incr);
      }
    }
    else  // kElementTypeElastic
    {
      // Elastic absorbing boundary: velocity-stress formulation
      float normal[3];
      mesh_local.faceNormal(elementNumber, static_cast<model::CubicFace>(fi),
                            normal);
      real_t const nx = normal[0], ny = normal[1], nz = normal[2];

      real_t density = 0.0f;
      real_t vp = 0.0f;
      real_t vs = 0.0f;

      if constexpr (!IS_MODEL_ON_NODES)
      {
        density = mesh_local.getModelRhoOnElement(elementNumber);
        vp = mesh_local.getModelVpOnElement(elementNumber);
        vs = mesh_local.getModelVsOnElement(elementNumber);
      }

      for (int q = 0; q < numNodesPerFace; ++q)
      {
        int const gn = mesh_local.getGlobalNodeFromFace(f, q);
        if (mesh_local.isFreeSurface(gn)) continue;

        if constexpr (IS_MODEL_ON_NODES)
        {
          density = mesh_local.getModelRhoOnNodes(gn);
          vp = mesh_local.getModelVpOnNodes(gn);
          vs = mesh_local.getModelVsOnNodes(gn);
        }

        real_t const aux = density * INTEGRAL_TYPE::computeDampingTerm(q, coords);
        real_t const ix =
            aux * (vp * std::fabs(nx) + vs * std::sqrt(ny * ny + nz * nz));
        real_t const iy =
            aux * (vp * std::fabs(ny) + vs * std::sqrt(nx * nx + nz * nz));
        real_t const iz =
            aux * (vp * std::fabs(nz) + vs * std::sqrt(nx * nx + ny * ny));

        ATOMICADD(elastic_d0[gn], ix);
        ATOMICADD(elastic_d1[gn], iy);
        ATOMICADD(elastic_d2[gn], iz);
      }
    }
  }
  MAINLOOPEND
  FENCE
}

//============================================================================
// ComputeAcousticElementContributions
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    ComputeAcousticElementContributions(const DataType& data)
{
  constexpr int kPts = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  int const nElem = m_mesh_.getNumberOfElements();

  auto mesh_local = m_mesh_;
  auto elem_type = m_element_type_;
  auto work_p = m_acoustic_solver_.getForceVector(0);

  // p_curr is field 0 of the acoustic wavefield.
  auto p_curr = data.m_wavefield.m_acoustic.m_pnGlobalCurr;

  MAINLOOPHEAD(nElem, elementNumber)
  if (elementNumber >= nElem) return;
  if (elem_type[elementNumber] != kElementTypeAcoustic) return;

  int const dim = mesh_local.getOrder() + 1;
  float localField[kPts] = {0.0f};
  float localWork[kPts] = {0.0f};

  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
      {
        int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
        localField[i + j * dim + k * dim * dim] = p_curr(gIdx);
      }

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, mesh_local,
                                                      transformData);

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
          int const gIdx = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
          inv_density = 1.0f / mesh_local.getModelRhoOnNodes(gIdx);
        }
      },
      [&](const int i, const int j, const real_t val) {
        localWork[i] += inv_density * val * localField[j];
      });

  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
      {
        int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
        ATOMICADD(work_p[gIdx], localWork[i + j * dim + k * dim * dim]);
      }
  MAINLOOPEND
}

//============================================================================
// ComputeElasticElementContributions  (isotropic)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    ComputeElasticElementContributions(const DataType& data)
{
  constexpr int kPts = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  int const nElem = m_mesh_.getNumberOfElements();

  auto mesh_local = m_mesh_;
  auto elem_type = m_element_type_;
  auto work_ux = m_elastic_solver_.getForceVector(0);
  auto work_uy = m_elastic_solver_.getForceVector(1);
  auto work_uz = m_elastic_solver_.getForceVector(2);

  auto ux_curr = data.m_wavefield.m_elastic.m_uxnGlobalCurr;
  auto uy_curr = data.m_wavefield.m_elastic.m_uynGlobalCurr;
  auto uz_curr = data.m_wavefield.m_elastic.m_uznGlobalCurr;

  MAINLOOPHEAD(nElem, elementNumber)
  if (elementNumber >= nElem) return;
  if (elem_type[elementNumber] != kElementTypeElastic) return;

  int const dim = mesh_local.getOrder() + 1;
  float localFields[3][kPts] = {{0.0f}};
  float localWork[3][kPts] = {{0.0f}};

  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
      {
        int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
        int const lIdx = i + j * dim + k * dim * dim;
        localFields[0][lIdx] = ux_curr(gIdx);
        localFields[1][lIdx] = uy_curr(gIdx);
        localFields[2][lIdx] = uz_curr(gIdx);
      }

  typename INTEGRAL_TYPE::TransformType transformData;
  model_discretization_interface::gatherTransformData(elementNumber, mesh_local,
                                                      transformData);

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
        float vp, vs, rho;
        if constexpr (IS_MODEL_ON_NODES)
        {
          int const gIdx = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
          vp = mesh_local.getModelVpOnNodes(gIdx);
          vs = mesh_local.getModelVsOnNodes(gIdx);
          rho = mesh_local.getModelRhoOnNodes(gIdx);
        }
        else
        {
          vp = mesh_local.getModelVpOnElement(elementNumber);
          vs = mesh_local.getModelVsOnElement(elementNumber);
          rho = mesh_local.getModelRhoOnElement(elementNumber);
        }

        float const mu = rho * vs * vs;
        float const lambda = rho * (vp * vp - 2.0f * vs * vs);
        float const lp2m = lambda + 2.0f * mu;

        for (int p = 0; p < 3; ++p)
        {
          float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
          for (int r = 0; r < 3; ++r)
          {
            float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
            int const idx = p * 3 + r;
            float const v0 = lp2m * Jp0 * Jr0 + mu * (Jp1 * Jr1 + Jp2 * Jr2);
            float const v1 = mu * Jp0 * Jr0 + lp2m * Jp1 * Jr1 + mu * Jp2 * Jr2;
            float const v2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) + lp2m * Jp2 * Jr2;
            float const v3 = lambda * Jp0 * Jr1 + mu * Jp1 * Jr0;
            float const v4 = lambda * Jp0 * Jr2 + mu * Jp2 * Jr0;
            float const v5 = lambda * Jp1 * Jr2 + mu * Jp2 * Jr1;
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
        int const idx = p * 3 + r;
#ifdef __CUDACC__
        float3 const u = make_float3(localFields[0][j], localFields[1][j],
                                     localFields[2][j]);
        float4 const a = CJflat[idx].a;
        float2 const b = CJflat[idx].b;
        localWork[0][i] +=
            fmaf(val * a.x, u.x, fmaf(val * a.w, u.y, val * b.x * u.z));
        localWork[1][i] +=
            fmaf(val * a.w, u.x, fmaf(val * a.y, u.y, val * b.y * u.z));
        localWork[2][i] +=
            fmaf(val * b.x, u.x, fmaf(val * b.y, u.y, val * a.z * u.z));
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
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
      {
        int const gIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
        int const lIdx = i + j * dim + k * dim * dim;
        ATOMICADD(work_ux[gIdx], localWork[0][lIdx]);
        ATOMICADD(work_uy[gIdx], localWork[1][lIdx]);
        ATOMICADD(work_uz[gIdx], localWork[2][lIdx]);
      }
  MAINLOOPEND
}

//============================================================================
// computeForces  (both domains, no coupling — for potential DD use)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    computeForces(const float& dt, const int& timeSample, DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);

  resetGlobalVectors(m_mesh_.getNumberOfNodes());
  FENCE

  // Apply acoustic source term using the acoustic sub-solver's method.
  // Constructs a temporary SEMsolverData<kAcoustic> view of the acoustic data.
  SEMsolverData<enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE

  ComputeAcousticElementContributions(myData);
  FENCE
  ComputeElasticElementContributions(myData);
  FENCE
}

//============================================================================
// updateSolution  (both domains simultaneously — delegates to sub-solvers)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::updateSolution(
    const float& dt, DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);

  // Delegate update to each sub-solver.  Their mass/damping matrices have been
  // corrected by the coupled-solver's computeGlobalMassMatrix /
  // computeDampingMatrix overrides (domain-masked), so updateFields() will
  // automatically skip nodes whose mass is zero (the other domain's nodes).
  SEMsolverData<enums::physicType::kElastic> elastic_data(
      myData.m_wavefield.m_elastic, RhsElastic{});
  m_elastic_solver_.updateFields(dt, elastic_data);
  FENCE

  SEMsolverData<enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
  m_acoustic_solver_.updateFields(dt, acoustic_data);
  FENCE
}

//============================================================================
// ApplyCouplingAcousticToElastic
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    ApplyCouplingAcousticToElastic(const DataType& data)
{
  // Pressure traction on the solid at interface node j:
  //   T = -p * n̂_s  (outward traction, Komatitsch 2000, eq. 3)
  //
  // Sign convention (see updateFields / applyRHSTerm):
  //   work = K·u − F_ext  →  update: u^{n+1} = … − dt²·work/M
  //
  // F_ext_z = T_z = −p·n_z  →  work_z = K·u_z − (−p·n_z) = K·u_z + p·n_z
  //   ∴  work_u[j] += p * n[j]  (ADD, not subtract)

  auto p_curr = data.m_wavefield.m_acoustic.m_pnGlobalCurr;
  auto work_ux = m_elastic_solver_.getForceVector(0);
  auto work_uy = m_elastic_solver_.getForceVector(1);
  auto work_uz = m_elastic_solver_.getForceVector(2);
  auto is_iface = m_is_interface_node_;
  auto cx = m_coupling_coeff_x_;
  auto cy = m_coupling_coeff_y_;
  auto cz = m_coupling_coeff_z_;
  int const nNode = m_mesh_.getNumberOfNodes();

  LOOPHEAD(nNode, j)
  {
    if (is_iface[j])
    {
      float const p = p_curr(j);
      work_ux[j] += p * cx[j];
      work_uy[j] += p * cy[j];
      work_uz[j] += p * cz[j];
    }
  }
  LOOPEND
}

//============================================================================
// ApplyCouplingElasticToAcoustic
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    ApplyCouplingElasticToAcoustic(float /*dt*/, const DataType& /*data*/)
{
  // Kinematic coupling: solid normal acceleration acts as source for acoustic.
  //
  // With corrected A→E sign, work_e = K·u + p·n  (see ApplyCouplingAcousticToElastic).
  // Newton (solid): M_e·ü_n = −work_e·n  →  ü_n = −(work·n)/M_e = a_n.
  //
  // Sign convention (see updateFields / applyRHSTerm):
  //   work = K·p − F_ext  →  update: p^{n+1} = … − dt²·work/M_f
  //
  // F_ext_acoustic = ü_n  →  work_p = K_f·p − ü_n
  //   ∴  work_p[j] -= a_n  (SUBTRACT, not add)

  auto work_ux = m_elastic_solver_.getForceVector(0);
  auto work_uy = m_elastic_solver_.getForceVector(1);
  auto work_uz = m_elastic_solver_.getForceVector(2);
  auto elastic_mass = m_elastic_solver_.getMassMatrix();
  auto work_p = m_acoustic_solver_.getForceVector(0);
  auto is_iface = m_is_interface_node_;
  auto cx = m_coupling_coeff_x_;
  auto cy = m_coupling_coeff_y_;
  auto cz = m_coupling_coeff_z_;
  int const nNode = m_mesh_.getNumberOfNodes();

  LOOPHEAD(nNode, j)
  {
    if (is_iface[j])
    {
      float const M_e = elastic_mass[j];
      if (M_e > 0.0f)
      {
        float const a_n =
            -(work_ux[j] * cx[j] + work_uy[j] * cy[j] + work_uz[j] * cz[j]) /
            M_e;
        work_p[j] -= a_n;
      }
    }
  }
  LOOPEND
}

//============================================================================
// computeOneStep  (staggered elasto-acoustic coupling scheme)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    computeOneStep(const float& dt, const int& timeSample, DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  // =========================================================================
  // ELASTIC STEP
  // =========================================================================

  // 1. Reset elastic work vectors.
  {
    auto w0 = m_elastic_solver_.getForceVector(0);
    auto w1 = m_elastic_solver_.getForceVector(1);
    auto w2 = m_elastic_solver_.getForceVector(2);
    LOOPHEAD(nNode, i)
    {
      w0[i] = 0.0f;
      w1[i] = 0.0f;
      w2[i] = 0.0f;
    }
    LOOPEND
    FENCE
  }

  // 2. Compute elastic stiffness contributions (elastic elements only).
  ComputeElasticElementContributions(myData);
  FENCE

  // 3. Apply coupling A→E: subtract acoustic pressure traction at interface.
  ApplyCouplingAcousticToElastic(myData);
  FENCE

  // 4. Update elastic displacement fields (u^{n+1} written into u_prev).
  {
    SEMsolverData<enums::physicType::kElastic> elastic_data(
        myData.m_wavefield.m_elastic, RhsElastic{});
    m_elastic_solver_.updateFields(dt, elastic_data);
  }
  FENCE

  // =========================================================================
  // ACOUSTIC STEP
  // =========================================================================

  // 5. Reset acoustic work vector.
  {
    auto w0 = m_acoustic_solver_.getForceVector(0);
    LOOPHEAD(nNode, i) { w0[i] = 0.0f; }
    LOOPEND
    FENCE
  }

  // 6. Apply acoustic source term.
  SEMsolverData<enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE

  // 7. Compute acoustic stiffness contributions (acoustic elements only).
  ComputeAcousticElementContributions(myData);
  FENCE

  // 8. Apply coupling E→A: add elastic normal acceleration at interface.
  //    Uses elastic work vectors (still valid from step 2–3, not yet reset).
  ApplyCouplingElasticToAcoustic(dt, myData);
  FENCE

  // 9. Update acoustic pressure field (p^{n+1} written into p_prev).
  {
    SEMsolverData<enums::physicType::kAcoustic> acoustic_data(
        myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
    m_acoustic_solver_.updateFields(dt, acoustic_data);
  }
  FENCE
}

//============================================================================
// outputSolutionValues
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                              IS_MODEL_ON_NODES>::
    outputSolutionValues(const int& t, int& e, const VECTOR_REAL_VIEW& field,
                         const char* fieldName)
{
  m_acoustic_solver_.outputSolutionValues(t, e, field, fieldName);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_
