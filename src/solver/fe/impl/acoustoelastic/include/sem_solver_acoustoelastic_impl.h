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
  computeDampingMatrix();

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
  // The coupling coefficient at each interface node is the area-weighted
  // integral of the solid→fluid unit normal n̂ over the interface faces
  // adjacent to that node.  n̂ is obtained from faceNormal() on the acoustic
  // (fluid) element: the inward normal of a fluid element at its interface
  // face points into the fluid interior, which is the solid→fluid direction.
  //
  // This works for any mesh geometry; the horizontal bicouche is a special
  // case where the result is (0, 0, 1) everywhere.
  constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);

  int const nElem = m_mesh_.getNumberOfElements();

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

      // Outward normal of the fluid element at this face = solid→fluid direction.
      float normal[3];
      m_mesh_.faceNormal(elementNumber, static_cast<model::CubicFace>(fi),
                         normal);

      // Gather 4 corner node coordinates for face quadrature.
      float coords[4][3];
      for (int j = 0; j < 4; ++j)
      {
        int const gn = m_mesh_.getGlobalNodeFromFace(
            f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
        for (int d = 0; d < 3; ++d) coords[j][d] = m_mesh_.nodeCoord(gn, d);
      }

      // Accumulate area-weighted normal at each interface face node.
      for (int q = 0; q < numNodesPerFace; ++q)
      {
        int const gn = m_mesh_.getGlobalNodeFromFace(f, q);
        float const aux =
            static_cast<float>(INTEGRAL_TYPE::computeDampingTerm(q, coords));
        m_coupling_coeff_x_[gn] += aux * normal[0];
        m_coupling_coeff_y_[gn] += aux * normal[1];
        m_coupling_coeff_z_[gn] += aux * normal[2];
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

  // Reset sub-solver mass matrices before re-accumulation.
  auto acoustic_mass = m_acoustic_solver_.getMassMatrixAcoustic();
  auto elastic_mass = m_elastic_solver_.getMassMatrixElastic();
  LOOPHEAD(nNode, i)
  {
    acoustic_mass[i] = 0.0f;
    elastic_mass[i] = 0.0f;
  }
  LOOPEND
  FENCE

  // Delegate to each sub-solver using the element type mask: each sub-solver
  // processes only its own elements (acoustic physics for fluid elements,
  // elastic physics for solid elements).
  m_acoustic_solver_.computeGlobalMassMatrixMasked(m_element_type_,
                                                   kElementTypeAcoustic);
  FENCE
  m_elastic_solver_.computeGlobalMassMatrixMasked(m_element_type_,
                                                  kElementTypeElastic);
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

  // Delegate to each sub-solver using the element type mask.
  m_acoustic_solver_.computeDampingMatrixMasked(m_element_type_,
                                                kElementTypeAcoustic);
  m_elastic_solver_.computeDampingMatrixMasked(m_element_type_,
                                               kElementTypeElastic);
  FENCE
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

  SEMsolverData<enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
  SEMsolverData<enums::physicType::kElastic> elastic_data(
      myData.m_wavefield.m_elastic, myData.m_rhs.m_rhs_elastic);

  resetGlobalVectors(m_mesh_.getNumberOfNodes());
  FENCE

  // Apply source terms, then stiffness contributions for both domains.
  // Sub-solvers skip nodes with zero mass (the other domain), so processing all
  // elements is correct: the wavefield is zero at "wrong domain" nodes.
  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE
  m_elastic_solver_.applyRHSTerm(timeSample, dt, elastic_data);
  FENCE

  m_acoustic_solver_.computeElementContributions(acoustic_data);
  FENCE
  m_elastic_solver_.computeElementContributions(elastic_data);
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
  auto elastic_mass = m_elastic_solver_.getMassMatrixElastic();
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

  // Sub-solver data views are constructed once and reused throughout the step.
  // Sub-solvers skip nodes whose mass is zero (the other domain), so calling
  // computeElementContributions over all elements is correct: the wavefield
  // stays zero at "wrong domain" nodes (guaranteed by the mass guard in
  // updateFields).
  SEMsolverData<enums::physicType::kElastic> elastic_data(
      myData.m_wavefield.m_elastic, myData.m_rhs.m_rhs_elastic);
  SEMsolverData<enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);

  // =========================================================================
  // ELASTIC STEP
  // =========================================================================

  // 1. Reset elastic work vectors.
  m_elastic_solver_.resetGlobalVectors(nNode);
  FENCE

  // 1b. Apply elastic source term (zero if no solid source).
  m_elastic_solver_.applyRHSTerm(timeSample, dt, elastic_data);
  FENCE

  // 2. Compute elastic stiffness contributions.
  m_elastic_solver_.computeElementContributions(elastic_data);
  FENCE

  // 3. Apply coupling A→E: add acoustic pressure traction at interface.
  ApplyCouplingAcousticToElastic(myData);
  FENCE

  // 4. Update elastic displacement fields (u^{n+1} written into u_prev).
  m_elastic_solver_.updateFields(dt, elastic_data);
  FENCE

  // =========================================================================
  // ACOUSTIC STEP
  // =========================================================================

  // 5. Reset acoustic work vector.
  m_acoustic_solver_.resetGlobalVectors(nNode);
  FENCE

  // 6. Apply acoustic source term.
  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE

  // 7. Compute acoustic stiffness contributions.
  m_acoustic_solver_.computeElementContributions(acoustic_data);
  FENCE

  // 8. Apply coupling E→A: add elastic normal acceleration at interface.
  //    Uses elastic work vectors from step 2 (not yet reset).
  ApplyCouplingElasticToAcoustic(dt, myData);
  FENCE

  // 9. Update acoustic pressure field (p^{n+1} written into p_prev).
  m_acoustic_solver_.updateFields(dt, acoustic_data);
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
