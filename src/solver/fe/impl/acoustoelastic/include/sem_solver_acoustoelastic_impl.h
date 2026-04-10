#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "data_type.h"
#include "Integrals.h"
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
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::computeFEInit(model::ModelApi<float, int>& mesh_in,
                                      const std::array<float, 3>& sponge_size,
                                      const bool surface_sponge,
                                      const float taper_delta)
{
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in))
  {
    m_mesh_ = *typed;
  }
  else
  {
    throw std::runtime_error(
        "SEMsolverAcoustoElastic: incompatible mesh type in computeFEInit");
  }

  // Initialise sub-solvers (mass/damping matrices are overridden below).
  m_acoustic_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge,
                                   taper_delta);
  m_elastic_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge,
                                  taper_delta);

  allocateFEarrays();

  TagElements();
  std::cout << "SEMsolverAcoustoElastic: " << num_acoustic_elements_
            << " acoustic elements, " << num_elastic_elements_
            << " elastic elements." << std::endl;

  computeGlobalMassMatrix();
  computeDampingMatrix();

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
  m_interface_node_index_ =
      allocateVector<VECTOR_INT_VIEW>(nNode, "interfaceNodeIndex");

  m_coupling_coeff_x_ =
      allocateVector<VECTOR_REAL_VIEW>(nNode, "couplingCoeffX");
  m_coupling_coeff_y_ =
      allocateVector<VECTOR_REAL_VIEW>(nNode, "couplingCoeffY");
  m_coupling_coeff_z_ =
      allocateVector<VECTOR_REAL_VIEW>(nNode, "couplingCoeffZ");

  LOOPHEAD(nNode, i)
  {
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
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::resetGlobalVectors(int numNodes)
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

  for (int e = 0; e < nElem; ++e)
  {
    float vs, rho;
    if constexpr (IS_MODEL_ON_NODES)
    {
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

  acoustic_elem_list_ = allocateVector<VECTOR_INT_VIEW>(num_acoustic_elements_,
                                                        "acousticElemList");
  elastic_elem_list_ =
      allocateVector<VECTOR_INT_VIEW>(num_elastic_elements_, "elasticElemList");
  int ia = 0;
  int ie = 0;
  for (int e = 0; e < nElem; ++e)
  {
    if (m_element_type_[e] == kElementTypeAcoustic)
      acoustic_elem_list_[ia++] = e;
    else
      elastic_elem_list_[ie++] = e;
  }
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
  {
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
  }
  MAINLOOPEND
  FENCE

  int n_interface = 0;
  for (int n = 0; n < nNode; ++n)
  {
    if (acoustic_count[n] > 0 && elastic_count[n] > 0) ++n_interface;
  }
  num_interface_nodes_ = n_interface;
  n_interface_nodes_ = n_interface;
  m_interface_node_indices_ = allocateVector<VECTOR_INT_VIEW>(
      n_interface_nodes_, "interfaceNodeIndices");

  int idx = 0;
  for (int n = 0; n < nNode; ++n)
  {
    if (acoustic_count[n] > 0 && elastic_count[n] > 0)
    {
      m_interface_node_index_[n] = idx;
      m_interface_node_indices_[idx] = n;
      ++idx;
    }
    else
    {
      m_interface_node_index_[n] = -1;
    }
  }
  FENCE

  // Allocate compact nm1 arrays (one entry per interface node).
  m_ux_nm1_iface_ =
      allocateVector<VECTOR_REAL_VIEW>(n_interface_nodes_, "uxNm1Iface");
  m_uy_nm1_iface_ =
      allocateVector<VECTOR_REAL_VIEW>(n_interface_nodes_, "uyNm1Iface");
  m_uz_nm1_iface_ =
      allocateVector<VECTOR_REAL_VIEW>(n_interface_nodes_, "uzNm1Iface");
  for (int i = 0; i < n_interface_nodes_; ++i)
  {
    m_ux_nm1_iface_[i] = 0.0f;
    m_uy_nm1_iface_[i] = 0.0f;
    m_uz_nm1_iface_[i] = 0.0f;
  }
}

//============================================================================
// ComputeInterfaceCouplingCoefficients
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::ComputeInterfaceCouplingCoefficients()
{
  constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);

  int const nElem = m_mesh_.getNumberOfElements();

  for (int elementNumber = 0; elementNumber < nElem; ++elementNumber)
  {
    if (m_element_type_[elementNumber] != kElementTypeAcoustic) continue;

    for (int fi = 0; fi < 6; ++fi)
    {
      int const f = m_mesh_.getGlobalFace(elementNumber,
                                          static_cast<model::CubicFace>(fi));

      // All nodes must be interface nodes (excludes lateral corner faces).
      int iface_count = 0;
      for (int q = 0; q < numNodesPerFace; ++q)
      {
        if (m_interface_node_index_[m_mesh_.getGlobalNodeFromFace(f, q)] >= 0)
          ++iface_count;
      }
      if (iface_count < numNodesPerFace) continue;

      float normal[3];
      m_mesh_.faceNormal(elementNumber, static_cast<model::CubicFace>(fi),
                         normal);

      float coords[4][3];
      for (int j = 0; j < 4; ++j)
      {
        int const gn = m_mesh_.getGlobalNodeFromFace(
            f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
        for (int d = 0; d < 3; ++d) coords[j][d] = m_mesh_.nodeCoord(gn, d);
      }

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

  auto acoustic_mass = m_acoustic_solver_.getMassMatrixAcoustic();
  auto elastic_mass = m_elastic_solver_.getMassMatrixElastic();
  LOOPHEAD(nNode, i)
  {
    acoustic_mass[i] = 0.0f;
    elastic_mass[i] = 0.0f;
  }
  LOOPEND
  FENCE

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
                             IS_MODEL_ON_NODES>::computeForces(const float& dt,
                                                               const int&
                                                                   timeSample,
                                                               DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);

    SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
    SEMsolverData<utils::enums::physicType::kElastic> elastic_data(
      myData.m_wavefield.m_elastic, myData.m_rhs.m_rhs_elastic);

  resetGlobalVectors(m_mesh_.getNumberOfNodes());
  FENCE

  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE
  m_elastic_solver_.applyRHSTerm(timeSample, dt, elastic_data);
  FENCE

  m_acoustic_solver_.computeElementContributionsFromList(
      acoustic_data, acoustic_elem_list_, num_acoustic_elements_);
  FENCE
  m_elastic_solver_.computeElementContributionsFromList(
      elastic_data, elastic_elem_list_, num_elastic_elements_);
  FENCE
}

//============================================================================
// updateSolution  (both domains simultaneously — delegates to sub-solvers)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::updateSolution(const float& dt,
                                                                DataStruct&
                                                                    data)
{
  auto& myData = dynamic_cast<DataType&>(data);

    SEMsolverData<utils::enums::physicType::kElastic> elastic_data(
      myData.m_wavefield.m_elastic, myData.m_rhs.m_rhs_elastic);
  m_elastic_solver_.updateFields(dt, elastic_data);
  FENCE

    SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(
      myData.m_wavefield.m_acoustic, myData.m_rhs.m_rhs_acoustic);
  m_acoustic_solver_.updateFields(dt, acoustic_data);
  FENCE
}

//============================================================================
// ApplyCouplingAcousticToElastic
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::ApplyCouplingAcousticToElastic(float dt,
                                                       const DataType& data)
{
  float const dt2 = dt * dt;
  auto p_curr = data.m_wavefield.m_acoustic.getCurrentField(0);    // p^n
  auto u_prev_x = data.m_wavefield.m_elastic.getPreviousField(0);  // u_x^{n+1}
  auto u_prev_y = data.m_wavefield.m_elastic.getPreviousField(1);  // u_y^{n+1}
  auto u_prev_z = data.m_wavefield.m_elastic.getPreviousField(2);  // u_z^{n+1}
  auto M_e = m_elastic_solver_.getMassMatrixElastic();
  auto cx = m_coupling_coeff_x_;
  auto cy = m_coupling_coeff_y_;
  auto cz = m_coupling_coeff_z_;
  auto iface_list = m_interface_node_indices_;
  int const n_iface = n_interface_nodes_;

  LOOPHEAD(n_iface, i)
  {
    int const j = iface_list[i];
    if (M_e[j] > 0.0f)
    {
      float const aux = -p_curr[j] / M_e[j];
      u_prev_x[j] += dt2 * cx[j] * aux;
      u_prev_y[j] += dt2 * cy[j] * aux;
      u_prev_z[j] += dt2 * cz[j] * aux;
    }
  }
  LOOPEND
}

//============================================================================
// ApplyCouplingElasticToAcoustic
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::ApplyCouplingElasticToAcoustic(const DataType& data)
{
  auto p_prev = data.m_wavefield.m_acoustic.getPreviousField(0);
  auto u_np1_x = data.m_wavefield.m_elastic.getPreviousField(0);
  auto u_np1_y = data.m_wavefield.m_elastic.getPreviousField(1);
  auto u_np1_z = data.m_wavefield.m_elastic.getPreviousField(2);
  auto u_n_x = data.m_wavefield.m_elastic.getCurrentField(0);
  auto u_n_y = data.m_wavefield.m_elastic.getCurrentField(1);
  auto u_n_z = data.m_wavefield.m_elastic.getCurrentField(2);
  auto u_nm1_x = m_ux_nm1_iface_;
  auto u_nm1_y = m_uy_nm1_iface_;
  auto u_nm1_z = m_uz_nm1_iface_;
  auto M_f = m_acoustic_solver_.getMassMatrixAcoustic();
  auto cx = m_coupling_coeff_x_;
  auto cy = m_coupling_coeff_y_;
  auto cz = m_coupling_coeff_z_;
  auto iface_list = m_interface_node_indices_;
  int const n_iface = n_interface_nodes_;

  LOOPHEAD(n_iface, i)
  {
    int const j = iface_list[i];
    if (M_f[j] > 0.0f)
    {
      float const fd_x = u_np1_x[j] - 2.0f * u_n_x[j] + u_nm1_x[i];
      float const fd_y = u_np1_y[j] - 2.0f * u_n_y[j] + u_nm1_y[i];
      float const fd_z = u_np1_z[j] - 2.0f * u_n_z[j] + u_nm1_z[i];
      p_prev[j] += (cx[j] * fd_x + cy[j] * fd_y + cz[j] * fd_z) / M_f[j];
    }
  }
  LOOPEND
}

//============================================================================
// computeOneStep  (staggered elasto-acoustic coupling scheme)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::computeOneStep(const float& dt, const int& timeSample,
                                       DataStruct& data)
{
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  // Sub-solver data views are constructed once and reused throughout the step.
    SEMsolverData<utils::enums::physicType::kElastic> elastic_data(
      myData.m_wavefield.m_elastic, myData.m_rhs.m_rhs_elastic);
    SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(
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

  // 2. Compute elastic stiffness (list: elastic elements only).
  m_elastic_solver_.computeElementContributionsFromList(
      elastic_data, elastic_elem_list_, num_elastic_elements_);
  FENCE

  // 2.5. Save u^{n-1} for interface nodes only (compact array, size
  // n_interface_nodes_).  getPreviousField() still holds u^{n-1} at this
  // point; it will be overwritten by the Verlet below.
  {
    auto ux_prev = elastic_data.getPreviousField(0);
    auto uy_prev = elastic_data.getPreviousField(1);
    auto uz_prev = elastic_data.getPreviousField(2);
    auto iface_list = m_interface_node_indices_;
    auto ux_nm1 = m_ux_nm1_iface_;
    auto uy_nm1 = m_uy_nm1_iface_;
    auto uz_nm1 = m_uz_nm1_iface_;
    int const n_iface = n_interface_nodes_;
    LOOPHEAD(n_iface, i)
    {
      int const j = iface_list[i];
      ux_nm1[i] = ux_prev[j];
      uy_nm1[i] = uy_prev[j];
      uz_nm1[i] = uz_prev[j];
    }
    LOOPEND
    FENCE
  }

  // 3. Elastic Verlet: u^{n+1} written into elastic_data.getPreviousField().
  m_elastic_solver_.updateFields(dt, elastic_data);
  FENCE

  // 4. A→E coupling (GEOS post-Verlet): u^{n+1} += dt²·c·(-p^n)/M_e.
  ApplyCouplingAcousticToElastic(dt, myData);
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

  // 7. Compute acoustic stiffness (list: acoustic elements only).
  m_acoustic_solver_.computeElementContributionsFromList(
      acoustic_data, acoustic_elem_list_, num_acoustic_elements_);
  FENCE

  // 8. Acoustic Verlet: p^{n+1} written into acoustic_data.getPreviousField().
  m_acoustic_solver_.updateFields(dt, acoustic_data);
  FENCE

  // 9. E→A coupling post-Verlet.
  ApplyCouplingElasticToAcoustic(myData);
  FENCE
}

//============================================================================
// outputSolutionValues
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::outputSolutionValues(const int& t, int& e,
                                             const VECTOR_REAL_VIEW& field,
                                             const char* fieldName)
{
  m_acoustic_solver_.outputSolutionValues(t, e, field, fieldName);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_
