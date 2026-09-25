#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "Integrals.h"
#include "data_type.h"
#include "sem_solver_acoustoelastic.h"
#include "sem_solver_data.h"

namespace solver {
namespace fe {

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& sponge_size, const bool surface_sponge,
    const float taper_delta) {
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh_ = *typed;
  } else {
    throw std::runtime_error("SEMsolverAcoustoElastic: incompatible mesh type in computeFEInit");
  }

  // The sub-solvers assemble full-mesh mass and damping matrices here; they are
  // zeroed and rebuilt per domain by computeGlobalMassMatrix / computeDampingMatrix below.
  m_acoustic_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  m_elastic_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);

  allocateFEarrays();

  TagElements();
  std::cout << "SEMsolverAcoustoElastic: " << num_acoustic_elements_ << " acoustic elements, " << num_elastic_elements_
            << " elastic elements." << std::endl;

  computeGlobalMassMatrix();
  computeDampingMatrix();

  TagNodes();
  std::cout << "SEMsolverAcoustoElastic: " << num_interface_nodes_ << " interface nodes, "
            << utils::enums::to_string(interface_property_convention_) << "." << std::endl;

  ComputeInterfaceCouplingCoefficients();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::allocateFEarrays() {
  int const nElem = m_mesh_.getNumberOfElements();
  int const nNode = m_mesh_.getNumberOfNodes();

  m_element_type_ = allocateVector<vectorInt>(nElem, "acoustoElasticElementType");
  m_interface_node_index_ = allocateVector<vectorInt>(nNode, "interfaceNodeIndex");

  m_coupling_coeff_x_ = allocateVector<vectorReal>(nNode, "couplingCoeffX");
  m_coupling_coeff_y_ = allocateVector<vectorReal>(nNode, "couplingCoeffY");
  m_coupling_coeff_z_ = allocateVector<vectorReal>(nNode, "couplingCoeffZ");

  auto interface_node_index = m_interface_node_index_;
  auto coupling_coeff_x = m_coupling_coeff_x_;
  auto coupling_coeff_y = m_coupling_coeff_y_;
  auto coupling_coeff_z = m_coupling_coeff_z_;

  Kokkos::parallel_for(
      "allocateFEarrays_init", nNode, KOKKOS_LAMBDA(const int i) {
        interface_node_index[i] = -1;
        coupling_coeff_x[i] = 0.0f;
        coupling_coeff_y[i] = 0.0f;
        coupling_coeff_z[i] = 0.0f;
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::initFEarrays() {}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::initSpongeValues() {
  m_acoustic_solver_.initSpongeValues();
  m_elastic_solver_.initSpongeValues();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::resetGlobalVectors(int numNodes) {
  auto acoustic_w0 = m_acoustic_solver_.getForceVector(0);
  auto elastic_w0 = m_elastic_solver_.getForceVector(0);
  auto elastic_w1 = m_elastic_solver_.getForceVector(1);
  auto elastic_w2 = m_elastic_solver_.getForceVector(2);

  Kokkos::parallel_for(
      "resetGlobalVectors_init", numNodes, KOKKOS_LAMBDA(const int i) {
        acoustic_w0[i] = 0.0f;
        elastic_w0[i] = 0.0f;
        elastic_w1[i] = 0.0f;
        elastic_w2[i] = 0.0f;
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::TagElements() {
  int const nElem = m_mesh_.getNumberOfElements();
  int n_acoustic = 0;
  int n_elastic = 0;

  for (int e = 0; e < nElem; ++e) {
    float vs, rho;
    if constexpr (IS_MODEL_ON_NODES) {
      // Read an interior node: a corner node may carry fluid properties when
      // the >= convention is used at the interface.
      int const mid = ORDER / 2;
      int const gIdx = m_mesh_.globalNodeIndex(e, mid, mid, mid);
      vs = m_mesh_.getModelVsOnNodes(gIdx);
      rho = m_mesh_.getModelRhoOnNodes(gIdx);
    } else {
      vs = m_mesh_.getModelVsOnElement(e);
      rho = m_mesh_.getModelRhoOnElement(e);
    }

    float const mu = rho * vs * vs;
    if (mu < kMuTolerance) {
      m_element_type_[e] = kElementTypeAcoustic;
      ++n_acoustic;
    } else {
      m_element_type_[e] = kElementTypeElastic;
      ++n_elastic;
    }
  }

  num_acoustic_elements_ = n_acoustic;
  num_elastic_elements_ = n_elastic;

  acoustic_elem_list_ = allocateVector<vectorInt>(num_acoustic_elements_, "acousticElemList");
  elastic_elem_list_ = allocateVector<vectorInt>(num_elastic_elements_, "elasticElemList");
  int ia = 0;
  int ie = 0;
  for (int e = 0; e < nElem; ++e) {
    if (m_element_type_[e] == kElementTypeAcoustic)
      acoustic_elem_list_[ia++] = e;
    else
      elastic_elem_list_[ie++] = e;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::TagNodes() {
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();
  int const dim = ORDER + 1;

  vectorInt acoustic_count = allocateVector<vectorInt>(nNode, "acousticCount");
  vectorInt elastic_count = allocateVector<vectorInt>(nNode, "elasticCount");

  Kokkos::parallel_for(
      "TagNodes_initCount", nNode, KOKKOS_LAMBDA(const int i) {
        acoustic_count[i] = 0;
        elastic_count[i] = 0;
      });
  FENCE

  auto elem_type = m_element_type_;
  auto mesh_local = m_mesh_;

  // An interface node is a node shared by at least one acoustic and one elastic element.
  Kokkos::parallel_for(
      "TagNodes_mainLoop", nElem, KOKKOS_LAMBDA(const int e) {
        if (e >= nElem) return;

        int const etype = elem_type[e];
        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const gIdx = mesh_local.globalNodeIndex(e, i, j, k);
              if (etype == kElementTypeAcoustic) {
                ATOMICADD(acoustic_count[gIdx], 1);
              } else {
                ATOMICADD(elastic_count[gIdx], 1);
              }
            }
      });
  FENCE

  int n_interface = 0;
  for (int n = 0; n < nNode; ++n) {
    if (acoustic_count[n] > 0 && elastic_count[n] > 0) ++n_interface;
  }
  num_interface_nodes_ = n_interface;
  n_interface_nodes_ = n_interface;
  m_interface_node_indices_ = allocateVector<vectorInt>(n_interface_nodes_, "interfaceNodeIndices");

  int idx = 0;
  for (int n = 0; n < nNode; ++n) {
    if (acoustic_count[n] > 0 && elastic_count[n] > 0) {
      m_interface_node_index_[n] = idx;
      m_interface_node_indices_[idx] = n;
      ++idx;
    } else {
      m_interface_node_index_[n] = -1;
    }
  }
  FENCE

  // Compact node lists per domain. Interface nodes appear in both lists: the
  // acoustic solver updates pressure there, the elastic solver displacement.
  {
    int n_acou = 0, n_elas = 0;
    for (int n = 0; n < nNode; ++n) {
      if (acoustic_count[n] > 0) ++n_acou;
      if (elastic_count[n] > 0) ++n_elas;
    }
    num_acoustic_nodes_ = n_acou;
    num_elastic_nodes_ = n_elas;
    acoustic_node_list_ = allocateVector<vectorInt>(n_acou, "acousticNodeList");
    elastic_node_list_ = allocateVector<vectorInt>(n_elas, "elasticNodeList");
    int ia = 0, ie = 0;
    for (int n = 0; n < nNode; ++n) {
      if (acoustic_count[n] > 0) acoustic_node_list_[ia++] = n;
      if (elastic_count[n] > 0) elastic_node_list_[ie++] = n;
    }
  }

  // u^{n-1} of the solid, one entry per interface node (indexed by the compact interface index).
  m_ux_nm1_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "uxNm1Iface");
  m_uy_nm1_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "uyNm1Iface");
  m_uz_nm1_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "uzNm1Iface");
  for (int i = 0; i < n_interface_nodes_; ++i) {
    m_ux_nm1_iface_[i] = 0.0f;
    m_uy_nm1_iface_[i] = 0.0f;
    m_uz_nm1_iface_[i] = 0.0f;
  }

  // With node-based models an interface node holds a single set of properties,
  // so the solid and fluid sides are stored separately and swapped into the
  // mesh around each domain kernel.
  if constexpr (IS_MODEL_ON_NODES) {
    // One adjacent elastic element per interface node, used to read the solid side.
    m_interface_adj_elastic_elem_ = allocateVector<vectorInt>(n_interface_nodes_, "interfaceAdjElasticElem");
    for (int i = 0; i < n_interface_nodes_; ++i) m_interface_adj_elastic_elem_[i] = -1;

    for (int ei = 0; ei < num_elastic_elements_; ++ei) {
      int const e = elastic_elem_list_[ei];
      for (int ii = 0; ii < dim; ++ii)
        for (int jj = 0; jj < dim; ++jj)
          for (int kk = 0; kk < dim; ++kk) {
            int const gn = m_mesh_.globalNodeIndex(e, ii, jj, kk);
            int const iface_idx = m_interface_node_index_[gn];
            if (iface_idx >= 0 && m_interface_adj_elastic_elem_[iface_idx] < 0)
              m_interface_adj_elastic_elem_[iface_idx] = e;
          }
    }

    m_vp_solid_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "vpSolidIface");
    m_vs_solid_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "vsSolidIface");
    m_rho_solid_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "rhoSolidIface");
    m_vp_fluid_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "vpFluidIface");
    m_rho_fluid_iface_ = allocateVector<vectorReal>(n_interface_nodes_, "rhoFluidIface");

    for (int i = 0; i < n_interface_nodes_; ++i) {
      int const j = m_interface_node_indices_[i];

      // Fluid side: the builder assigns interface nodes with the >= convention,
      // so the values currently in the mesh are the fluid ones.
      m_vp_fluid_iface_[i] = m_mesh_.getModelVpOnNodes(j);
      m_rho_fluid_iface_[i] = m_mesh_.getModelRhoOnNodes(j);

      if (interface_property_convention_ == utils::enums::interfacePropertyConvention::kSharedOnInterfaceNodes) {
        // The builder gave the interface node a single state meant for both
        // sides, so there is nothing to rebuild and the mass fix below is a
        // no-op.
        m_vp_solid_iface_[i] = m_vp_fluid_iface_[i];
        m_vs_solid_iface_[i] = m_mesh_.getModelVsOnNodes(j);
        m_rho_solid_iface_[i] = m_rho_fluid_iface_[i];
        continue;
      }

      // Solid side: read from a non-interface node of the adjacent elastic
      // element to avoid picking up fluid-contaminated corner properties.
      int const e_adj = m_interface_adj_elastic_elem_[i];
      bool found = false;
      if (e_adj >= 0) {
        for (int ii = 0; ii < dim && !found; ++ii)
          for (int jj = 0; jj < dim && !found; ++jj)
            for (int kk = 0; kk < dim && !found; ++kk) {
              int const g = m_mesh_.globalNodeIndex(e_adj, ii, jj, kk);
              if (m_interface_node_index_[g] < 0) {
                m_vp_solid_iface_[i] = m_mesh_.getModelVpOnNodes(g);
                m_vs_solid_iface_[i] = m_mesh_.getModelVsOnNodes(g);
                m_rho_solid_iface_[i] = m_mesh_.getModelRhoOnNodes(g);
                found = true;
              }
            }
      }
      if (!found) {
        // Degenerate case: the adjacent element has no non-interface node.
        m_vp_solid_iface_[i] = m_vp_fluid_iface_[i];
        m_vs_solid_iface_[i] = 0.0f;
        m_rho_solid_iface_[i] = m_rho_fluid_iface_[i];
      }
    }

    // computeGlobalMassMatrix ran before TagNodes and integrated rho_fluid over
    // the elastic elements. The SEM lumped mass is M_e[j] = rho[j] * J * w_j,
    // so rescaling by rho_solid / rho_fluid is exact.
    auto elastic_mass = m_elastic_solver_.getMassMatrixElastic();
    for (int i = 0; i < n_interface_nodes_; ++i) {
      int const j = m_interface_node_indices_[i];
      float const rho_fluid = m_rho_fluid_iface_[i];
      float const rho_solid = m_rho_solid_iface_[i];
      if (rho_fluid > 0.0f && rho_solid != rho_fluid) elastic_mass[j] *= rho_solid / rho_fluid;
    }
    FENCE
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                             IS_MODEL_ON_NODES>::ComputeInterfaceCouplingCoefficients() {
  constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);

  int const nElem = m_mesh_.getNumberOfElements();

  for (int elementNumber = 0; elementNumber < nElem; ++elementNumber) {
    if (m_element_type_[elementNumber] != kElementTypeAcoustic) continue;

    for (int fi = 0; fi < 6; ++fi) {
      int const f = m_mesh_.getGlobalFace(elementNumber, static_cast<model::CubicFace>(fi));

      // A coupling face has all its nodes on the interface; this excludes lateral faces
      // that only touch it at a corner or an edge.
      int iface_count = 0;
      for (int q = 0; q < numNodesPerFace; ++q) {
        if (m_interface_node_index_[m_mesh_.getGlobalNodeFromFace(f, q)] >= 0) ++iface_count;
      }
      if (iface_count < numNodesPerFace) continue;

      float normal[3];
      m_mesh_.faceNormal(elementNumber, static_cast<model::CubicFace>(fi), normal);

      float coords[4][3];
      for (int j = 0; j < 4; ++j) {
        int const gn = m_mesh_.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
        for (int d = 0; d < 3; ++d) coords[j][d] = m_mesh_.nodeCoord(gn, d);
      }

      for (int q = 0; q < numNodesPerFace; ++q) {
        int const gn = m_mesh_.getGlobalNodeFromFace(f, q);
        float const aux = static_cast<float>(INTEGRAL_TYPE::computeDampingTerm(q, coords));
        m_coupling_coeff_x_[gn] += aux * normal[0];
        m_coupling_coeff_y_[gn] += aux * normal[1];
        m_coupling_coeff_z_[gn] += aux * normal[2];
      }
    }
  }
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::SaveInterfaceUnm1(
    const DataType& data) {
  auto ux_prev = data.m_wavefield.m_elastic.getPreviousField(0);
  auto uy_prev = data.m_wavefield.m_elastic.getPreviousField(1);
  auto uz_prev = data.m_wavefield.m_elastic.getPreviousField(2);
  auto iface_list = m_interface_node_indices_;
  auto ux_nm1 = m_ux_nm1_iface_;
  auto uy_nm1 = m_uy_nm1_iface_;
  auto uz_nm1 = m_uz_nm1_iface_;
  int const n_iface = n_interface_nodes_;

  Kokkos::parallel_for(
      "SaveUnm1Interface_Loop", n_iface, KOKKOS_LAMBDA(const int i) {
        int const j = iface_list[i];
        ux_nm1[i] = ux_prev[j];
        uy_nm1[i] = uy_prev[j];
        uz_nm1[i] = uz_prev[j];
      });
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeGlobalMassMatrix() {
  int const nNode = m_mesh_.getNumberOfNodes();

  auto acoustic_mass = m_acoustic_solver_.getMassMatrixAcoustic();
  auto elastic_mass = m_elastic_solver_.getMassMatrixElastic();
  Kokkos::parallel_for(
      "computeGlobalMassMatrix_init", nNode, KOKKOS_LAMBDA(const int i) {
        acoustic_mass[i] = 0.0f;
        elastic_mass[i] = 0.0f;
      });
  FENCE

  m_acoustic_solver_.computeGlobalMassMatrixMasked(m_element_type_, kElementTypeAcoustic);
  FENCE
  m_elastic_solver_.computeGlobalMassMatrixMasked(m_element_type_, kElementTypeElastic);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeDampingMatrix() {
  int const nNode = m_mesh_.getNumberOfNodes();

  auto acoustic_d0 = m_acoustic_solver_.getDampingMatrix(0);
  auto elastic_d0 = m_elastic_solver_.getDampingMatrix(0);
  auto elastic_d1 = m_elastic_solver_.getDampingMatrix(1);
  auto elastic_d2 = m_elastic_solver_.getDampingMatrix(2);
  Kokkos::parallel_for(
      "computeDampingMatrix_init", nNode, KOKKOS_LAMBDA(const int i) {
        acoustic_d0[i] = 0.0f;
        elastic_d0[i] = 0.0f;
        elastic_d1[i] = 0.0f;
        elastic_d2[i] = 0.0f;
      });
  FENCE

  m_acoustic_solver_.computeDampingMatrixMasked(m_element_type_, kElementTypeAcoustic);
  m_elastic_solver_.computeDampingMatrixMasked(m_element_type_, kElementTypeElastic);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeForces(const float& dt,
                                                                                                const int& timeSample,
                                                                                                DataStruct& data) {
  // The interface coupling is not applied here but in updateSolutionForward,
  // because it acts on the updated fields rather than on the right-hand side.
  // This split lets a distributed driver assemble the force vector at
  // partition boundaries in between.
  auto& myData = dynamic_cast<DataType&>(data);

  SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(myData.m_wavefield.m_acoustic,
                                                                   myData.m_rhs.m_rhs_acoustic);
  SEMsolverData<utils::enums::physicType::kElastic> elastic_data(myData.m_wavefield.m_elastic,
                                                                 myData.m_rhs.m_rhs_elastic);

  resetGlobalVectors(m_mesh_.getNumberOfNodes());
  FENCE

  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE
  m_elastic_solver_.applyRHSTerm(timeSample, dt, elastic_data);
  FENCE

  m_acoustic_solver_.computeElementContributionsFromList(acoustic_data, acoustic_elem_list_, num_acoustic_elements_);
  FENCE
  // Swap in the solid properties at interface nodes for the elastic kernel, then restore the fluid ones.
  if constexpr (IS_MODEL_ON_NODES) {
    for (int i = 0; i < n_interface_nodes_; ++i) {
      int const j = m_interface_node_indices_[i];
      m_mesh_.setModelNodeProps(j, m_vp_solid_iface_[i], m_vs_solid_iface_[i], m_rho_solid_iface_[i]);
    }
  }
  m_elastic_solver_.computeElementContributionsFromList(elastic_data, elastic_elem_list_, num_elastic_elements_);
  FENCE
  if constexpr (IS_MODEL_ON_NODES) {
    for (int i = 0; i < n_interface_nodes_; ++i) {
      int const j = m_interface_node_indices_[i];
      m_mesh_.setModelNodeProps(j, m_vp_fluid_iface_[i], 0.0f, m_rho_fluid_iface_[i]);
    }
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::updateSolutionForward(
    const float& dt, DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);

  SEMsolverData<utils::enums::physicType::kElastic> elastic_data(myData.m_wavefield.m_elastic,
                                                                 myData.m_rhs.m_rhs_elastic);
  if (myData.m_wavefield.hasPrevPrev()) {
    throw std::runtime_error(
        "updateSolutionForward called with 3-buffer wavefield. "
        "Use updateSolutionBackward() for adjoint mode.");
  }
  SaveInterfaceUnm1(myData);
  m_elastic_solver_.updateFieldsFromListForward(dt, elastic_data, elastic_node_list_, num_elastic_nodes_);
  FENCE

  SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(myData.m_wavefield.m_acoustic,
                                                                   myData.m_rhs.m_rhs_acoustic);
  m_acoustic_solver_.updateFieldsFromListForward(dt, acoustic_data, acoustic_node_list_, num_acoustic_nodes_);
  FENCE

  ApplyInterfaceCoupling(dt, myData);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::ApplyInterfaceCoupling(
    float dt, const DataType& data) {
  // Correct the solid with p^n, then the fluid with the discrete solid
  // acceleration that correction has just produced: both corrections are then
  // centred on time n. Moving the traction to p^{n+1} instead breaks that
  // symmetry and slowly injects energy, so the order below matters.
  ApplyCouplingAcousticToElastic(dt, data);
  FENCE
  ApplyCouplingElasticToAcoustic(dt, data);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::updateSolutionBackward(
    const float& dt, DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);

  if (!myData.m_wavefield.hasPrevPrev()) {
    throw std::runtime_error(
        "updateSolutionBackward called with 2-buffer wavefield. "
        "Use updateSolutionForward() for forward mode.");
  }
  SEMsolverData<utils::enums::physicType::kElastic> elastic_data(myData.m_wavefield.m_elastic,
                                                                 myData.m_rhs.m_rhs_elastic);
  // The interface coupling is not applied here: in backward mode the Verlet
  // writes u^{n-1} into the prevPrev buffer, whereas ApplyCoupling* read and
  // correct the previous buffer. The adjoint is therefore uncoupled.
  m_elastic_solver_.updateFieldsFromListBackward(dt, elastic_data, elastic_node_list_, num_elastic_nodes_);
  FENCE

  SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(myData.m_wavefield.m_acoustic,
                                                                   myData.m_rhs.m_rhs_acoustic);
  m_acoustic_solver_.updateFieldsFromListBackward(dt, acoustic_data, acoustic_node_list_, num_acoustic_nodes_);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::ApplyCouplingAcousticToElastic(
    float dt, const DataType& data) {
  float const dt2 = dt * dt;
  float const half_dt = 0.5f * dt;
  auto p_curr = data.m_wavefield.m_acoustic.getCurrentField(0);    // p^n
  auto u_prev_x = data.m_wavefield.m_elastic.getPreviousField(0);  // u_x^{n+1}
  auto u_prev_y = data.m_wavefield.m_elastic.getPreviousField(1);  // u_y^{n+1}
  auto u_prev_z = data.m_wavefield.m_elastic.getPreviousField(2);  // u_z^{n+1}
  auto M_e = m_elastic_solver_.getMassMatrixElastic();
  auto C_ex = m_elastic_solver_.getDampingMatrix(0);
  auto C_ey = m_elastic_solver_.getDampingMatrix(1);
  auto C_ez = m_elastic_solver_.getDampingMatrix(2);
  auto taper_e = m_elastic_solver_.getSpongeTaperCoeff();
  auto mesh_local = m_mesh_;
  auto cx = m_coupling_coeff_x_;
  auto cy = m_coupling_coeff_y_;
  auto cz = m_coupling_coeff_z_;
  auto iface_list = m_interface_node_indices_;
  int const n_iface = n_interface_nodes_;

  Kokkos::parallel_for(
      "ApplyCouplingAcousticToElastic_Loop", n_iface, KOKKOS_LAMBDA(const int i) {
        int const j = iface_list[i];
        if (M_e[j] > 0.0f && !mesh_local.isFreeSurface(j)) {
          // Same denominator and taper the Verlet update applied to the physical
          // RHS: without them this correction is O(dt) inconsistent in the sponge.
          float const aux = -dt2 * p_curr[j] * taper_e[j];
          u_prev_x[j] += cx[j] * aux / (M_e[j] + half_dt * C_ex[j]);
          u_prev_y[j] += cy[j] * aux / (M_e[j] + half_dt * C_ey[j]);
          u_prev_z[j] += cz[j] * aux / (M_e[j] + half_dt * C_ez[j]);
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::ApplyCouplingElasticToAcoustic(
    float dt, const DataType& data) {
  float const half_dt = 0.5f * dt;
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
  auto C_f = m_acoustic_solver_.getDampingMatrix(0);
  auto taper_f = m_acoustic_solver_.getSpongeTaperCoeff();
  auto mesh_local = m_mesh_;
  auto cx = m_coupling_coeff_x_;
  auto cy = m_coupling_coeff_y_;
  auto cz = m_coupling_coeff_z_;
  auto iface_list = m_interface_node_indices_;
  int const n_iface = n_interface_nodes_;

  Kokkos::parallel_for(
      "ApplyCouplingElasticToAcoustic_Loop", n_iface, KOKKOS_LAMBDA(const int i) {
        int const j = iface_list[i];
        if (M_f[j] > 0.0f && !mesh_local.isFreeSurface(j)) {
          // Second time difference of the solid displacement; u_nm1_* are indexed by
          // the compact interface index i, the other fields by the global node j.
          float const fd_x = u_np1_x[j] - 2.0f * u_n_x[j] + u_nm1_x[i];
          float const fd_y = u_np1_y[j] - 2.0f * u_n_y[j] + u_nm1_y[i];
          float const fd_z = u_np1_z[j] - 2.0f * u_n_z[j] + u_nm1_z[i];
          // Same denominator and taper the Verlet update applied to the physical RHS.
          p_prev[j] += taper_f[j] * (cx[j] * fd_x + cy[j] * fd_y + cz[j] * fd_z) / (M_f[j] + half_dt * C_f[j]);
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeOneStep(const float& dt,
                                                                                                 const int& timeSample,
                                                                                                 DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  SEMsolverData<utils::enums::physicType::kElastic> elastic_data(myData.m_wavefield.m_elastic,
                                                                 myData.m_rhs.m_rhs_elastic);
  SEMsolverData<utils::enums::physicType::kAcoustic> acoustic_data(myData.m_wavefield.m_acoustic,
                                                                   myData.m_rhs.m_rhs_acoustic);

  // Elastic step.
  m_elastic_solver_.resetGlobalVectors(nNode);
  FENCE

  m_elastic_solver_.applyRHSTerm(timeSample, dt, elastic_data);
  FENCE

  // With node-based models, swap in the solid properties at interface nodes so
  // the elastic kernel uses the correct lambda, mu and rho.
  if constexpr (IS_MODEL_ON_NODES) {
    for (int i = 0; i < n_interface_nodes_; ++i) {
      int const j = m_interface_node_indices_[i];
      m_mesh_.setModelNodeProps(j, m_vp_solid_iface_[i], m_vs_solid_iface_[i], m_rho_solid_iface_[i]);
    }
  }
  m_elastic_solver_.computeElementContributionsFromList(elastic_data, elastic_elem_list_, num_elastic_elements_);
  FENCE
  if constexpr (IS_MODEL_ON_NODES) {
    for (int i = 0; i < n_interface_nodes_; ++i) {
      int const j = m_interface_node_indices_[i];
      m_mesh_.setModelNodeProps(j, m_vp_fluid_iface_[i], 0.0f, m_rho_fluid_iface_[i]);
    }
  }

  // The previous buffer still holds u^{n-1} here; the Verlet update below overwrites it.
  SaveInterfaceUnm1(myData);

  // u^{n+1} is written into the previous buffer.
  m_elastic_solver_.updateFieldsFromListForward(dt, elastic_data, elastic_node_list_, num_elastic_nodes_);
  FENCE

  // Acoustic step.
  m_acoustic_solver_.resetGlobalVectors(nNode);
  FENCE

  m_acoustic_solver_.applyRHSTerm(timeSample, dt, acoustic_data);
  FENCE

  m_acoustic_solver_.computeElementContributionsFromList(acoustic_data, acoustic_elem_list_, num_acoustic_elements_);
  FENCE

  // p^{n+1} is written into the previous buffer.
  m_acoustic_solver_.updateFieldsFromListForward(dt, acoustic_data, acoustic_node_list_, num_acoustic_nodes_);
  FENCE

  // Enforce the fluid/solid interface conditions on the two predictors.
  ApplyInterfaceCoupling(dt, myData);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void SEMsolverAcoustoElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::outputSolutionValues(
    const int& t, int& e, const vectorReal& field, const char* fieldName) {
  m_acoustic_solver_.outputSolutionValues(t, e, field, fieldName);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_SEM_SOLVER_ACOUSTOELASTIC_IMPL_H_
