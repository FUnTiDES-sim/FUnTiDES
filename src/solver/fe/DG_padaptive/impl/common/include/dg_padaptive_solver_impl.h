#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_IMPL_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Integrals.h"
#include "dg_padaptive_solver.h"
#include "dg_padaptive_solver_data.h"
#include "dg_solver.h"
#include "dg_solver_impl.h"

namespace solver {
namespace fe {

//============================================================================
// computeFEInit
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::computeFEInit(model::ModelApi<float, int>& mesh_in,
                                               const std::array<float, 3>& sponge_size, const bool surface_sponge,
                                               const float taper_delta) {
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh_ = *typed;
  } else {
    throw std::runtime_error("DGPAdaptiveSolver: incompatible mesh type in computeFEInit");
  }

  m_face_connectivity_.build(m_mesh_);

  // Initialise sub-solvers (mass/damping matrices are overridden below).
  m_pMin_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  m_pMax_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);

  // pMin's own face_connectivity_ was built at ORDER_MIN, but the shared mesh is templated at
  // ORDER_MAX: without this, pMin's own interior "Plus"-type faces (kXPlus/kYPlus/kZPlus)
  // silently fail neighbor node-ID matching (face-normal coordinate stuck at ORDER_MIN instead
  // of the mesh's true far edge). No-op for pMax (already ORDER_MAX == mesh order).
  if constexpr (ORDER_MIN != ORDER_MAX) {
    m_pMin_solver_.rebuildFaceConnectivityGeometry(ORDER_MAX);
  }

  m_penalty_factor_ = m_pMax_solver_.getPenaltyFactor();  // both solvers have the same penalty factor

  // Two dedicated execution-space instances (CUDA streams on the CUDA backend), created once so
  // pMin's and pMax's independent per-step kernel chains can be launched concurrently instead of
  // serializing on the single default stream (see computeOneStep). Equal weight for now -- pMin
  // and pMax element counts differ, but so does per-element cost; revisit if profiling shows one
  // stream idling relative to the other.
  {
    auto exec_instances = Kokkos::Experimental::partition_space(Kokkos::DefaultExecutionSpace{}, 1, 1);
    m_pMin_exec_ = exec_instances[0];
    m_pMax_exec_ = exec_instances[1];
  }

  allocateFEarrays();

  ComputeMortarProjection();

  std::cout << "DGPAdaptiveSolver: ORDER_MIN=" << ORDER_MIN << ", ORDER_MAX=" << ORDER_MAX << std::endl;

  TagElements();
  std::cout << "DGPAdaptiveSolver: " << num_pMin_elements_ << " pMin elements, " << num_pMax_elements_
            << " pMax elements." << std::endl;

  TagNodes();
  std::cout << "DGPAdaptiveSolver: " << num_interface_faces_ << " interface faces." << std::endl;
}

//============================================================================
// allocateFEarrays
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::allocateFEarrays() {
  int const nElem = m_mesh_.getNumberOfElements();
  m_element_type_ = allocateVector<vectorInt>(nElem, "pMinpMaxElementType");
  m_mortar_projection =
      allocateArray2D<arrayReal>(pMaxSolver::knumNodesPerFace, pMinSolver::knumNodesPerFace, "mortarProjectionMatrix");
}

//============================================================================
// TagElements
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::TagElements() {
  int const nElem = m_mesh_.getNumberOfElements();
  int n_pMin = 0;
  int n_pMax = 0;

  for (int e = 0; e < nElem; ++e) {
    int const gIdx = m_mesh_.globalNodeIndex(e, 0, 0, 0);
    float const zCoord = m_mesh_.nodeCoord(gIdx, 2);
    if (zCoord < pAdaptive_interface_z_) {
      m_element_type_[e] = kElementTypePMin;
      ++n_pMin;
    } else {
      m_element_type_[e] = kElementTypePMax;
      ++n_pMax;
    }
  }

  num_pMin_elements_ = n_pMin;
  num_pMax_elements_ = n_pMax;

  pMin_elem_list_ = allocateVector<vectorInt>(num_pMin_elements_, "pMinElemList");
  pMax_elem_list_ = allocateVector<vectorInt>(num_pMax_elements_, "pMaxElemList");
  int ipMin = 0;
  int ipMax = 0;
  for (int e = 0; e < nElem; ++e) {
    if (m_element_type_[e] == kElementTypePMin)
      pMin_elem_list_[ipMin++] = e;
    else
      pMax_elem_list_[ipMax++] = e;
  }
}

//============================================================================
// TagNodes
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::TagNodes() {
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();
  int const dim = ORDER_MAX + 1;

  vectorInt pMin_count = allocateVector<vectorInt>(nNode, "pMinCount");
  vectorInt pMax_count = allocateVector<vectorInt>(nNode, "pMaxCount");

  Kokkos::parallel_for(
      "TagNodes_initCount", nNode, KOKKOS_LAMBDA(const int i) {
        pMin_count[i] = 0;
        pMax_count[i] = 0;
      });
  FENCE

  auto elem_type = m_element_type_;
  auto mesh_local = m_mesh_;

  Kokkos::parallel_for(
      "TagNodes_mainLoop", nElem, KOKKOS_LAMBDA(const int e) {
        if (e >= nElem) return;

        int const etype = elem_type[e];
        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const gIdx = mesh_local.globalNodeIndex(e, i, j, k);
              if (etype == kElementTypePMin) {
                ATOMICADD(pMin_count[gIdx], 1);
              } else {
                ATOMICADD(pMax_count[gIdx], 1);
              }
            }
      });
  FENCE

  // Host mirrors: pMin_count/pMax_count are filled by device kernel; need host access for face loop.
  auto h_pMin_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, pMin_count);
  auto h_pMax_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, pMax_count);

  // Iterate over m_face_connectivity_ faces — same index space as coupling kernels.
  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  int n_interface = 0;
  for (int f = 0; f < num_faces_fc; ++f) {
    if (m_face_connectivity_.isBoundaryFace(f)) continue;
    bool face_on_interface = true;
    for (int j = 0; j < pMaxSolver::knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_pMin_count(gn) == 0 || h_pMax_count(gn) == 0) {
        face_on_interface = false;
        break;
      }
    }
    if (face_on_interface) ++n_interface;
  }

  num_interface_faces_ = n_interface;
  m_interface_face_indices_ = allocateVector<vectorInt>(num_interface_faces_, "interfaceFaceIndices");

  int idx = 0;
  for (int f = 0; f < num_faces_fc; ++f) {
    if (m_face_connectivity_.isBoundaryFace(f)) continue;
    bool face_on_interface = true;
    for (int j = 0; j < pMaxSolver::knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_pMin_count(gn) == 0 || h_pMax_count(gn) == 0) {
        face_on_interface = false;
        break;
      }
    }
    if (face_on_interface) m_interface_face_indices_[idx++] = f;
  }

  BuildInteriorFaceLists();
  std::cout << "DGPAdaptiveSolver: " << m_n_pMin_interior_faces_ << " pMin interior faces." << std::endl;
  std::cout << "DGPAdaptiveSolver: " << m_n_pMax_interior_faces_ << " pMax interior faces." << std::endl;
}

//============================================================================
// BuildDGInteriorFaceList
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::BuildInteriorFaceLists() {
  auto h_elem_type = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_element_type_);
  auto h_iface = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_interface_face_indices_);

  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());

  std::vector<bool> is_iface(num_faces_fc, false);
  for (int i = 0; i < num_interface_faces_; ++i) is_iface[h_iface(i)] = true;

  // Collect faces adjacent to at least one pMin element and not on the pMin-pMax interface.
  std::vector<int> result_pMin;
  result_pMin.reserve(num_faces_fc / 2);
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool pMin_adj = (h_elem_type(oe) == kElementTypePMin);
    if (pMin_adj) result_pMin.push_back(f);
  }

  m_n_pMin_interior_faces_ = static_cast<int>(result_pMin.size());
  m_pMin_interior_face_list_ = allocateVector<vectorInt>(m_n_pMin_interior_faces_, "pMinInteriorFaceList");
  auto h_pMin_list = Kokkos::create_mirror_view(m_pMin_interior_face_list_);
  for (int i = 0; i < m_n_pMin_interior_faces_; ++i) h_pMin_list(i) = result_pMin[i];
  Kokkos::deep_copy(m_pMin_interior_face_list_, h_pMin_list);

  // Collect faces adjacent to at least one pMax element and not on the pMin-pMax interface.
  std::vector<int> result_pMax;
  result_pMax.reserve(num_faces_fc / 2);
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool pMax_adj = (h_elem_type(oe) == kElementTypePMax);
    if (pMax_adj) result_pMax.push_back(f);
  }

  m_n_pMax_interior_faces_ = static_cast<int>(result_pMax.size());
  m_pMax_interior_face_list_ = allocateVector<vectorInt>(m_n_pMax_interior_faces_, "pMaxInteriorFaceList");
  auto h_pMax_list = Kokkos::create_mirror_view(m_pMax_interior_face_list_);
  for (int i = 0; i < m_n_pMax_interior_faces_; ++i) h_pMax_list(i) = result_pMax[i];
  Kokkos::deep_copy(m_pMax_interior_face_list_, h_pMax_list);
}

//============================================================================
// ComputeMortarProjection
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::ComputeMortarProjection() {
  for (int k = 0; k < ORDER_MAX + 1; ++k) {
    for (int l = 0; l < ORDER_MAX + 1; ++l) {
      for (int m = 0; m < ORDER_MIN + 1; ++m) {
        for (int n = 0; n < ORDER_MIN + 1; ++n) {
          m_mortar_projection(k + l * (ORDER_MAX + 1), m + n * (ORDER_MIN + 1)) =
              INTEGRAL_TYPE_MIN::BasisType::value(m, INTEGRAL_TYPE_MAX::BasisType::parentSupportCoord(k)) *
              INTEGRAL_TYPE_MIN::BasisType::value(n, INTEGRAL_TYPE_MAX::BasisType::parentSupportCoord(l));
        }
      }
    }
  }
}

//============================================================================
// ApplyCoupling - SIPG flux:  pMin pressure → pMax pressure
//                             pMax pressure → pMin pressure
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::ApplyCoupling(const DataType& data) {
  auto mesh_local = m_mesh_;
  auto face_connectivity_local = m_face_connectivity_;
  auto const pField_pMin = data.m_wavefield.m_pMinAcoustic.getCurrentField(0);
  auto const pField_pMax = data.m_wavefield.m_pMaxAcoustic.getCurrentField(0);

  auto iface_list = m_interface_face_indices_;
  int const n_iface = num_interface_faces_;
  auto element_type = m_element_type_;
  arrayReal stiff_pMin = m_pMin_solver_.m_stiff_local_;
  arrayReal stiff_pMax = m_pMax_solver_.m_stiff_local_;
  auto const min_face_to_elem_dof = pMinSolver::kFaceToElemDof;
  auto const max_face_to_elem_dof = pMaxSolver::kFaceToElemDof;
  real_t const penalty_local = m_penalty_factor_;
  arrayReal mortar_projection_local = m_mortar_projection;

  Kokkos::parallel_for(
      "ApplyCoupling", n_iface, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = iface_list(_loop_idx);

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        bool const owner_is_pMin = (element_type(owner_e) == kElementTypePMin);
        int const pMin_e = (owner_is_pMin) ? owner_e : neighbor_e;
        int const pMax_e = (owner_is_pMin) ? neighbor_e : owner_e;
        int const fid_pMin = (owner_is_pMin) ? fid_o : fid_n;
        int const fid_pMax = (owner_is_pMin) ? fid_n : fid_o;

        auto pMin_to_pMax = [&](int i) {
          return owner_is_pMin ? face_connectivity_local.getNeighborFaceDof(f, i)
                               : face_connectivity_local.getOwnerFaceDof(f, i);
        };
        auto pMax_to_pMin = [&](int i) {
          return owner_is_pMin ? face_connectivity_local.getOwnerFaceDof(f, i)
                               : face_connectivity_local.getNeighborFaceDof(f, i);
        };

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni =
              face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE_MAX::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float pMin_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(pMin_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        pMin_coords[iv + 2 * jv + 4 * kv]);
        }
        float pMax_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(pMax_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        pMax_coords[iv + 2 * jv + 4 * kv]);
        }

        real_t const inv_rho_min = 1.0f / mesh_local.getModelRhoOnElement(pMin_e);
        real_t const inv_rho_max = 1.0f / mesh_local.getModelRhoOnElement(pMax_e);

        float normal_pMin[3];
        mesh_local.faceNormal(pMin_e, static_cast<model::CubicFace>(fid_pMin), normal_pMin);

        real_t const gamma_min = computeSIPGPenalty<ORDER_MIN>(faceCoords, pMin_coords, penalty_local);
        real_t const gamma_max = computeSIPGPenalty<ORDER_MAX>(faceCoords, pMax_coords, penalty_local);
        // Symmetric penalty: same gamma on both sides of the interface, else the SIPG
        // bilinear form loses symmetry across the hp-nonconforming face (spurious reflection).
        real_t const gamma_iface = (gamma_min > gamma_max) ? gamma_min : gamma_max;

        float stiff_min[pMinSolver::knumNodesPerFace] = {0};
        float stiff_max[pMaxSolver::knumNodesPerFace] = {0};
        float stiff_min_mortar[pMaxSolver::knumNodesPerFace] = {0};

        // Map the pressure field of pMin into pMax on the face
        float pField_mortar[pMaxSolver::knumNodesPerFace] = {0};
        for (int i = 0; i < pMaxSolver::knumNodesPerFace; ++i) {
          for (int j = 0; j < pMinSolver::knumNodesPerFace; ++j) {
            pField_mortar[i] += mortar_projection_local(i, j) * pField_pMin(pMin_e, min_face_to_elem_dof[fid_pMin][j]);
          }
        }

        // Flux computation on pMin element with INTEGRAL_TYPE_MAX
        INTEGRAL_TYPE_MAX::computeInterfaceFluxTerm(
            faceCoords, pMin_coords, fid_pMin, [&](const int i, const int j, const int k, const real_t val) {
              int const j_perm = pMin_to_pMax(j);
              int const ej_perm = max_face_to_elem_dof[fid_pMax][j_perm];
              float const nk = normal_pMin[k];
              stiff_min_mortar[i] +=
                  inv_rho_min * nk * (-0.5f * val * pField_mortar[j] + 0.5f * val * pField_pMax(pMax_e, ej_perm));
              stiff_min_mortar[j] += inv_rho_min * nk * (-0.5f * val * pField_mortar[i]);
              stiff_max[j_perm] += inv_rho_min * nk * (0.5f * val * pField_mortar[i]);
            });

        for (int i = 0; i < pMaxSolver::knumNodesPerFace; ++i) {
          int const ei_perm = max_face_to_elem_dof[fid_pMax][pMin_to_pMax(i)];
          stiff_min_mortar[i] += gamma_iface * INTEGRAL_TYPE_MAX::computeDampingTerm(i, faceCoords) *
                                 (pField_mortar[i] - pField_pMax(pMax_e, ei_perm));
        }

        // Flux computation on pMax element with INTEGRAL_TYPE_MAX
        INTEGRAL_TYPE_MAX::computeInterfaceFluxTerm(
            faceCoords, pMax_coords, fid_pMax, [&](const int i, const int j, const int k, const real_t val) {
              int const ei = max_face_to_elem_dof[fid_pMax][i];
              int const ej = max_face_to_elem_dof[fid_pMax][j];
              int const j_perm = pMax_to_pMin(j);
              float const nk = -normal_pMin[k];
              stiff_max[i] +=
                  inv_rho_max * nk * (-0.5f * val * pField_pMax(pMax_e, ej) + 0.5f * val * pField_mortar[j_perm]);
              stiff_max[j] += inv_rho_max * nk * (-0.5f * val * pField_pMax(pMax_e, ei));
              stiff_min_mortar[j_perm] += inv_rho_max * nk * (0.5f * val * pField_pMax(pMax_e, ei));
            });

        for (int i = 0; i < pMaxSolver::knumNodesPerFace; ++i) {
          int const ei = max_face_to_elem_dof[fid_pMax][i];
          int const i_perm = pMax_to_pMin(i);
          stiff_max[i] += gamma_iface * INTEGRAL_TYPE_MAX::computeDampingTerm(i, faceCoords) *
                          (pField_pMax(pMax_e, ei) - pField_mortar[i_perm]);
        }

        // Map the stiffness matrix back into pMin space
        for (int i = 0; i < pMinSolver::knumNodesPerFace; ++i) {
          for (int j = 0; j < pMaxSolver::knumNodesPerFace; ++j) {
            stiff_min[i] += mortar_projection_local(j, i) * stiff_min_mortar[j];
          }
        }

        for (int i = 0; i < pMinSolver::knumNodesPerFace; ++i)
          ATOMICADD(stiff_pMin(pMin_e, min_face_to_elem_dof[fid_pMin][i]), stiff_min[i]);
        for (int i = 0; i < pMaxSolver::knumNodesPerFace; ++i)
          ATOMICADD(stiff_pMax(pMax_e, max_face_to_elem_dof[fid_pMax][i]), stiff_max[i]);
      });
}

//============================================================================
// computeOneStep  (staggered DG-SEM coupling scheme)
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::computeOneStep(const float& dt, const int& timeSample, DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  if (myData.isDistributed) {
    throw std::runtime_error(
        "computeOneStep called in distributed mode. Use computeForces() -> "
        "synchronize() -> updateSolutionForward().");
  }

  // Sub-solver data views are constructed once and reused throughout the step.
  DGsolverDataAcoustic pMin_data(myData.m_wavefield.m_pMinAcoustic, myData.m_rhs.m_rhs_pMinAcoustic);

  DGsolverDataAcoustic pMax_data(myData.m_wavefield.m_pMaxAcoustic, myData.m_rhs.m_rhs_pMaxAcoustic);

  // ====================================================================================
  // pMin DG: volume + pMin-pMin interior flux (interface faces excluded from face list)
  // ====================================================================================

  m_pMin_solver_.m_list_mode_ = true;
  m_pMin_solver_.m_elem_list_ = pMin_elem_list_;
  m_pMin_solver_.m_n_elem_list_ = num_pMin_elements_;
  m_pMin_solver_.m_face_list_ = m_pMin_interior_face_list_;
  m_pMin_solver_.m_n_face_list_ = m_n_pMin_interior_faces_;

  // pMin and pMax read/write completely disjoint per-sub-solver arrays (separate DGsolver
  // instances, separate m_stiff_local_/m_mass_local_/m_damp_local_/m_rhs_elem_) until
  // ApplyCoupling below, so their phases are launched on two dedicated execution-space
  // instances (m_pMin_exec_/m_pMax_exec_, see computeFEInit) to let both chains co-reside on
  // the GPU instead of serializing on the single default stream. No FENCE between phases
  // within a chain: kernels on the same stream are already ordered by construction.
  m_pMin_solver_.applyRHSTerm(timeSample, dt, pMin_data, m_pMin_exec_);
  m_pMin_solver_.computeVolumeAndBoundary(num_pMin_elements_, pMin_data.getCurrentField(0), m_pMin_exec_);
  m_pMin_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_pMin_interior_faces_, pMin_data.getCurrentField(0),
                                                        m_pMin_exec_);

  // ====================================================================================
  // pMax DG: volume + pMax-pMax interior flux (interface faces excluded from face list)
  // ====================================================================================

  m_pMax_solver_.m_list_mode_ = true;
  m_pMax_solver_.m_elem_list_ = pMax_elem_list_;
  m_pMax_solver_.m_n_elem_list_ = num_pMax_elements_;
  m_pMax_solver_.m_face_list_ = m_pMax_interior_face_list_;
  m_pMax_solver_.m_n_face_list_ = m_n_pMax_interior_faces_;

  m_pMax_solver_.applyRHSTerm(timeSample, dt, pMax_data, m_pMax_exec_);
  m_pMax_solver_.computeVolumeAndBoundary(num_pMax_elements_, pMax_data.getCurrentField(0), m_pMax_exec_);
  m_pMax_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_pMax_interior_faces_, pMax_data.getCurrentField(0),
                                                        m_pMax_exec_);

  // =========================================================================
  // Symmetric SIPG interface coupling: both sides read p^n (no temporal lag).
  // ApplyCoupling reads both m_stiff_local_ arrays, so both streams must be
  // drained first (explicit per-instance fence, not the bare FENCE macro,
  // since a bare Kokkos::fence() is not guaranteed to target these named
  // instances specifically).
  // =========================================================================

  m_pMin_exec_.fence();
  m_pMax_exec_.fence();
  ApplyCoupling(myData);
  FENCE

  // =========================================================================
  // Both Verlets -- still independent (disjoint arrays), overlap them too.
  // =========================================================================

  m_pMin_solver_.applyVerlet(num_pMin_elements_, dt, pMin_data.getCurrentField(0), pMin_data.getPreviousField(0),
                             m_pMin_exec_);
  m_pMin_solver_.m_list_mode_ = false;

  m_pMax_solver_.applyVerlet(num_pMax_elements_, dt, pMax_data.getCurrentField(0), pMax_data.getPreviousField(0),
                             m_pMax_exec_);
  m_pMax_solver_.m_list_mode_ = false;

  // Device fully done (both named streams) before returning to the caller, e.g. Python numpy reads.
  m_pMin_exec_.fence();
  m_pMax_exec_.fence();
}

//============================================================================
// outputSolutionValues (pMin or pMax solution output: p[nElem][nDof_pMax])
//============================================================================

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::outputSolutionValues(const int& t, int& e, const arrayReal& field,
                                                      const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << field(e, 0) << endl;
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_IMPL_H_
