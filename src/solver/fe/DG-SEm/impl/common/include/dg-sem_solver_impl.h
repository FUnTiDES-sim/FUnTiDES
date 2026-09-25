#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_

#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Integrals.h"
#include "data_type.h"
#include "dg-sem_solver.h"
#include "dg_penalty.h"

namespace solver {
namespace fe {

/**
 * @brief Builds the connectivity, initializes both sub-solvers and splits the mesh into DG and SEM parts.
 * @param[in,out] mesh_in Mesh; must be of type MESH_TYPE.
 * @param[in] sponge_size Forwarded unchanged to the sub-solvers.
 * @param[in] surface_sponge Forwarded unchanged to the sub-solvers.
 * @param[in] taper_delta Forwarded unchanged to the sub-solvers.
 * @throws std::runtime_error If @p mesh_in is not a MESH_TYPE.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& sponge_size, const bool surface_sponge,
    const float taper_delta) {
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh_ = *typed;
  } else {
    throw std::runtime_error("DGSEMsolver: incompatible mesh type in computeFEInit");
  }

  m_face_connectivity_.build(m_mesh_);

  m_SEm_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  // Must precede the DG init: BuildDGInteriorFaceList() below emits face ids in this numbering.
  m_DG_solver_.setFaceConnectivity(m_face_connectivity_);
  m_DG_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);

  m_penalty_factor_ = m_DG_solver_.getPenaltyFactor();

  allocateFEarrays();

  TagElements();
  std::cout << "DGSEMsolver: " << num_SEm_elements_ << " SEm elements, " << num_DG_elements_ << " DG elements."
            << std::endl;

  // The full-mesh assembly done by m_SEm_solver_.computeFEInit() also accumulated the
  // contributions of DG-tagged elements, which doubled (about 2x) the mass of the shared interface
  // nodes and produced a spurious partial reflection along the interface. The SEM weak form only
  // owns the SEM elements, so both matrices are rebuilt restricted to them.
  m_SEm_solver_.computeGlobalMassMatrixMasked(m_element_type_, kElementTypeSEM);
  m_SEm_solver_.computeDampingMatrixMasked(m_element_type_, kElementTypeSEM);

  TagNodes();
  std::cout << "DGSEMsolver: " << num_interface_faces_ << " interface faces." << std::endl;
}

/// @brief Allocates the per-element type tag, size numberOfElements.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  int const nElem = m_mesh_.getNumberOfElements();
  m_element_type_ = allocateVector<vectorInt>(nElem, "DGSEMElementType");
}

/**
 * @brief Tags each element as DG or SEM and builds the two element lists.
 *
 * Uses the split given through setElementTags() when its size matches the mesh; otherwise an
 * element is DG when the z coordinate of its central node is below DG_SEM_interface_z_.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::TagElements() {
  int const nElem = m_mesh_.getNumberOfElements();
  int n_dg = 0;
  int n_sem = 0;

  if (m_external_element_type_.size() == static_cast<size_t>(nElem)) {
    // The z-threshold heuristic is skipped: it only cuts the intended plane on a flat mesh.
    for (int e = 0; e < nElem; ++e) {
      m_element_type_[e] = m_external_element_type_[e];
      if (m_element_type_[e] == kElementTypeDG)
        ++n_dg;
      else
        ++n_sem;
    }
  } else {
    int const mid = ORDER / 2;
    for (int e = 0; e < nElem; ++e) {
      int const gIdx = m_mesh_.globalNodeIndex(e, mid, mid, mid);
      float const zCoord = m_mesh_.nodeCoord(gIdx, 2);
      if (zCoord < DG_SEM_interface_z_) {
        m_element_type_[e] = kElementTypeDG;
        ++n_dg;
      } else {
        m_element_type_[e] = kElementTypeSEM;
        ++n_sem;
      }
    }
  }

  num_DG_elements_ = n_dg;
  num_SEm_elements_ = n_sem;

  DG_elem_list_ = allocateVector<vectorInt>(num_DG_elements_, "DGElemList");
  SEm_elem_list_ = allocateVector<vectorInt>(num_SEm_elements_, "SEmElemList");
  int idg = 0;
  int isem = 0;
  for (int e = 0; e < nElem; ++e) {
    if (m_element_type_[e] == kElementTypeDG)
      DG_elem_list_[idg++] = e;
    else
      SEm_elem_list_[isem++] = e;
  }
}

/**
 * @brief Builds the interface face list, the SEM node list and the DG interior face list.
 *
 * A face is on the DG-SEM interface when it is not a boundary face and each of its nodes belongs
 * to at least one DG element and at least one SEM element.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::TagNodes() {
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();
  int const dim = ORDER + 1;

  vectorInt dg_count = allocateVector<vectorInt>(nNode, "dgCount");
  vectorInt sem_count = allocateVector<vectorInt>(nNode, "semCount");

  Kokkos::parallel_for(
      "TagNodes_initCount", nNode, KOKKOS_LAMBDA(const int i) {
        dg_count[i] = 0;
        sem_count[i] = 0;
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
              if (etype == kElementTypeDG) {
                ATOMICADD(dg_count[gIdx], 1);
              } else {
                ATOMICADD(sem_count[gIdx], 1);
              }
            }
      });
  FENCE

  // The counts are filled on the device; the face loops below run on the host.
  auto h_dg_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dg_count);
  auto h_sem_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, sem_count);

  // Face ids are those of m_face_connectivity_, the same index space as the coupling kernels.
  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  int n_interface = 0;
  for (int f = 0; f < num_faces_fc; ++f) {
    if (m_face_connectivity_.isBoundaryFace(f)) continue;
    bool face_on_interface = true;
    for (int j = 0; j < knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_dg_count(gn) == 0 || h_sem_count(gn) == 0) {
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
    for (int j = 0; j < knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_dg_count(gn) == 0 || h_sem_count(gn) == 0) {
        face_on_interface = false;
        break;
      }
    }
    if (face_on_interface) m_interface_face_indices_[idx++] = f;
  }

  // Compact list of the nodes owned by at least one SEM element.
  {
    int n_sem = 0;
    for (int n = 0; n < nNode; ++n)
      if (h_sem_count(n) > 0) ++n_sem;

    num_SEm_nodes_ = n_sem;
    SEm_node_list_ = allocateVector<vectorInt>(n_sem, "SEmNodeList");

    int isem = 0;
    for (int n = 0; n < nNode; ++n)
      if (h_sem_count(n) > 0) SEm_node_list_[isem++] = n;
  }

  BuildDGInteriorFaceList();
  std::cout << "DGSEMsolver: " << m_n_DG_interior_faces_ << " DG interior faces." << std::endl;
}

/**
 * @brief Builds the list of faces handled by the DG interior flux kernel.
 *
 * Selects every face whose owner element is DG and which is not on the DG-SEM interface.
 * Requires m_interface_face_indices_ to be built.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::BuildDGInteriorFaceList() {
  auto h_elem_type = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_element_type_);
  auto h_iface = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_interface_face_indices_);

  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());

  std::vector<bool> is_iface(num_faces_fc, false);
  for (int i = 0; i < num_interface_faces_; ++i) is_iface[h_iface(i)] = true;

  std::vector<int> result;
  result.reserve(num_faces_fc / 2);
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool dg_adj = (h_elem_type(oe) == kElementTypeDG);
    if (dg_adj) result.push_back(f);
  }

  m_n_DG_interior_faces_ = static_cast<int>(result.size());
  m_DG_interior_face_list_ = allocateVector<vectorInt>(m_n_DG_interior_faces_, "DGInteriorFaceList");
  auto h_list = Kokkos::create_mirror_view(m_DG_interior_face_list_);
  for (int i = 0; i < m_n_DG_interior_faces_; ++i) h_list(i) = result[i];
  Kokkos::deep_copy(m_DG_interior_face_list_, h_list);
}

/**
 * @brief Adds the symmetric SIPG interface flux between the SEM pressure and the DG pressure.
 *
 * Both sides read the current pressure. The DG contribution is accumulated (atomically) into the
 * DG stiffness array of the DG sub-solver and the SEM contribution into the SEM force vector.
 * @param[in] data Coupled wavefield; only the current fields are read.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::ApplyCoupling(const DataType& data) {
  auto mesh_local = m_mesh_;
  auto face_connectivity_local = m_face_connectivity_;
  auto const p_DG = data.m_wavefield.m_DGacoustic.getCurrentField(0);
  auto const p_SEM = data.m_wavefield.m_SEMacoustic.getCurrentField(0);

  auto iface_list = m_interface_face_indices_;
  int const n_iface = num_interface_faces_;
  auto element_type = m_element_type_;
  vectorReal work_sem = m_SEm_solver_.getForceVector(0);
  arrayReal stiff_dg = m_DG_solver_.m_stiff_local_;
  auto const face_to_elem_dof = dgSolver::kFaceToElemDof;
  auto const face_to_elem_dof_depth = dgSolver::kFaceToElemDofAtDepth;
  real_t const penalty_local = m_penalty_factor_;

  Kokkos::parallel_for(
      "ApplyCouplingSEMToDG", n_iface, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = iface_list(_loop_idx);

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        bool const owner_is_dg = (element_type(owner_e) == kElementTypeDG);
        int const dg_e = owner_is_dg ? owner_e : neighbor_e;
        int const sem_e = owner_is_dg ? neighbor_e : owner_e;
        int const fid_dg = owner_is_dg ? fid_o : fid_n;
        int const fid_sem = owner_is_dg ? fid_n : fid_o;

        auto dg_to_sem = [&](int i) {
          return owner_is_dg ? face_connectivity_local.getNeighborFaceDof(f, i)
                             : face_connectivity_local.getOwnerFaceDof(f, i);
        };
        auto sem_to_dg = [&](int i) {
          return owner_is_dg ? face_connectivity_local.getOwnerFaceDof(f, i)
                             : face_connectivity_local.getNeighborFaceDof(f, i);
        };

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni = face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float dg_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(dg_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        dg_coords[iv + 2 * jv + 4 * kv]);
        }

        float sem_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(sem_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        sem_coords[iv + 2 * jv + 4 * kv]);
        }

        float normal_dg[3];
        mesh_local.faceNormal(dg_e, static_cast<model::CubicFace>(fid_dg), normal_dg);

        real_t const face_area = computeFaceArea(faceCoords);
        real_t const inv_rho_dg = 1.0f / mesh_local.getModelRhoOnElement(dg_e);
        real_t const gamma_dg = computeSIPGPenaltyFromArea<ORDER>(face_area, dg_coords, penalty_local);
        real_t const inv_rho_sem = 1.0f / mesh_local.getModelRhoOnElement(sem_e);
        real_t const gamma_sem = computeSIPGPenaltyFromArea<ORDER>(face_area, sem_coords, penalty_local);

        // Face-sized accumulator indexed by DG-side face dof. The coupling flux only touches the
        // (ORDER+1)^2 dofs of the shared face, so an element-sized array would cost a much larger
        // per-thread local-memory footprint and an all-element atomic flush of mostly zeros.
        float stiff_dg_local[knumNodesPerFace] = {0};

        float const neg_normal_dg[3] = {-normal_dg[0], -normal_dg[1], -normal_dg[2]};
        real_t const half_dg = 0.5f * inv_rho_dg;
        real_t const half_sem = 0.5f * inv_rho_sem;

        // The SEM side numbers its dofs globally, so an element dof at depth is resolved through the
        // mesh. The element dof layout is i + j*n + k*n^2 (as in faceLocalToElemLocalAtDepth() and
        // globalNodeIndex()); this only undoes that packing.
        auto sem_node_at = [&](int const elem_dof) {
          constexpr int n = ORDER + 1;
          return mesh_local.globalNodeIndex(sem_e, elem_dof % n, (elem_dof / n) % n, elem_dof / (n * n));
        };

        // One face quadrature point q at a time, both sides fused. The contracted callbacks fold
        // sum_k C_ijk * n_k, so each contribution fires once instead of once per physical
        // direction. They always fire with j == q, so everything derived from q is computed here.
        for (int q = 0; q < knumNodesPerFace; ++q) {
          // Face-normal accumulators (SIPG consistency term). These dofs are off the face, so they
          // bypass the face-sized array.
          float norm_dg[ORDER + 1] = {0};
          float norm_sem[ORDER + 1] = {0};

          int const sem_q = dg_to_sem(q);
          int const gn_dg_q = face_connectivity_local.getGlobalNodeFromFace(f, sem_q);
          int const ej_dg = face_to_elem_dof[fid_dg][q];
          int const gn_sem_q = face_connectivity_local.getGlobalNodeFromFace(f, q);
          int const dg_q = sem_to_dg(q);
          int const ej_perm = face_to_elem_dof[fid_dg][dg_q];

          // Half-weighted pressure jump at q seen from each side; the first contribution of every
          // callback is val times it.
          real_t const dp_dg = half_dg * (p_SEM(gn_dg_q) - p_DG(dg_e, ej_dg));
          real_t const dp_sem = half_sem * (p_DG(dg_e, ej_perm) - p_SEM(gn_sem_q));

          // The two other contributions of each side land on fixed slots with equal magnitude and
          // opposite sign, so one register carries both and is flushed once below.
          float acc_dg = 0.0f;
          float acc_sem = 0.0f;

          // DG side, outward normal normal_dg.
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, dg_coords, fid_dg, normal_dg,
              [&](const int i, const int, const real_t val) {
                stiff_dg_local[i] += val * dp_dg;
                acc_dg -= half_dg * val * p_DG(dg_e, face_to_elem_dof[fid_dg][i]);
              },
              [&](const int m, const int, const real_t val) {
                norm_dg[m] += val * dp_dg;
                acc_dg -= half_dg * val * p_DG(dg_e, face_to_elem_dof_depth[fid_dg][q][m]);
              });

          // SEM side, outward normal -normal_dg.
          INTEGRAL_TYPE::computeInterfaceFluxTermAt(
              q, faceCoords, sem_coords, fid_sem, neg_normal_dg,
              [&](const int i, const int, const real_t val) {
                int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, i);
                ATOMICADD(work_sem(gn_i), val * dp_sem);
                acc_sem -= half_sem * val * p_SEM(gn_i);
              },
              [&](const int m, const int, const real_t val) {
                norm_sem[m] += val * dp_sem;
                acc_sem -= half_sem * val * p_SEM(sem_node_at(face_to_elem_dof_depth[fid_sem][q][m]));
              });

          // Negation is exact in IEEE-754, so each register is mirrored onto the opposite side.
          stiff_dg_local[q] += acc_dg;
          ATOMICADD(work_sem(gn_dg_q), -acc_dg);
          ATOMICADD(work_sem(gn_sem_q), acc_sem);
          stiff_dg_local[dg_q] -= acc_sem;

          // Off-face dofs: written directly, not through the face-sized array.
          for (int m = 0; m <= ORDER; ++m) {
            ATOMICADD(stiff_dg(dg_e, face_to_elem_dof_depth[fid_dg][q][m]), norm_dg[m]);
            ATOMICADD(work_sem(sem_node_at(face_to_elem_dof_depth[fid_sem][q][m])), norm_sem[m]);
          }
        }

        // SIPG penalty and atomic write-back. Both sides use the same damping term at face dof i,
        // so it is computed once. The DG row is complete once its penalty is added, so its flush
        // is done in the same loop.
        for (int i = 0; i < knumNodesPerFace; ++i) {
          real_t const damping_i = INTEGRAL_TYPE::computeDampingTerm(i, faceCoords);

          int const ei = face_to_elem_dof[fid_dg][i];
          int const gn_dg_i = face_connectivity_local.getGlobalNodeFromFace(f, dg_to_sem(i));
          stiff_dg_local[i] += gamma_dg * damping_i * (p_DG(dg_e, ei) - p_SEM(gn_dg_i));
          ATOMICADD(stiff_dg(dg_e, ei), stiff_dg_local[i]);

          int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, i);
          int const ei_perm = face_to_elem_dof[fid_dg][sem_to_dg(i)];
          ATOMICADD(work_sem(gn_i), inv_rho_sem * gamma_sem * damping_i * (p_SEM(gn_i) - p_DG(dg_e, ei_perm)));
        }
      });
}

/**
 * @brief Advances the coupled DG and SEM acoustic fields by one time step.
 *
 * Both sub-solvers read the pressure at step n, so the interface coupling has no temporal lag.
 * @param[in] dt Time step.
 * @param[in] timeSample Index of the current time sample.
 * @param[in,out] data Must be a DataType.
 * @throws std::runtime_error If @p data is in distributed mode.
 * @throws std::bad_cast If @p data is not a DataType.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeOneStep(const float& dt,
                                                                                              const int& timeSample,
                                                                                              DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  if (myData.isDistributed) {
    throw std::runtime_error(
        "computeOneStep called in distributed mode. Use computeForces() -> "
        "synchronize() -> updateSolutionForward().");
  }

  DGsolverDataAcoustic DG_data(myData.m_wavefield.m_DGacoustic, myData.m_rhs.m_rhs_DGacoustic);

  SEMsolverData<utils::enums::physicType::kAcoustic> SEm_data(myData.m_wavefield.m_SEMacoustic,
                                                              myData.m_rhs.m_rhs_SEMacoustic);

  // DG: volume terms and DG-DG interior flux. The interface faces are excluded from the face list.
  m_DG_solver_.m_list_mode_ = true;
  m_DG_solver_.m_elem_list_ = DG_elem_list_;
  m_DG_solver_.m_n_elem_list_ = num_DG_elements_;
  m_DG_solver_.m_face_list_ = m_DG_interior_face_list_;
  m_DG_solver_.m_n_face_list_ = m_n_DG_interior_faces_;

  m_DG_solver_.applyRHSTerm(timeSample, dt, DG_data);
  FENCE
  m_DG_solver_.computeVolumeAndBoundary(num_DG_elements_, DG_data.getCurrentField(0));
  FENCE
  m_DG_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_DG_interior_faces_, DG_data.getCurrentField(0));
  FENCE

  // SEM: source and stiffness; the interface contributes nothing until ApplyCoupling().
  m_SEm_solver_.resetGlobalVectors(nNode);
  FENCE
  m_SEm_solver_.applyRHSTerm(timeSample, dt, SEm_data);
  FENCE
  m_SEm_solver_.computeElementContributionsFromList(SEm_data, SEm_elem_list_, num_SEm_elements_);
  FENCE

  ApplyCoupling(myData);
  FENCE

  // Time integration of both sides (Verlet).
  m_DG_solver_.applyVerlet(num_DG_elements_, dt, DG_data.getCurrentField(0), DG_data.getPreviousField(0));
  m_DG_solver_.m_list_mode_ = false;
  FENCE

  m_SEm_solver_.updateFieldsFromListForward(dt, SEm_data, SEm_node_list_, num_SEm_nodes_);
  FENCE
}

/// @brief Writes a SEM field, one value per node, through the SEM sub-solver.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const vectorReal& field, const char* fieldName) {
  m_SEm_solver_.outputSolutionValues(t, e, field, fieldName);
}

/// @brief Writes a DG field, indexed by (element, dof), through the DG sub-solver.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const arrayReal& field, const char* fieldName) {
  m_DG_solver_.outputSolutionValues(t, e, field, fieldName);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_
