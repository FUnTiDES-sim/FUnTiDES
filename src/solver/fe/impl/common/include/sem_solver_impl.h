#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#include <data_type.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>

#include "Integrals.h"
#include "elastic_flux.h"
#include "mesh_type_traits.h"
#include "sem_solver.h"

namespace solver {
namespace fe {

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& sponge_size, const bool surface_sponge,
    const float taper_delta) {
  if (auto* typed_mesh = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh = *typed_mesh;
  } else {
    throw std::runtime_error("Incompatible mesh type in solver");
  }

  // The caches below are derived from the model, which has just been replaced.
  cttiNodesReady_ = false;
  gemmMetricsReady_ = false;
  zDeformedReady_ = false;

  sponge_size_[0] = sponge_size[0];
  sponge_size_[1] = sponge_size[1];
  sponge_size_[2] = sponge_size[2];
  surface_sponge_ = surface_sponge;
  taper_delta_ = taper_delta;

  allocateFEarrays();
  initFEarrays();

  computeGlobalMassMatrix();
  computeDampingMatrix();

  if (attenuationEnabled_ && nSls_ > 0) {
    // Smallest quality factor over the elements, sampled at local node (0,0,0) of each element.
    float minQVal = std::numeric_limits<float>::max();
    for (int e = 0; e < m_mesh.getNumberOfElements(); ++e) {
      if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
        float q = IS_MODEL_ON_NODES ? m_mesh.getModelQpOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                                    : m_mesh.getModelQpOnElement(e);
        minQVal = std::min(minQVal, q);
      } else {
        float qp = IS_MODEL_ON_NODES ? m_mesh.getModelQpOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                                     : m_mesh.getModelQpOnElement(e);
        float qs = IS_MODEL_ON_NODES ? m_mesh.getModelQsOnNodes(m_mesh.globalNodeIndex(e, 0, 0, 0))
                                     : m_mesh.getModelQsOnElement(e);
        minQVal = std::min(minQVal, std::min(qp, qs));
      }
    }

    // A negative coefficient means "not set by the user": it is replaced by a default derived from minQVal.
    // @todo VERIFY: origin and meaning of the default 2*Qmin/(max(1.0001, Qmin) - 1).
    for (int l = 0; l < nSls_; ++l) {
      if (slsAnelasticityCoefficients_.extent(0) > 0 && slsAnelasticityCoefficients_[l] < 0.0f) {
        slsAnelasticityCoefficients_[l] = 2.0f * minQVal / (std::max(1.0001f, minQVal) - 1.0f);
      }
    }
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeForces(const float& dt,
                                                                                           const int& timeSample,
                                                                                           Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);

  resetGlobalVectors(m_mesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, myData);
  FENCE
  computeElementContributions(myData);
  FENCE
  if (attenuationEnabled_ && nSls_ > 0) {
    computeAttenuationContributions(myData);
    FENCE
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolutionForward(
    const float& dt, Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  if (myData.getPrevPrevField(0).extent(0) > 0) {
    throw std::runtime_error(
        "updateSolutionForward called with 3-buffer wavefield. "
        "Use updateSolutionBackward() for adjoint mode.");
  }
  updateFieldsForward(dt, myData);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateSolutionBackward(
    const float& dt, Solver::DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  if (myData.getPrevPrevField(0).extent(0) == 0) {
    throw std::runtime_error(
        "updateSolutionBackward requires 3-buffer wavefield. "
        "Construct wavefield with prevprev buffer (3 args for acoustic, 9 for elastic).");
  }
  updateFieldsBackward(dt, myData);
  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::resetGlobalVectors(int numNodes) {
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  std::array<std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>, kNumFields>
      local_attenuationWorkVectorsGlobal;

  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    if (has_attenuation) {
      local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];
    }
  }

  Kokkos::parallel_for(
      "Solver Reset GVector", numNodes, KOKKOS_LAMBDA(const int i) {
        for (int f = 0; f < kNumFields; ++f) {
          local_workVectorsGlobal[f][i] = 0;
          if (has_attenuation) {
            local_attenuationWorkVectorsGlobal[f][i] = 0;
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::resetGlobalVectorsFromList(
    const vectorInt& node_list, int n_nodes) {
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  std::array<std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>, kNumFields>
      local_attenuationWorkVectorsGlobal;

  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    if (has_attenuation) {
      local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];
    }
  }
  auto list_local = node_list;

  Kokkos::parallel_for(
      "Solver Reset GVector From List", n_nodes, KOKKOS_LAMBDA(const int n) {
        int const i = list_local[n];
        for (int f = 0; f < kNumFields; ++f) {
          local_workVectorsGlobal[f][i] = 0;
          if (has_attenuation) {
            local_attenuationWorkVectorsGlobal[f][i] = 0;
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyRHSTerm(int timeSample, float dt,
                                                                                          const DataType& data) {
  constexpr int kDim = ORDER + 1;
  constexpr int kPointsPerElem = kDim * kDim * kDim;
  auto const rhs_element = data.getRhsElement();
  int const nb_rhs_element = rhs_element.extent(0);
  auto mesh_local = m_mesh;

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }
  std::array<arrayReal, kNumRhs> rhs_term;
  std::array<arrayReal, kNumRhs> rhs_weights;
  for (int f = 0; f < kNumRhs; ++f) {
    rhs_term[f] = data.getRhsTerm(f);
    rhs_weights[f] = data.getRhsWeights(f);
  }

  // One thread per (source element, node): source elements may share nodes, hence the atomics.
  Kokkos::parallel_for(
      "Solver Apply RHSTerm", detail::lightWeightRange(nb_rhs_element * kPointsPerElem), KOKKOS_LAMBDA(const int t) {
        int const i = t / kPointsPerElem;
        int const localNodeId = t - i * kPointsPerElem;
        int const x = localNodeId % kDim;
        int const y = (localNodeId / kDim) % kDim;
        int const z = localNodeId / (kDim * kDim);
        int const nodeRHS = mesh_local.globalNodeIndex(rhs_element[i], x, y, z);

        for (int f = 0; f < kNumRhs; ++f) {
          float const source = rhs_term[f](i, timeSample) * rhs_weights[f](i, localNodeId);
          ATOMICADD(local_workVectorsGlobal[f](nodeRHS), -source);
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions(
    const DataType& data) {
  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    if (anisotropyType_ == model::AnisotropyType::kIso) {
      computeElementContributions_Iso(data);
    } else if (anisotropyType_ == model::AnisotropyType::kVTI) {
      computeElementContributions_Vti(data);
    } else {
      computeElementContributions_Tti(data);
    }
  } else {
    computeElementContributions_Acoustic(data);
  }
}

namespace detail {
/// @brief True when T declares a nested type `TeamGemm`, i.e. the integral back-end provides the team GEMM kernels.
template <typename, typename = void>
struct has_team_gemm : std::false_type {};
template <typename T>
struct has_team_gemm<T, std::void_t<typename T::TeamGemm>> : std::true_type {};

/// @brief True when T declares a nested type `ZDeformedSumFact`, i.e. the integral back-end provides
/// the z-deformed sum-factorization kernels.
template <typename, typename = void>
struct has_z_deformed_sumfact : std::false_type {};
template <typename T>
struct has_z_deformed_sumfact<T, std::void_t<typename T::ZDeformedSumFact>> : std::true_type {};
}  // namespace detail

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_Gemm(
    const DataType& data) {
  if constexpr (detail::has_team_gemm<INTEGRAL_TYPE>::value) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    constexpr int kStride = INTEGRAL_TYPE::numNodes * 6;

    auto mesh_local = m_mesh;
    int const dim = mesh_local.getOrder() + 1;
    int const pointsPerElem = dim * dim * dim;
    int const nElemsFull = mesh_local.getNumberOfElements();

    // W = w * alpha * B is precomputed once for every element, kStride values per element.
    // It runs on the first call (the benchmark warmup), so it stays out of the timed loop.
    if (!gemmMetricsReady_) {
      gemmMetrics_ = allocateVector<vectorReal>(static_cast<size_t>(nElemsFull) * kStride, "gemmMetrics");
      auto W_global = gemmMetrics_;
      auto mesh_pc = m_mesh;
      Kokkos::parallel_for(
          "Gemm Precompute Metrics", Kokkos::RangePolicy<ExecSpace>(0, nElemsFull), KOKKOS_LAMBDA(const int e) {
            float cornerCoords[8][3];
            auto const eIdx = mesh_pc.elementIndex(e);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_pc.vertexCoords(mesh_pc.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);

            real_t inv_density = 0.0f;
            if constexpr (!IS_MODEL_ON_NODES) {
              inv_density = 1.0f / mesh_pc.getModelRhoOnElement(e);
            }

            INTEGRAL_TYPE::computeElementMetrics(
                cornerCoords,
                [&](const int qa, const int qb, const int qc) -> real_t {
                  if constexpr (IS_MODEL_ON_NODES) {
                    int const g = mesh_pc.globalNodeIndex(e, qa, qb, qc);
                    return 1.0f / mesh_pc.getModelRhoOnNodes(g);
                  } else {
                    return inv_density;
                  }
                },
                W_global.data() + static_cast<size_t>(e) * kStride);
          });
      Kokkos::fence();
      gemmMetricsReady_ = true;
    }

    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : nElemsFull;

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) local_workVectorsGlobal[f] = workVectorsGlobal_[f];

    auto W_global = gemmMetrics_;

    // 1D derivative operator, built on host, captured by value (no device View).
    Kokkos::Array<real_t, INTEGRAL_TYPE::num1dNodes * INTEGRAL_TYPE::num1dNodes> D_arr;
    INTEGRAL_TYPE::fillDerivativeMatrix(D_arr.data());

    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView2D = Kokkos::View<float**, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    TeamPolicyType policy(n_iter, Kokkos::AUTO);
    // Scratch holds the gather/scatter buffers and the 12 [N][N^2] GEMM tensors.
    // W lives in global memory, so it needs no scratch: less shared memory than
    // the streaming variant, which also relaxes the occupancy limit slightly.
    size_t bytes_fields = ScratchView2D::shmem_size(kNumFields, pointsPerElem) * 2;
    size_t bytes_gemm = ScratchView2D::shmem_size(INTEGRAL_TYPE::num1dNodes, INTEGRAL_TYPE::numNodesPerFace) * 12;
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_gemm));

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic Gemm", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          ScratchView2D localFields(team.team_scratch(0), kNumFields, pointsPerElem);
          ScratchView2D localWork(team.team_scratch(0), kNumFields, pointsPerElem);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pointsPerElem), [&](const int localIdx) {
            int i = localIdx % dim;
            int j = (localIdx / dim) % dim;
            int k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            for (int f = 0; f < kNumFields; ++f) {
              localFields(f, localIdx) = data.getCurrentField(f)(globalIdx);
              localWork(f, localIdx) = 0.0f;
            }
          });
          team.team_barrier();

          real_t const* W_ptr = W_global.data() + static_cast<size_t>(elementNumber) * kStride;
          INTEGRAL_TYPE::computeStiffnessOperatorTeamVector(team, &localFields(0, 0), &localWork(0, 0), W_ptr,
                                                            D_arr.data());

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pointsPerElem), [&](const int localIdx) {
            int i = localIdx % dim;
            int j = (localIdx / dim) % dim;
            int k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork(f, localIdx));
            }
          });
        });
  } else {
    throw std::runtime_error(
        "computeElementContributions_Acoustic_Gemm: INTEGRAL_TYPE has no GEMM path "
        "(use Qk_Hexahedron_Tensorial).");
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic(
    const DataType& data) {
  if constexpr (detail::has_team_gemm<INTEGRAL_TYPE>::value) {
    computeElementContributions_Acoustic_Gemm(data);
  } else {
    if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
      prepareZDeformedGeometry();
      if (zDeformedEnabled_) {
        computeElementContributions_Acoustic_TeamZ(data);
        return;
      }
    }
    computeElementContributions_Acoustic_Flat(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::prepareZDeformedGeometry() {
  if (zDeformedReady_) return;
  zDeformedReady_ = true;
  zDeformedEnabled_ = false;
  zDeformedGeom_ = vectorReal();
  zDeformedInvRho_ = vectorReal();
  zDeformedNodes_ = vectorInt();
  zDeformedNodesList_ = vectorInt();
  zDeformedInvRhoNodes_ = vectorReal();
  zDeformedBasisTab_ = vectorReal();

  bool const is_acoustic = (PHYSICS == utils::enums::physicType::kAcoustic);

  auto mesh_pc = m_mesh;
  int const nElems = mesh_pc.getNumberOfElements();
  vectorReal geom = allocateVector<vectorReal>(nElems * kZGeomStride, "zDeformedGeom");
  vectorReal invRho;
  if (is_acoustic && IS_MODEL_ON_NODES)
    invRho = allocateVector<vectorReal>(nElems * kPointsPerElement, "zDeformedInvRho");
  float* geomPtr = geom.data();
  float* invRhoPtr = invRho.data();

  int nonConforming = 0;
  Kokkos::parallel_reduce(
      "ZDeformed Geometry", Kokkos::RangePolicy<>(0, nElems),
      KOKKOS_LAMBDA(const int e, int& bad) {
        constexpr int dim = ORDER + 1;
        float X[8][3];
        auto const eIdx = mesh_pc.elementIndex(e);
        int I = 0;
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv) mesh_pc.vertexCoords(mesh_pc.globalVertexIndex(eIdx, iv, jv, kv), X[I++]);

        // Vertex k has indices (k % 2, (k % 4) / 2, k / 4): x may follow the first only, y the second only.
        float const hx = X[1][0] - X[0][0];
        float const hy = X[2][1] - X[0][1];
        float const tolX = 1e-4f * Kokkos::fabs(hx);
        float const tolY = 1e-4f * Kokkos::fabs(hy);
        for (int k = 0; k < 8; ++k) {
          if (Kokkos::fabs(X[k][0] - X[k % 2][0]) > tolX || Kokkos::fabs(X[k][1] - X[2 * ((k % 4) / 2)][1]) > tolY) {
            ++bad;
            break;
          }
        }

        geomPtr[e] = 0.5f * hx;
        geomPtr[nElems + e] = 0.5f * hy;
        for (int k = 0; k < 8; ++k) geomPtr[(2 + k) * nElems + e] = X[k][2];

        if (!is_acoustic) {
          geomPtr[(kZGeomStride - 1) * nElems + e] = 0.0f;
        } else if (IS_MODEL_ON_NODES) {
          geomPtr[(kZGeomStride - 1) * nElems + e] = 0.0f;
          for (int q = 0; q < kPointsPerElement; ++q) {
            int const g = mesh_pc.globalNodeIndex(e, q % dim, (q / dim) % dim, q / (dim * dim));
            invRhoPtr[q * nElems + e] = 1.0f / mesh_pc.getModelRhoOnNodes(g);
          }
        } else {
          geomPtr[(kZGeomStride - 1) * nElems + e] = 1.0f / mesh_pc.getModelRhoOnElement(e);
        }
      },
      nonConforming);

  if (nonConforming > 0) return;

  zDeformedGeom_ = geom;
  zDeformedInvRho_ = invRho;
  zDeformedEnabled_ = true;
  if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    vectorReal basisTab = allocateVector<vectorReal>(INTEGRAL_TYPE::kBasisTableSize, "zDeformedBasisTab");
    auto h_basisTab = Kokkos::create_mirror_view(basisTab);
    for (int idx = 0; idx < INTEGRAL_TYPE::kBasisTableSize; ++idx) h_basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx);
    Kokkos::deep_copy(basisTab, h_basisTab);
    zDeformedBasisTab_ = basisTab;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
int const* SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::zDeformedNodeTable(
    bool element_major) {
  bool const list_on = m_list_mode_;
  int const n_visited = list_on ? m_n_elem_list_ : m_mesh.getNumberOfElements();
  int const* list_key = list_on ? m_elem_list_.data() : nullptr;
  if (zDeformedNodes_.extent(0) == static_cast<size_t>(n_visited) * kPointsPerElement &&
      zDeformedNodesList_.data() == list_key && zDeformedNodesElementMajor_ == element_major)
    return zDeformedNodes_.data();

  auto mesh_pc = m_mesh;
  auto list_local = m_elem_list_;
  // Filled on the device, so that the kernels never read the mesh connectivity, whose pages the
  // host may pull back under unified memory.
  vectorInt nodes = allocateVector<vectorInt>(n_visited * kPointsPerElement, "zDeformedNodes");
  int* nodesPtr = nodes.data();
  Kokkos::parallel_for(
      "ZDeformed Node Table", Kokkos::RangePolicy<>(0, n_visited), KOKKOS_LAMBDA(const int i) {
        constexpr int dim = ORDER + 1;
        int const e = list_on ? list_local[i] : i;
        // The team kernels read the indices element by element, the flat kernel node by node.
        for (int q = 0; q < kPointsPerElement; ++q)
          nodesPtr[element_major ? i * kPointsPerElement + q : q * n_visited + i] =
              mesh_pc.globalNodeIndex(e, q % dim, (q / dim) % dim, q / (dim * dim));
      });

  zDeformedNodes_ = nodes;
  zDeformedNodesList_ = list_on ? m_elem_list_ : vectorInt();
  zDeformedNodesElementMajor_ = element_major;
  return zDeformedNodes_.data();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
float const* SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::zDeformedInvRhoNodeTable() {
  if constexpr (IS_MODEL_ON_NODES) {
    int const n_nodes = m_mesh.getNumberOfNodes();
    if (zDeformedInvRhoNodes_.extent(0) == static_cast<size_t>(n_nodes)) return zDeformedInvRhoNodes_.data();
    auto mesh_pc = m_mesh;
    vectorReal table = allocateVector<vectorReal>(n_nodes, "zDeformedInvRhoNodes");
    float* tablePtr = table.data();
    Kokkos::parallel_for(
        "ZDeformed InvRho Node Table", Kokkos::RangePolicy<>(0, n_nodes),
        KOKKOS_LAMBDA(const int g) { tablePtr[g] = 1.0f / mesh_pc.getModelRhoOnNodes(g); });
    zDeformedInvRhoNodes_ = table;
  }
  return zDeformedInvRhoNodes_.data();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_FlatZ(
    const DataType& data) {
  if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    auto mesh_local = m_mesh;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = mesh_local.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto geom = zDeformedGeom_;
    float const* invRhoNodes = zDeformedInvRho_.data();
    int const* elemNodes = zDeformedNodeTable();
    auto force = workVectorsGlobal_[0];

    // Input and output arrays of one element live in registers: 3 blocks per SM leave room for both.
    using Policy = Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, 3>>;

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic FlatZ", detail::lightWeight(Policy(0, n_iter)), KOKKOS_LAMBDA(const int _loop_idx) {
          int const e = list_on ? list_local[_loop_idx] : _loop_idx;
          constexpr int dim = ORDER + 1;

          float u[kPointsPerElement];
          for (int k = 0; k < dim; ++k)
            for (int j = 0; j < dim; ++j)
              for (int i = 0; i < dim; ++i)
                u[i + j * dim + k * dim * dim] =
                    data.getCurrentField(0)(elemNodes[(i + j * dim + k * dim * dim) * n_iter + _loop_idx]);

          real_t const J00 = geom(e);
          real_t const J11 = geom(nElems + e);
          real_t Z[8];
          for (int k = 0; k < 8; ++k) Z[k] = geom((2 + k) * nElems + e);
          real_t const invRhoElem = geom((kZGeomStride - 1) * nElems + e);

          auto get_alpha = [&](const int q) -> real_t {
            if (IS_MODEL_ON_NODES) return invRhoNodes[q * nElems + e];
            return invRhoElem;
          };

          float v[kPointsPerElement];
          INTEGRAL_TYPE::computeStiffnessTermSumFactZDeformedScatter(J00, J11, Z, u, v, get_alpha);

          for (int k = 0; k < dim; ++k)
            for (int j = 0; j < dim; ++j)
              for (int i = 0; i < dim; ++i)
                ATOMICADD(force[elemNodes[(i + j * dim + k * dim * dim) * n_iter + _loop_idx]],
                          v[i + j * dim + k * dim * dim]);
        });
  } else {
    throw std::runtime_error(
        "computeElementContributions_Acoustic_FlatZ: INTEGRAL_TYPE has no z-deformed sum-factorization path.");
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_TeamZ(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kAcoustic || !detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    throw std::runtime_error("computeElementContributions_Acoustic_TeamZ: needs an acoustic z-deformed setup.");
  } else {
    constexpr int n = ORDER + 1;
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    constexpr int kMinTeamsPerSM = std::clamp(1024 / kPreferredTeamSize, 1, 16);

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<kPreferredTeamSize, kMinTeamsPerSM>>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    // One thread per node and one per geometry entry are required. The last geometry entry, 1/rho of
    // the element, is only read when the model lives on elements.
    constexpr int kGeomLoaded = IS_MODEL_ON_NODES ? kZGeomStride - 1 : kZGeomStride;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    if (team_size < std::max(kPointsPerElement, kGeomLoaded)) {
      computeElementContributions_Acoustic_FlatZ(data);
      return;
    }

    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = m_mesh.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto zgeom = zDeformedGeom_;
    int const* elemNodes = zDeformedNodeTable(true);
    float const* invRhoNodes = zDeformedInvRhoNodeTable();
    float const* basisTabGlobal = zDeformedBasisTab_.data();
    auto force = workVectorsGlobal_[0];

    // Same choice as the TTI team kernel: up to 8 elements per team, at least four waves of teams.
    constexpr int kMaxElemsPerTeam = 8;
    constexpr int kMinWaves = 4;
    int const resident_teams = std::max<int>(1, ExecSpace::concurrency() / kPreferredTeamSize);
    int const elems_per_team = std::clamp(n_iter / (kMinWaves * resident_teams), 1, kMaxElemsPerTeam);
    TeamPolicyType policy((n_iter + elems_per_team - 1) / elems_per_team, team_size);

    // Nodal values, then the three fluxes, then the basis table, then the element geometry.
    constexpr int kFluxOffset = kPointsPerElement;
    constexpr int kTabOffset = kFluxOffset + 3 * kPointsPerElement;
    constexpr int kGeomOffset = kTabOffset + INTEGRAL_TYPE::kBasisTableSize;
    constexpr int kScratchFloats = kGeomOffset + kZGeomStride;
    policy.set_scratch_size(0, Kokkos::PerTeam(ScratchView1D::shmem_size(kScratchFloats)));

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic TeamZ", detail::lightWeight(policy),
        KOKKOS_LAMBDA(const TeamMember& team) {
          ScratchView1D scratch(team.team_scratch(0), kScratchFloats);
          // Captured here rather than first inside an if constexpr branch, which nvcc rejects.
          float const* inv_rho_nodes = invRhoNodes;
          auto const geom = zgeom;
          float* u = scratch.data();
          float* G = scratch.data() + kFluxOffset;
          float* tab = scratch.data() + kTabOffset;
          float* geom_sh = scratch.data() + kGeomOffset;

          // Thread q owns point q, which is also node q.
          int const q = team.team_rank();
          bool const owns = q < kPointsPerElement;
          int const qa = q % n;
          int const qb = (q / n) % n;
          int const qc = q / (n * n);

          // Filled before the first barrier, which covers it.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { tab[idx] = basisTabGlobal[idx]; });

          // The same for every element of the team: basis derivatives along the three lines through q,
          // the quadrature weight and the interpolation coordinates of q.
          float d_xi[n];
          float d_eta[n];
          float d_zeta[n];
          float weight = 0.0f;
          float alpha[3] = {0.0f, 0.0f, 0.0f};
          if (owns) {
            for (int i = 0; i < n; ++i) {
              d_xi[i] = basisTabGlobal[i * n + qa];
              d_eta[i] = basisTabGlobal[i * n + qb];
              d_zeta[i] = basisTabGlobal[i * n + qc];
            }
            weight = basisTabGlobal[n * n + qa] * basisTabGlobal[n * n + qb] * basisTabGlobal[n * n + qc];
            alpha[0] = basisTabGlobal[n * n + n + qa];
            alpha[1] = basisTabGlobal[n * n + n + qb];
            alpha[2] = basisTabGlobal[n * n + n + qc];
          }

          int const first_rank = team.league_rank() * elems_per_team;
          int const n_here = Kokkos::min(elems_per_team, n_iter - first_rank);

          // Two-stage pipeline held in registers: the node index and element number two elements
          // ahead, so that loading the next element's values never waits on its index. Each of the
          // first kGeomLoaded threads also carries one geometry entry of the next element.
          int ahead_node = 0;
          int ahead_e = 0;
          int next_node = 0;
          float next_u = 0.0f;
          float next_inv_rho = 0.0f;
          float next_geom = 0.0f;
          auto const load_index = [&](int rank) {
            if (owns) ahead_node = elemNodes[rank * kPointsPerElement + q];
            if (q < kGeomLoaded) ahead_e = list_on ? list_local[rank] : rank;
          };
          auto const load_values = [&]() {
            if (owns) {
              next_node = ahead_node;
              next_u = data.getCurrentField(0)(next_node);
              if constexpr (IS_MODEL_ON_NODES) next_inv_rho = inv_rho_nodes[next_node];
            }
            if (q < kGeomLoaded) next_geom = geom(q * nElems + ahead_e);
          };
          load_index(first_rank);
          load_values();
          if (n_here > 1) load_index(first_rank + 1);

          for (int k = 0; k < n_here; ++k) {
            int const node = next_node;
            float inv_rho = next_inv_rho;

            if (owns) u[q] = next_u;
            if (q < kGeomLoaded) geom_sh[q] = next_geom;
            team.team_barrier();

            if (k + 1 < n_here) load_values();
            if (k + 2 < n_here) load_index(first_rank + k + 2);

            float const J00 = geom_sh[0];
            float const J11 = geom_sh[1];
            float Z[8];
            for (int j = 0; j < 8; ++j) Z[j] = geom_sh[2 + j];
            if constexpr (!IS_MODEL_ON_NODES) inv_rho = geom_sh[kZGeomStride - 1];

            if (owns) {
              float du[3] = {0.0f, 0.0f, 0.0f};
              for (int i = 0; i < n; ++i) {
                du[0] += d_xi[i] * u[i + qb * n + qc * n * n];
                du[1] += d_eta[i] * u[qa + i * n + qc * n * n];
                du[2] += d_zeta[i] * u[qa + qb * n + i * n * n];
              }

              float dZ[3][4];
              INTEGRAL_TYPE::zDeformedHalfDifferences(Z, dZ);
              float invJ[3][3];
              float const detJ = INTEGRAL_TYPE::invJacobianZDeformedBilinear(alpha, J00, J11, dZ, invJ);
              float const a = invJ[0][0];
              float const b = invJ[1][1];
              float const c = invJ[2][0];
              float const d = invJ[2][1];
              float const ez = invJ[2][2];
              // B = det(J) J^-1 J^-T, whose xy entry is zero.
              float const scale = weight * inv_rho * detJ;
              float const B0 = a * a;
              float const B1 = b * b;
              float const B2 = c * c + d * d + ez * ez;
              float const B3 = b * d;
              float const B4 = a * c;
              G[q] = scale * (B0 * du[0] + B4 * du[2]);
              G[kPointsPerElement + q] = scale * (B1 * du[1] + B3 * du[2]);
              G[2 * kPointsPerElement + q] = scale * (B4 * du[0] + B3 * du[1] + B2 * du[2]);
            }
            team.team_barrier();

            // Transpose: node (qa, qb, qc) gathers D(qa, p) G_xi(p, qb, qc), and so on.
            if (owns) {
              float v = 0.0f;
              for (int p = 0; p < n; ++p) {
                v += tab[qa * n + p] * G[p + qb * n + qc * n * n];
                v += tab[qb * n + p] * G[kPointsPerElement + qa + p * n + qc * n * n];
                v += tab[qc * n + p] * G[2 * kPointsPerElement + qa + qb * n + p * n * n];
              }
              ATOMICADD(force[node], v);
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  using Policy = Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>;

  Kokkos::parallel_for(
      "Solver Element Contribution Acoustic", Policy(0, n_iter), KOKKOS_LAMBDA(const int _loop_idx) {
        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        constexpr int dim = ORDER + 1;

        float localFields[kNumFields][kPointsPerElement];
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int k = 0; k < dim; ++k) {
          for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        real_t inv_density = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density = 1.0f / mesh_local.getModelRhoOnElement(elementNumber);
        }

        INTEGRAL_TYPE::computeStiffnessTermSumFact(
            cornerCoords, localFields[0], localWork[0], [&](const int qa, const int qb, const int qc) -> real_t {
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                return 1.0f / mesh_local.getModelRhoOnNodes(gIndex);
              } else {
                return inv_density;
              }
            });

        for (int k = 0; k < dim; ++k) {
          for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) {
                ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributions(
    const DataType& data) {
  if (!attenuationEnabled_ || nSls_ <= 0) return;

  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic)
    computeAttenuationContributionsAcoustic(data);
  else
    computeAttenuationContributionsElastic(data);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsAcoustic(
    const DataType& data) {
  auto mesh_local = m_mesh;
  Kokkos::Array<vectorReal, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Acoustic",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWorkA[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
            }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        real_t inv_density_q = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density_q =
              1.0f / (mesh_local.getModelRhoOnElement(elementNumber) * mesh_local.getModelQpOnElement(elementNumber));
        }

        INTEGRAL_TYPE::computeStiffnessTerm(
            cornerCoords,
            [&](const int qa, const int qb, const int qc) {
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                inv_density_q = 1.0f / (mesh_local.getModelRhoOnNodes(gIndex) * mesh_local.getModelQpOnNodes(gIndex));
              }
            },
            [&](const int i, const int j, const real_t val) {
              localWorkA[0][i] += inv_density_q * val * localFields[0][j];
            });

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
                ATOMICADD(local_attenuationWorkVectorsGlobal[f][globalIdx], localWorkA[f][localIdx]);
            }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsElastic(
    const DataType& data) {
  if (anisotropyType_ != model::AnisotropyType::kIso) return;

  auto mesh_local = m_mesh;

  Kokkos::Array<vectorReal, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Elastic",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        constexpr int dim = ORDER + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWorkA[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
            }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        // Per (p, r) pair of reference directions: the six distinct coefficients of the
        // isotropic stiffness contracted with the inverse Jacobian.
        // a0, a1, a2 multiply the diagonal displacement terms; a3, b0, b1 the xy, xz, yz couplings.
        struct CJPacked {
          float a0, a1, a2, a3;
          float b0, b1;
        };
        CJPacked CJflat[3 * 3];

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            cornerCoords,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              float vp, vs, rho, qp, qs;
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                vp = mesh_local.getModelVpOnNodes(gIndex);
                vs = mesh_local.getModelVsOnNodes(gIndex);
                rho = mesh_local.getModelRhoOnNodes(gIndex);
                qp = mesh_local.getModelQpOnNodes(gIndex);
                qs = mesh_local.getModelQsOnNodes(gIndex);
              } else {
                vp = mesh_local.getModelVpOnElement(elementNumber);
                vs = mesh_local.getModelVsOnElement(elementNumber);
                rho = mesh_local.getModelRhoOnElement(elementNumber);
                qp = mesh_local.getModelQpOnElement(elementNumber);
                qs = mesh_local.getModelQsOnElement(elementNumber);
              }

              float mu = rho * vs * vs;
              float lambda = rho * (vp * vp - 2.0f * vs * vs);
              float lambdap2mua = (lambda + 2.0f * mu) / qp;
              mu = mu / qs;
              lambda = lambdap2mua - 2.0f * mu;
              float lambda_plus_2mu = lambda + 2.0f * mu;

              for (int p = 0; p < 3; ++p) {
                float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
                for (int r = 0; r < 3; ++r) {
                  float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                  int const idx = p * 3 + r;
                  CJflat[idx].a0 = lambda_plus_2mu * Jp0 * Jr0 + mu * (Jp1 * Jr1 + Jp2 * Jr2);
                  CJflat[idx].a1 = mu * Jp0 * Jr0 + lambda_plus_2mu * Jp1 * Jr1 + mu * Jp2 * Jr2;
                  CJflat[idx].a2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) + lambda_plus_2mu * Jp2 * Jr2;
                  CJflat[idx].a3 = lambda * Jp0 * Jr1 + mu * Jp1 * Jr0;
                  CJflat[idx].b0 = lambda * Jp0 * Jr2 + mu * Jp2 * Jr0;
                  CJflat[idx].b1 = lambda * Jp1 * Jr2 + mu * Jp2 * Jr1;
                }
              }
            },
            [&](int i, int j, float val, const int p, const int r) {
              int const idx = p * 3 + r;
              float const uxj = localFields[0][j];
              float const uyj = localFields[1][j];
              float const uzj = localFields[2][j];
              localWorkA[0][i] += val * (CJflat[idx].a0 * uxj + CJflat[idx].a3 * uyj + CJflat[idx].b0 * uzj);
              localWorkA[1][i] += val * (CJflat[idx].a3 * uxj + CJflat[idx].a1 * uyj + CJflat[idx].b1 * uzj);
              localWorkA[2][i] += val * (CJflat[idx].b0 * uxj + CJflat[idx].b1 * uyj + CJflat[idx].a2 * uzj);
            });

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
                ATOMICADD(local_attenuationWorkVectorsGlobal[f][globalIdx], localWorkA[f][localIdx]);
            }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Element Contribution Iso Flat",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;

        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float mu_e = 0.0f, lambda_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
          }

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork,
                                                        [&](int qa, int qb, int qc, float const(&J_inv)[3][3],
                                                            float const(&grad_u_ref)[3][3], float(&flux)[3][3]) {
                                                          float mu, lambda;
                                                          if constexpr (IS_MODEL_ON_NODES) {
                                                            int const gIndex =
                                                                mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                                                            float const vp = mesh_local.getModelVpOnNodes(gIndex);
                                                            float const vs = mesh_local.getModelVsOnNodes(gIndex);
                                                            float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                                                            mu = rho * vs * vs;
                                                            lambda = rho * (vp * vp - 2.0f * vs * vs);
                                                          } else {
                                                            mu = mu_e;
                                                            lambda = lambda_e;
                                                          }

                                                          flux::elasticFluxIso(J_inv, mu, lambda, grad_u_ref, flux);
                                                        });

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Iso_Flat(data);
  } else {
    computeElementContributions_Iso_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso_Team(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  // avoid extended __host__ __device__ lambda cannot first-capture variable in
  // constexpr-if context
  (void)mesh_local;
  (void)list_on;
  (void)list_local;
  (void)n_iter;
  (void)local_workVectorsGlobal;

  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;

    // A team owns one element: gather and scatter are spread over its threads,
    // which makes them coalesced, and the per-element buffers move from
    // per-thread local memory to shared memory.
    //
    // Every team range covers kPointsPerElement, so a larger team only adds
    // threads that idle at the barriers. Round up to a whole warp so no partial
    // warp idles where kPointsPerElement is not a multiple of 32, then clamp to
    // what the backend can actually provide: a host backend caps the team size
    // at its thread count, and asking for more aborts at launch. A
    // TeamThreadRange is correct at any team size.
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    // Constant-Jacobian meshes carry the element geometry (9 inverse-Jacobian
    // entries + the determinant) in a small per-team scratch buffer.
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;
    TeamPolicyType policy(n_iter, team_size);
    // One buffer carries the displacements in and the forces out: the stiffness
    // kernel overwrites it after the barrier that ends its read phase. Halving
    // this buffer is what lifts the occupancy off its shared-memory limit.
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Iso Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float mu_e = 0.0f, lambda_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
          }

          auto const iso_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            float mu, lambda;
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float const vp = mesh_local.getModelVpOnNodes(gIndex);
              float const vs = mesh_local.getModelVsOnNodes(gIndex);
              float const rho = mesh_local.getModelRhoOnNodes(gIndex);
              mu = rho * vs * vs;
              lambda = rho * (vp * vp - 2.0f * vs * vs);
            } else {
              mu = mu_e;
              lambda = lambda_e;
            }
            flux::elasticFluxIso(J_inv, mu, lambda, grad_u_ref, flux);
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              // invJacobianTransformation accumulates into J, so it must start at zero.
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), iso_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), iso_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Element Contribution Vti Flat",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;
        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float c11_e = 0, c12_e = 0, c13_e = 0, c33_e = 0, c44_e = 0, c66_e = 0;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            float const delta_e = mesh_local.getModelDeltaOnElement(elementNumber);
            float const epsilon_e = mesh_local.getModelEpsilonOnElement(elementNumber);
            float const gamma_e = mesh_local.getModelGammaOnElement(elementNumber);
            float const rho_vp2 = rho_e * vp_e * vp_e;
            float const rho_vs2 = rho_e * vs_e * vs_e;
            c33_e = rho_vp2;
            c44_e = rho_vs2;
            c11_e = rho_vp2 * (1.0f + 2.0f * epsilon_e);
            c66_e = rho_vs2 * (1.0f + 2.0f * gamma_e);
            float const vp2_vs2 = vp_e * vp_e - vs_e * vs_e;
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp_e * vp_e * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(
              cornerCoords, localFields, localWork,
              [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                  float(&flux)[3][3]) {
                float c11, c12, c13, c33, c44, c66;
                if constexpr (IS_MODEL_ON_NODES) {
                  int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  float const vp = mesh_local.getModelVpOnNodes(gIndex);
                  float const vs = mesh_local.getModelVsOnNodes(gIndex);
                  float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                  float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
                  float const epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
                  float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
                  float const rho_vp2 = rho * vp * vp;
                  float const rho_vs2 = rho * vs * vs;
                  c33 = rho_vp2;
                  c44 = rho_vs2;
                  c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
                  c66 = rho_vs2 * (1.0f + 2.0f * gamma);
                  float const vp2_vs2 = vp * vp - vs * vs;
                  c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;
                  c12 = c11 - 2.0f * c66;
                } else {
                  c11 = c11_e;
                  c12 = c12_e;
                  c13 = c13_e;
                  c33 = c33_e;
                  c44 = c44_e;
                  c66 = c66_e;
                }

                flux::elasticFluxVti(J_inv, c11, c12, c13, c33, c44, c66, grad_u_ref, flux);
              });

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Vti_Flat(data);
  } else {
    computeElementContributions_Vti_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti_Team(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  // avoid extended __host__ __device__ lambda cannot first-capture variable in
  // constexpr-if context
  (void)mesh_local;
  (void)list_on;
  (void)list_local;
  (void)local_workVectorsGlobal;

  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;
    // Warp-aligned preferred size, clamped to what the backend can provide (a
    // host backend caps it at its thread count; asking for more aborts).
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;

    TeamPolicyType policy(n_iter, team_size);
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Vti Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float c11_e = 0, c12_e = 0, c13_e = 0, c33_e = 0, c44_e = 0, c66_e = 0;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            float const delta_e = mesh_local.getModelDeltaOnElement(elementNumber);
            float const epsilon_e = mesh_local.getModelEpsilonOnElement(elementNumber);
            float const gamma_e = mesh_local.getModelGammaOnElement(elementNumber);
            float const rho_vp2 = rho_e * vp_e * vp_e;
            float const rho_vs2 = rho_e * vs_e * vs_e;
            c33_e = rho_vp2;
            c44_e = rho_vs2;
            c11_e = rho_vp2 * (1.0f + 2.0f * epsilon_e);
            c66_e = rho_vs2 * (1.0f + 2.0f * gamma_e);
            float const vp2_vs2 = vp_e * vp_e - vs_e * vs_e;
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp_e * vp_e * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          auto const vti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            float c11, c12, c13, c33, c44, c66;
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float const vp = mesh_local.getModelVpOnNodes(gIndex);
              float const vs = mesh_local.getModelVsOnNodes(gIndex);
              float const rho = mesh_local.getModelRhoOnNodes(gIndex);
              float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
              float const epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
              float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
              float const rho_vp2 = rho * vp * vp;
              float const rho_vs2 = rho * vs * vs;
              c33 = rho_vp2;
              c44 = rho_vs2;
              c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
              c66 = rho_vs2 * (1.0f + 2.0f * gamma);
              float const vp2_vs2 = vp * vp - vs * vs;
              c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;
              c12 = c11 - 2.0f * c66;
            } else {
              c11 = c11_e;
              c12 = c12_e;
              c13 = c13_e;
              c33 = c33_e;
              c44 = c44_e;
              c66 = c66_e;
            }
            flux::elasticFluxVti(J_inv, c11, c12, c13, c33, c44, c66, grad_u_ref, flux);
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), vti_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), vti_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

// Stores the compact TTI description of every node (see computeTtiCompact()) into cttiNodes_.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::precomputeTtiTensorsOnNodes() {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic || !IS_MODEL_ON_NODES) {
    return;
  } else {
    if (cttiNodesReady_) {
      return;
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    int const numNodes = m_mesh.getNumberOfNodes();
    cttiNodes_ = allocateArray2D<CttiView>(numNodes, flux::kTtiCompactSize, "cttiNodes");

    auto ctti = cttiNodes_;
    auto mesh_pc = m_mesh;
    Kokkos::parallel_for(
        "Tti Precompute Compact On Nodes", Kokkos::RangePolicy<ExecSpace>(0, numNodes), KOKKOS_LAMBDA(const int g) {
          float p[flux::kTtiCompactSize];
          computeTtiCompact(mesh_pc.getModelVpOnNodes(g), mesh_pc.getModelVsOnNodes(g), mesh_pc.getModelRhoOnNodes(g),
                            mesh_pc.getModelDeltaOnNodes(g), mesh_pc.getModelEpsilonOnNodes(g),
                            mesh_pc.getModelGammaOnNodes(g), mesh_pc.getModelPhiOnNodes(g),
                            mesh_pc.getModelThetaOnNodes(g), p);
          for (int k = 0; k < flux::kTtiCompactSize; ++k) ctti(g, k) = p[k];
        });
    Kokkos::fence();

    cttiNodesReady_ = true;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_Flat(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    precomputeTtiTensorsOnNodes();

    auto mesh_local = m_mesh;
    [[maybe_unused]] auto ctti_local = cttiNodes_;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    Kokkos::parallel_for(
        "Solver Element Contribution Tti Flat",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
        KOKKOS_LAMBDA(const int _loop_idx) {
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          int const dim = mesh_local.getOrder() + 1;
          float localFields[kNumFields][kPointsPerElement] = {{0}};
          float localWork[kNumFields][kPointsPerElement] = {{0}};

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
                }
              }
            }
          }

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          // Captured here, not at first use: nvcc rejects an extended lambda that first-captures
          // a variable inside an if-constexpr branch.
          [[maybe_unused]] auto const ctti = ctti_local;

          // Model on elements only: one tensor for the whole element.
          [[maybe_unused]] float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          auto ttiFlux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                             float(&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
              float p[flux::kTtiCompactSize];
              for (int k = 0; k < flux::kTtiCompactSize; ++k) p[k] = ctti(gIndex, k);
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            } else {
              flux::elasticFluxTti(J_inv, CTTI, grad_u_ref, flux);
            }
          };

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork, ttiFlux);

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Tti_Flat(data);
  } else {
    computeElementContributions_Tti_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_Team(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    if constexpr (!HasConstantJacobian<MESH_TYPE>::value && IS_MODEL_ON_NODES &&
                  detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
      prepareZDeformedGeometry();
      if (zDeformedEnabled_) {
        precomputeTtiTensorsOnNodes();
        computeElementContributions_Tti_TeamZ(data);
        return;
      }
    }
    precomputeTtiTensorsOnNodes();

    auto mesh_local = m_mesh;
    [[maybe_unused]] auto ctti_local = cttiNodes_;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;
    // Warp-aligned preferred size, clamped to what the backend can provide (a
    // host backend caps it at its thread count; asking for more aborts).
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;

    TeamPolicyType policy(n_iter, team_size);
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Tti Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          // Captured here, not at first use: nvcc rejects an extended lambda that first-captures
          // a variable inside an if-constexpr branch.
          [[maybe_unused]] auto const ctti = ctti_local;

          // Model on elements only: one tensor for the whole element.
          [[maybe_unused]] float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          auto const tti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float p[flux::kTtiCompactSize];
              for (int k = 0; k < flux::kTtiCompactSize; ++k) p[k] = ctti(gIndex, k);
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            } else {
              flux::elasticFluxTti(J_inv, CTTI, grad_u_ref, flux);
            }
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), tti_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), tti_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_TeamZ(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic || !IS_MODEL_ON_NODES ||
                !detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    throw std::runtime_error("computeElementContributions_Tti_TeamZ: needs an elastic on-nodes z-deformed setup.");
  } else {
    auto ctti_local = cttiNodes_;
    // Two 16-byte loads per node instead of eight scalar ones: each scalar load of a warp
    // touches the same cache lines again.
    static_assert(flux::kTtiCompactSize == 8, "two 16-byte loads per node");
    if (ctti_local.stride(0) != flux::kTtiCompactSize ||
        reinterpret_cast<std::uintptr_t>(ctti_local.data()) % alignof(detail::Float4) != 0)
      throw std::runtime_error("computeElementContributions_Tti_TeamZ: cttiNodes_ not packed for 16-byte loads.");
    detail::Float4 const* ctti_quads = reinterpret_cast<detail::Float4 const*>(ctti_local.data());
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = m_mesh.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto zgeom = zDeformedGeom_;
    int const* elemNodes = zDeformedNodeTable();
    float const* basisTabGlobal = zDeformedBasisTab_.data();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    constexpr int dim = ORDER + 1;
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    // Enough resident teams for half the maximum warps per SM, which caps registers at 64 per thread
    // without spilling.
    // Capped at 16, the lowest resident-block limit per SM among the targeted GPUs (24 on Ada, 16 on GA10x).
    constexpr int kMinTeamsPerSM = std::clamp(1024 / kPreferredTeamSize, 1, 16);

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<kPreferredTeamSize, kMinTeamsPerSM>>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    // With one thread per node, a team handles several consecutive elements and loads the next one's fields
    // into registers while it computes the current one, so the gather latency overlaps the computation.
    // Fewer elements per team on small meshes, so that the league still fills the device several times over.
    constexpr int kMaxElemsPerTeam = 8;
    constexpr int kMinWaves = 4;
    int const resident_teams = std::max<int>(1, ExecSpace::concurrency() / kPreferredTeamSize);
    bool const pipelined = team_size >= kPointsPerElement;
    int const elems_per_team =
        pipelined ? std::clamp(n_iter / (kMinWaves * resident_teams), 1, kMaxElemsPerTeam) : 1;
    TeamPolicyType policy((n_iter + elems_per_team - 1) / elems_per_team, team_size);
    // One float block: fields (then forces), fluxes and basis table.
    constexpr int kFluxOffset = kNumFields * kPointsPerElement;
    constexpr int kTabOffset = kFluxOffset + 9 * kPointsPerElement;
    constexpr int kScratchFloats = kTabOffset + INTEGRAL_TYPE::kBasisTableSize;
    size_t const bytes = ScratchView1D::shmem_size(kScratchFloats) + ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes));

    Kokkos::parallel_for(
        "Solver Element Contribution Tti TeamZ", detail::lightWeight(policy), KOKKOS_LAMBDA(const TeamMember& team) {
          ScratchView1D scratch(team.team_scratch(0), kScratchFloats);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Displacements in, forces out, same storage.
          float* localFields = scratch.data();
          float* fluxScratch = scratch.data() + kFluxOffset;
          float* basisTab = scratch.data() + kTabOffset;

          int const first_rank = team.league_rank() * elems_per_team;
          int const n_here = Kokkos::min(elems_per_team, n_iter - first_rank);
          int const point = team.team_rank();
          bool const owns_point = pipelined && point < kPointsPerElement;

          // Next element's node index and fields, held in registers across the current element's computation.
          int next_node = 0;
          float next_fields[kNumFields] = {};
          auto const prefetch = [&](int rank) {
            if (owns_point) {
              // Indexed by the league position, so the gather does not wait for the element list.
              next_node = elemNodes[rank * kPointsPerElement + point];
              for (int f = 0; f < kNumFields; ++f) next_fields[f] = data.getCurrentField(f)(next_node);
            }
          };
          prefetch(first_rank);

          for (int k = 0; k < n_here; ++k) {
            int const rank = first_rank + k;
            int const elementNumber = list_on ? list_local[rank] : rank;

            if (k == 0) {
              // Filled before the gather barrier, so it needs no barrier of its own.
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                                   [&](const int idx) { basisTab[idx] = basisTabGlobal[idx]; });
            }

            if (pipelined) {
              if (owns_point) {
                nodeIdx(point) = next_node;
                for (int f = 0; f < kNumFields; ++f) localFields[f * kPointsPerElement + point] = next_fields[f];
              }
            } else {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
                int const globalIdx = elemNodes[rank * kPointsPerElement + localIdx];
                nodeIdx(localIdx) = globalIdx;
                for (int f = 0; f < kNumFields; ++f)
                  localFields[f * kPointsPerElement + localIdx] = data.getCurrentField(f)(globalIdx);
              });
            }
            team.team_barrier();

            if (k + 1 < n_here) prefetch(rank + 1);

            // Same address for the whole team: one broadcast load per entry.
            float Z[8];
            for (int j = 0; j < 8; ++j) Z[j] = zgeom((2 + j) * nElems + elementNumber);

            auto const tti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3],
                                      float const(&grad_u_ref)[3][3], float(&flux)[3][3]) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float p[flux::kTtiCompactSize];
              detail::Float4 const lo = ctti_quads[2 * static_cast<size_t>(gIndex)];
              detail::Float4 const hi = ctti_quads[2 * static_cast<size_t>(gIndex) + 1];
              for (int j = 0; j < 4; ++j) {
                p[j] = lo.v[j];
                p[j + 4] = hi.v[j];
              }
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            };

            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeamZDeformed(
                team, zgeom(elementNumber), zgeom(nElems + elementNumber), Z, localFields, localFields, fluxScratch,
                tti_flux, basisTab);
            team.team_barrier();

            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
              int const globalIdx = nodeIdx(localIdx);
              for (int f = 0; f < kNumFields; ++f) {
                ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields[f * kPointsPerElement + localIdx]);
              }
            });
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsForward(float dt,
                                                                                                 const DataType& data) {
  // Members are copied into locals so that the device lambdas capture values, not this.
  float const dt_local = dt;
  float const dt2_local = dt * dt;
  int const n_sls = nSls_;
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);

  auto mesh_local = m_mesh;
  auto mass_matrix = massMatrixGlobal_;
  auto taper_coeff = spongeTaperCoeff_;

  auto sls_w = slsReferenceAngularFrequencies_;
  auto sls_beta = slsAnelasticityCoefficients_;

  std::array<std::remove_reference_t<decltype(data.getCurrentField(0))>, kNumFields> current_field;
  std::array<std::remove_reference_t<decltype(data.getPreviousField(0))>, kNumFields> prev_field;
  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> damping_matrix;
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> work_vector;
  std::array<std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>, kNumFields> atten_work_vec;
  std::array<std::remove_reference_t<decltype(attenuationMemoryVariables_[0])>, kNumFields> atten_mem_vars;

  for (int f = 0; f < kNumFields; ++f) {
    current_field[f] = data.getCurrentField(f);
    prev_field[f] = data.getPreviousField(f);
    damping_matrix[f] = dampingMatrixGlobal_[f];
    work_vector[f] = workVectorsGlobal_[f];
    if (has_attenuation) {
      atten_work_vec[f] = attenuationWorkVectorsGlobal_[f];
      atten_mem_vars[f] = attenuationMemoryVariables_[f];
    }
  }

  bool const list_on = m_node_list_mode_;
  auto list_local = m_node_list_;
  // The forces are zeroed right after being read, which saves the reset at the start of the next step.
  bool const reset_forces = m_reset_forces_in_update_;

  // The new value is written into the previous-field buffer (leapfrog), which the caller swaps afterwards.
  // Every per-node value is loaded before the first test and the first store: the loads then go out
  // together instead of waiting for each other.
  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
    int const n_iter = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();
    // A node only needs 40 bytes, so each thread takes two to keep more loads in flight. The second node
    // is n_threads further, which keeps every load of a warp contiguous.
    constexpr int kNodesPerThread = 2;
    int const n_threads = (n_iter + kNodesPerThread - 1) / kNodesPerThread;
    Kokkos::parallel_for(
        "Solver Update Field Acoustic", detail::lightWeightRange(n_threads), KOKKOS_LAMBDA(const int _thread_idx) {
          int I[kNodesPerThread];
          bool active[kNodesPerThread];
          float mass[kNodesPerThread];
          bool free_surface[kNodesPerThread];
          float taper[kNodesPerThread];
          float cur[kNodesPerThread];
          float prev[kNodesPerThread];
          float damp[kNodesPerThread];
          float work[kNodesPerThread];
          for (int n = 0; n < kNodesPerThread; ++n) {
            int const node_idx = _thread_idx + n * n_threads;
            active[n] = node_idx < n_iter;
            I[n] = active[n] ? (list_on ? list_local[node_idx] : node_idx) : 0;
          }
          for (int n = 0; n < kNodesPerThread; ++n) {
            if (!active[n]) continue;
            mass[n] = mass_matrix(I[n]);
            free_surface[n] = mesh_local.isFreeSurface(I[n]);
            taper[n] = taper_coeff(I[n]);
            cur[n] = current_field[0](I[n]);
            prev[n] = prev_field[0](I[n]);
            damp[n] = damping_matrix[0](I[n]);
            work[n] = work_vector[0](I[n]);
          }

          for (int n = 0; n < kNodesPerThread; ++n) {
            if (!active[n]) continue;
            int const node = I[n];
            if (mass[n] <= 0.0f) {
            } else if (free_surface[n]) {
              current_field[0](node) = 0.0f;
              prev_field[0](node) = 0.0f;
            } else {
              float next_val =
                  (2.0f * mass[n] * cur[n] - (mass[n] - 0.5f * dt_local * damp[n]) * prev[n] - dt2_local * work[n]);

              if (has_attenuation) {
                for (int l = 0; l < n_sls; ++l) {
                  float const w = sls_w[l];
                  float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[0](node, l) + beta_p * atten_work_vec[0](node));

                  atten_mem_vars[0](node, l) = gamma * atten_mem_vars[0](node, l) + beta * atten_work_vec[0](node);
                }
              }

              prev_field[0](node) = next_val / (mass[n] + 0.5f * dt_local * damp[n]) * taper[n];
              current_field[0](node) = cur[n] * taper[n];
            }
            if (reset_forces) {
              work_vector[0](node) = 0.0f;
              if (has_attenuation) atten_work_vec[0](node) = 0.0f;
            }
          }
        });
  } else {
    int const n_iter_el = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();

    Kokkos::parallel_for(
        "Solver Update Field Elastic", detail::lightWeightRange(n_iter_el), KOKKOS_LAMBDA(const int _node_idx) {
          if (_node_idx >= n_iter_el) return;
          int const I = list_on ? list_local[_node_idx] : _node_idx;
          float const mass = mass_matrix(I);
          bool const free_surface = mesh_local.isFreeSurface(I);
          float const taper = taper_coeff(I);
          float cur[kNumFields];
          float prev[kNumFields];
          float damp[kNumFields];
          float work[kNumFields];
          for (int f = 0; f < kNumFields; ++f) {
            cur[f] = current_field[f](I);
            prev[f] = prev_field[f](I);
            damp[f] = damping_matrix[f](I);
            work[f] = work_vector[f](I);
          }

          // No early return: the stores all come after the loads.
          if (mass > 0.0f) {
            for (int f = 0; f < kNumFields; ++f) {
              // No damping on the free surface.
              float const d = free_surface ? 0.0f : damp[f];
              float next_val = (2.0f * mass * cur[f] - (mass - 0.5f * dt_local * d) * prev[f] - dt2_local * work[f]);

              if (has_attenuation) {
                for (int l = 0; l < n_sls; ++l) {
                  float const w = sls_w[l];
                  float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) + beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) + beta * atten_work_vec[f](I);
                }
              }

              prev_field[f](I) = next_val / (mass + 0.5f * dt_local * d) * taper;
              current_field[f](I) = cur[f] * taper;
            }
          }
          if (reset_forces) {
            for (int f = 0; f < kNumFields; ++f) {
              work_vector[f](I) = 0.0f;
              if (has_attenuation) atten_work_vec[f](I) = 0.0f;
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsBackward(
    float dt, const DataType& data) {
  // Members are copied into locals so that the device lambdas capture values, not this.
  float const dt_local = dt;
  float const dt2_local = dt * dt;
  int const n_sls = nSls_;
  bool const has_attenuation = (attenuationEnabled_ && nSls_ > 0);

  auto mesh_local = m_mesh;
  auto mass_matrix = massMatrixGlobal_;
  auto taper_coeff = spongeTaperCoeff_;

  auto sls_w = slsReferenceAngularFrequencies_;
  auto sls_beta = slsAnelasticityCoefficients_;

  std::array<std::remove_reference_t<decltype(data.getCurrentField(0))>, kNumFields> current_field;
  std::array<std::remove_reference_t<decltype(data.getPreviousField(0))>, kNumFields> prev_field;
  std::array<std::remove_reference_t<decltype(data.getPrevPrevField(0))>, kNumFields> prevprev_field;
  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> damping_matrix;
  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> work_vector;
  std::array<std::remove_reference_t<decltype(attenuationWorkVectorsGlobal_[0])>, kNumFields> atten_work_vec;
  std::array<std::remove_reference_t<decltype(attenuationMemoryVariables_[0])>, kNumFields> atten_mem_vars;

  for (int f = 0; f < kNumFields; ++f) {
    current_field[f] = data.getCurrentField(f);
    prev_field[f] = data.getPreviousField(f);
    prevprev_field[f] = data.getPrevPrevField(f);
    damping_matrix[f] = dampingMatrixGlobal_[f];
    work_vector[f] = workVectorsGlobal_[f];
    if (has_attenuation) {
      atten_work_vec[f] = attenuationWorkVectorsGlobal_[f];
      atten_mem_vars[f] = attenuationMemoryVariables_[f];
    }
  }

  bool const list_on = m_node_list_mode_;
  auto list_local = m_node_list_;

  // Same scheme as updateFieldsForward, but the new value goes to the third buffer (prevprev),
  // so the previous-field buffer is kept for the caller.
  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
    int const n_iter = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();
    Kokkos::parallel_for(
        "Solver Update Field Acoustic Backward", n_iter, KOKKOS_LAMBDA(const int _node_idx) {
          if (_node_idx >= n_iter) return;
          int const I = list_on ? list_local[_node_idx] : _node_idx;
          if (mass_matrix[I] <= 0.0f) return;

          if (mesh_local.isFreeSurface(I)) {
            current_field[0](I) = 0.0f;
            prev_field[0](I) = 0.0f;
            prevprev_field[0](I) = 0.0f;
          } else {
            float next_val = (2.0f * mass_matrix(I) * current_field[0](I) -
                              (mass_matrix(I) - 0.5f * dt_local * damping_matrix[0](I)) * prev_field[0](I) -
                              dt2_local * work_vector[0](I));

            if (has_attenuation) {
              for (int l = 0; l < n_sls; ++l) {
                float const w = sls_w[l];
                float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                float const gamma_p = 0.5f + 0.5f * gamma;
                float const beta_p = 0.5f * beta;

                next_val += dt2_local * (gamma_p * atten_mem_vars[0](I, l) + beta_p * atten_work_vec[0](I));

                atten_mem_vars[0](I, l) = gamma * atten_mem_vars[0](I, l) + beta * atten_work_vec[0](I);
              }
            }

            prevprev_field[0](I) = next_val / (mass_matrix(I) + 0.5f * dt_local * damping_matrix[0](I));
            prevprev_field[0](I) *= taper_coeff(I);
            current_field[0](I) *= taper_coeff(I);
          }
        });
  } else {
    int const n_iter_el = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();

    Kokkos::parallel_for(
        "Solver Update Field Elastic Backward", n_iter_el, KOKKOS_LAMBDA(const int _node_idx) {
          if (_node_idx >= n_iter_el) return;
          int const I = list_on ? list_local[_node_idx] : _node_idx;
          if (mass_matrix[I] <= 0.0f) return;
          if (mesh_local.isFreeSurface(I)) {
            for (int f = 0; f < kNumFields; ++f) {
              float next_val = (2.0f * mass_matrix(I) * current_field[f](I) - mass_matrix(I) * prev_field[f](I) -
                                dt2_local * work_vector[f](I));

              if (has_attenuation) {
                for (int l = 0; l < n_sls; ++l) {
                  float const w = sls_w[l];
                  float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) + beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) + beta * atten_work_vec[f](I);
                }
              }

              prevprev_field[f](I) = next_val / mass_matrix(I);
              prevprev_field[f](I) *= taper_coeff(I);
              current_field[f](I) *= taper_coeff(I);
            }
          } else {
            for (int f = 0; f < kNumFields; ++f) {
              float next_val = (2.0f * mass_matrix(I) * current_field[f](I) -
                                (mass_matrix(I) - 0.5f * dt_local * damping_matrix[f](I)) * prev_field[f](I) -
                                dt2_local * work_vector[f](I));

              if (has_attenuation) {
                for (int l = 0; l < n_sls; ++l) {
                  float const w = sls_w[l];
                  float const gamma = (2.0f - w * dt_local) / (2.0f + w * dt_local);
                  float const beta = sls_beta[l] * w * 2.0f * dt_local / (2.0f + w * dt_local);
                  float const gamma_p = 0.5f + 0.5f * gamma;
                  float const beta_p = 0.5f * beta;

                  next_val += dt2_local * (gamma_p * atten_mem_vars[f](I, l) + beta_p * atten_work_vec[f](I));

                  atten_mem_vars[f](I, l) = gamma * atten_mem_vars[f](I, l) + beta * atten_work_vec[f](I);
                }
              }

              prevprev_field[f](I) = next_val / (mass_matrix(I) + 0.5f * dt_local * damping_matrix[f](I));
              prevprev_field[f](I) *= taper_coeff(I);
              current_field[f](I) *= taper_coeff(I);
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrix() {
  auto mesh_local = m_mesh;
  auto local_massMatrixGlobal = massMatrixGlobal_;
  bool const mask_enabled = m_mask_enabled_;
  auto element_mask = m_element_mask_;
  int const mask_active_value = m_mask_active_value_;

  Kokkos::parallel_for(
      "Solver Compute GMatrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (mask_enabled && element_mask[elementNumber] != mask_active_value) return;
        float massMatrixLocal[kPointsPerElement] = {0};
        int const dim = mesh_local.getOrder() + 1;

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        INTEGRAL_TYPE::computeMassTerm(cornerCoords, [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

        // Acoustic: 1 / (vp^2 rho). Elastic: rho.
        real_t model_factor = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            model_factor =
                1.0f / (mesh_local.getModelVpOnElement(elementNumber) * mesh_local.getModelVpOnElement(elementNumber) *
                        mesh_local.getModelRhoOnElement(elementNumber));
          } else {
            model_factor = mesh_local.getModelRhoOnElement(elementNumber);
          }
        }

        for (int i = 0; i < mesh_local.getNumberOfPointsPerElement(); ++i) {
          int x = i % dim;
          int z = (i / dim) % dim;
          int y = i / (dim * dim);
          int const gIndex = mesh_local.globalNodeIndex(elementNumber, x, y, z);

          if constexpr (IS_MODEL_ON_NODES) {
            if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
              model_factor = 1.0f / (mesh_local.getModelVpOnNodes(gIndex) * mesh_local.getModelVpOnNodes(gIndex) *
                                     mesh_local.getModelRhoOnNodes(gIndex));
            } else {
              model_factor = mesh_local.getModelRhoOnNodes(gIndex);
            }
          }

          massMatrixLocal[i] *= model_factor;
          ATOMICADD(local_massMatrixGlobal[gIndex], massMatrixLocal[i]);
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrix() {
  auto mesh_local = m_mesh;
  bool const mask_enabled = m_mask_enabled_;
  auto element_mask = m_element_mask_;
  int const mask_active_value = m_mask_active_value_;

  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> local_dampingMatrixGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_dampingMatrixGlobal[f] = dampingMatrixGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Compute Damping Matrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        (void)local_dampingMatrixGlobal;
        if (mask_enabled && element_mask[elementNumber] != mask_active_value) return;
        for (int i = 0; i < 6; ++i) {
          int f = mesh_local.getGlobalFace(elementNumber, static_cast<model::CubicFace>(i));

          // Only boundary faces carry an absorbing term.
          if (!mesh_local.isBoundaryFace(f)) continue;

          float coords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) {
              coords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
            }
          }

          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            real_t model_rho = 0.0f;
            real_t model_vp = 0.0f;
            real_t alpha = 0.0f;

            if constexpr (!IS_MODEL_ON_NODES) {
              model_rho = mesh_local.getModelRhoOnElement(elementNumber);
              model_vp = mesh_local.getModelVpOnElement(elementNumber);
              alpha = 1.0 / (model_rho * model_vp);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Free-surface nodes are not damped.
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                model_rho = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                model_vp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                alpha = 1.0 / (model_rho * model_vp);
              }

              real_t localIncrement = alpha * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrement);
            }
          } else {
            float normal[3];
            mesh_local.faceNormal(elementNumber, static_cast<model::CubicFace>(i), normal);
            real_t nx = normal[0], ny = normal[1], nz = normal[2];

            real_t density, velocityVp, velocityVs;

            if constexpr (!IS_MODEL_ON_NODES) {
              density = mesh_local.getModelRhoOnElement(elementNumber);
              velocityVp = mesh_local.getModelVpOnElement(elementNumber);
              velocityVs = mesh_local.getModelVsOnElement(elementNumber);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Free-surface nodes are not damped.
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                density = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                velocityVp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                velocityVs = mesh_local.getModelVsOnNodes(globalNodeIndex);
              }

              real_t aux = density * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              real_t localIncrementx = aux * (velocityVp * fabs(nx) + velocityVs * sqrt(ny * ny + nz * nz));
              real_t localIncrementy = aux * (velocityVp * fabs(ny) + velocityVs * sqrt(nx * nx + nz * nz));
              real_t localIncrementz = aux * (velocityVp * fabs(nz) + velocityVs * sqrt(nx * nx + ny * ny));

              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrementx);
              ATOMICADD(local_dampingMatrixGlobal[1][globalNodeIndex], localIncrementy);
              ATOMICADD(local_dampingMatrixGlobal[2][globalNodeIndex], localIncrementz);
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  massMatrixGlobal_ = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), "massMatrixGlobal");

  static constexpr const char* dampingNames[3] = {"dampingX", "dampingY", "dampingZ"};
  for (int f = 0; f < kNumFields; ++f) {
    dampingMatrixGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), dampingNames[f]);
  }

  static constexpr const char* workVectorNames[3] = {"workVec0", "workVec1", "workVec2"};
  for (int f = 0; f < kNumFields; ++f) {
    workVectorsGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), workVectorNames[f]);
  }

  if (attenuationEnabled_ && nSls_ > 0) {
    static constexpr const char* attWorkNames[3] = {"attWorkVec0", "attWorkVec1", "attWorkVec2"};
    static constexpr const char* attMemNames[3] = {"attMemory0", "attMemory1", "attMemory2"};
    for (int f = 0; f < kNumFields; ++f) {
      attenuationWorkVectorsGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), attWorkNames[f]);
      attenuationMemoryVariables_[f] = allocateArray2D<arrayReal>(m_mesh.getNumberOfNodes(), nSls_, attMemNames[f]);
    }
  }

  spongeTaperCoeff_ = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), "spongeTaperCoeff");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initFEarrays() {
  initSpongeValues();

  if (attenuationEnabled_ && nSls_ > 0) {
    for (int n = 0; n < m_mesh.getNumberOfNodes(); ++n) {
      for (int f = 0; f < kNumFields; ++f) {
        attenuationWorkVectorsGlobal_[f](n) = 0.0f;
        for (int l = 0; l < nSls_; ++l) {
          attenuationMemoryVariables_[f](n, l) = 0.0f;
        }
      }
    }
    FENCE
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initSpongeValues() {
  // @todo VERIFY: unit and origin of the peak damping factor sigma_max = 0.15 (dimensionless?).
  const double sigma_max = 0.15;

  for (int n = 0; n < m_mesh.getNumberOfNodes(); n++) {
    const double x = m_mesh.nodeCoord(n, 0);
    const double y = m_mesh.nodeCoord(n, 1);
    const double z = m_mesh.nodeCoord(n, 2);

    const double distToFrontierX = (surface_sponge_) ? m_mesh.domainSize(0) - x : min(m_mesh.domainSize(0) - x, x);
    const double distToFrontierY = min(m_mesh.domainSize(1) - y, y);
    const double distToFrontierZ = min(m_mesh.domainSize(2) - z, z);

    double minDistToFrontier = max(m_mesh.domainSize(0), max(m_mesh.domainSize(1), m_mesh.domainSize(2)));

    bool is_sponge = false;
    if (distToFrontierX < sponge_size_[0]) {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierX);
    }
    if (distToFrontierY < sponge_size_[1]) {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierY);
    }
    if (distToFrontierZ < sponge_size_[2]) {
      is_sponge = true;
      minDistToFrontier = min(minDistToFrontier, distToFrontierZ);
    }

    // Gaussian taper in the distance to the closest sponge frontier.
    if (is_sponge) {
      double d = minDistToFrontier;
      double delta = taper_delta_;
      double sigma = sigma_max * std::exp(-((d / delta) * (d / delta)));
      spongeTaperCoeff_(n) = 1.0 / (1.0 + sigma);
    } else {
      spongeTaperCoeff_(n) = 1.0;
    }
  }

  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const vectorReal& fieldGlobal, const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << fieldGlobal(m_mesh.globalNodeIndex(e, 0, 0, 0)) << endl;
}

// Builds the TTI tensor CTTI = M * CVTI * M^T, with CVTI the VTI tensor (symmetry axis z) in Voigt storage
// and M the rotation of the symmetry axis given by theta and phi.
// theta and phi are converted from degrees here.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
template <physicType P, typename>
PROXY_HOST_DEVICE void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeCMatrix(
    float const vp, float const vs, float const rho, float const delta, float const epsilon, float const gamma,
    float const phi, float const theta, float (&CTTI)[6][6]) {
  const float rho_vp2 = rho * vp * vp;
  const float rho_vs2 = rho * vs * vs;
  const float two_eps = 2.0f * epsilon;
  const float two_gam = 2.0f * gamma;

  float CVTI[6][6] = {0.0f};
  CVTI[0][0] = rho_vp2 * (1.0f + two_eps);
  CVTI[1][1] = CVTI[0][0];
  CVTI[2][2] = rho_vp2;
  CVTI[3][3] = rho_vs2;
  CVTI[4][4] = CVTI[3][3];
  CVTI[5][5] = rho_vs2 * (1.0f + two_gam);
  CVTI[0][1] = CVTI[0][0] - 2.0f * CVTI[5][5];
  CVTI[1][0] = CVTI[0][1];

  const float vp2_vs2 = vp * vp - vs * vs;
  const float sqrt_arg = vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2;
  CVTI[0][2] = rho * sqrtf(sqrt_arg) - rho_vs2;
  CVTI[1][2] = CVTI[0][2];
  CVTI[2][0] = CVTI[0][2];
  CVTI[2][1] = CVTI[0][2];

  constexpr float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  const float theta_rad = theta * DEG_TO_RAD;
  const float phi_rad = phi * DEG_TO_RAD;

  const float ctheta = cosf(theta_rad);
  const float stheta = sinf(theta_rad);
  const float cphi = cosf(phi_rad);
  const float sphi = sinf(phi_rad);

  const float ct_cp = ctheta * cphi;
  const float ct_sp = ctheta * sphi;
  const float st_cp = stheta * cphi;
  const float st_sp = stheta * sphi;

  float R[3][3];
  R[0][0] = ct_cp;
  R[0][1] = ct_sp;
  R[0][2] = -stheta;
  R[1][0] = -sphi;
  R[1][1] = cphi;
  R[1][2] = 0.0f;
  R[2][0] = st_cp;
  R[2][1] = st_sp;
  R[2][2] = ctheta;

  const float R00_2 = R[0][0] * R[0][0];
  const float R01_2 = R[0][1] * R[0][1];
  const float R02_2 = R[0][2] * R[0][2];
  const float R10_2 = R[1][0] * R[1][0];
  const float R11_2 = R[1][1] * R[1][1];
  const float R12_2 = R[1][2] * R[1][2];
  const float R20_2 = R[2][0] * R[2][0];
  const float R21_2 = R[2][1] * R[2][1];
  const float R22_2 = R[2][2] * R[2][2];

  float M[6][6] = {0.0f};
  M[0][0] = R00_2;
  M[0][1] = R01_2;
  M[0][2] = R02_2;
  M[1][0] = R10_2;
  M[1][1] = R11_2;
  M[1][2] = R12_2;
  M[2][0] = R20_2;
  M[2][1] = R21_2;
  M[2][2] = R22_2;
  // M * CVTI * M^T needs the Bond matrix for stress, whose factors of two sit in
  // the upper right block; the lower left placement builds the one for strain.
  M[0][3] = 2.0f * R[0][1] * R[0][2];
  M[0][4] = 2.0f * R[0][0] * R[0][2];
  M[0][5] = 2.0f * R[0][0] * R[0][1];
  M[1][3] = 2.0f * R[1][1] * R[1][2];
  M[1][4] = 2.0f * R[1][0] * R[1][2];
  M[1][5] = 2.0f * R[1][0] * R[1][1];
  M[2][3] = 2.0f * R[2][1] * R[2][2];
  M[2][4] = 2.0f * R[2][0] * R[2][2];
  M[2][5] = 2.0f * R[2][0] * R[2][1];
  M[3][0] = R[1][0] * R[2][0];
  M[3][1] = R[1][1] * R[2][1];
  M[3][2] = R[1][2] * R[2][2];
  M[3][3] = R[1][1] * R[2][2] + R[1][2] * R[2][1];
  M[3][4] = R[1][0] * R[2][2] + R[1][2] * R[2][0];
  M[3][5] = R[1][0] * R[2][1] + R[1][1] * R[2][0];
  M[4][0] = R[0][0] * R[2][0];
  M[4][1] = R[0][1] * R[2][1];
  M[4][2] = R[0][2] * R[2][2];
  M[4][3] = R[0][1] * R[2][2] + R[0][2] * R[2][1];
  M[4][4] = R[0][0] * R[2][2] + R[0][2] * R[2][0];
  M[4][5] = R[0][0] * R[2][1] + R[0][1] * R[2][0];
  M[5][0] = R[0][0] * R[1][0];
  M[5][1] = R[0][1] * R[1][1];
  M[5][2] = R[0][2] * R[1][2];
  M[5][3] = R[0][1] * R[1][2] + R[0][2] * R[1][1];
  M[5][4] = R[0][0] * R[1][2] + R[0][2] * R[1][0];
  M[5][5] = R[0][0] * R[1][1] + R[0][1] * R[1][0];

  float temp[6][6];

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++) {
        sum += M[i][k] * CVTI[k][j];
      }
      temp[i][j] = sum;
    }
  }

  for (int i = 0; i < 6; i++) {
    for (int j = i; j < 6; j++) {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++) {
        sum += temp[i][k] * M[j][k];
      }
      CTTI[i][j] = sum;
      if (i != j) CTTI[j][i] = sum;
    }
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
template <physicType P, typename>
PROXY_HOST_DEVICE void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeTtiCompact(
    float const vp, float const vs, float const rho, float const delta, float const epsilon, float const gamma,
    [[maybe_unused]] float const phi, float const theta, float (&p)[flux::kTtiCompactSize]) {
  // Same VTI coefficients as computeCMatrix().
  const float rho_vp2 = rho * vp * vp;
  const float rho_vs2 = rho * vs * vs;
  const float c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
  const float c66 = rho_vs2 * (1.0f + 2.0f * gamma);
  const float vp2_vs2 = vp * vp - vs * vs;
  const float c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;

  constexpr float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  const float theta_rad = theta * DEG_TO_RAD;
  float const n[3] = {-sinf(theta_rad), 0.0f, cosf(theta_rad)};
  flux::ttiCompactFromVti(c11, c13, rho_vp2, rho_vs2, c66, n, p);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrixMasked(
    const vectorInt& elem_mask, int active_value) {
  Kokkos::deep_copy(massMatrixGlobal_, 0.0f);
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeGlobalMassMatrix();
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrixMasked(
    const vectorInt& elem_mask, int active_value) {
  for (int f = 0; f < kNumFields; ++f) Kokkos::deep_copy(dampingMatrixGlobal_[f], 0.0f);
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeDampingMatrix();
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsMasked(
    const DataType& data, const vectorInt& elem_mask, int active_value) {
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeElementContributions(data);
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsFromList(
    const DataType& data, const vectorInt& elem_list, int n_elems) {
  m_elem_list_ = elem_list;
  m_n_elem_list_ = n_elems;
  m_list_mode_ = true;
  computeElementContributions(data);
  m_list_mode_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForward(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFieldsForward(dt, data);
  m_node_list_mode_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForwardAndReset(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_reset_forces_in_update_ = true;
  updateFieldsFromListForward(dt, data, node_list, n_nodes);
  m_reset_forces_in_update_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListBackward(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFieldsBackward(dt, data);
  m_node_list_mode_ = false;
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
