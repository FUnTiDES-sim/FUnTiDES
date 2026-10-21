#ifndef FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_STEP_H_
#define FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_STEP_H_

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
    // The coefficients are completed on the host (where setSLSAttenuation() put them), then copied to the device.
    if (slsAnelasticityCoefficients_.extent(0) > 0) {
      slsAnelasticityCoefficients_.sync_host();
      auto coefficients = slsAnelasticityCoefficients_.view_host();
      for (int l = 0; l < nSls_; ++l) {
        if (coefficients(l) < 0.0f) {
          coefficients(l) = 2.0f * minQVal / (std::max(1.0001f, minQVal) - 1.0f);
        }
      }
      slsAnelasticityCoefficients_.modify_host();
      slsAnelasticityCoefficients_.sync_device();
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


}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_STEP_H_
