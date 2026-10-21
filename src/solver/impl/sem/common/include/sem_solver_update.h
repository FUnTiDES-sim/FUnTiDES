#ifndef FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_UPDATE_H_
#define FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_UPDATE_H_

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
  auto taper_coeff = spongeTaperCoeff_.view_device();

  auto sls_w = slsReferenceAngularFrequencies_.view_device();
  auto sls_beta = slsAnelasticityCoefficients_.view_device();

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
  // Every per-node value is loaded before the first test and the first store.
  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
    int const n_iter = list_on ? m_n_node_list_ : mesh_local.getNumberOfNodes();
    // Two nodes per thread, n_threads apart so that the loads of a warp stay contiguous.
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
  auto taper_coeff = spongeTaperCoeff_.view_device();

  auto sls_w = slsReferenceAngularFrequencies_.view_device();
  auto sls_beta = slsAnelasticityCoefficients_.view_device();

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


}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_UPDATE_H_
