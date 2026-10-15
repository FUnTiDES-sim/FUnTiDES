#ifndef TESTS_BENCHMARKS_CPP_INCLUDE_PROXY_SEM_COMMON_H_
#define TESTS_BENCHMARKS_CPP_INCLUDE_PROXY_SEM_COMMON_H_

/**
 * @file proxy_sem_common.h
 * @brief Shared boilerplate of the SEM proxy harnesses
 * (proxy_sem_acoustic.cc, proxy_sem_elastic.cc).
 *
 * These harnesses are driven by the MIMIR proxy server: they run one fixed,
 * fully deterministic problem and (1) print a metrics block with one
 * key=value per line (time_s, finite, max_abs, n_values), (2) write the final
 * wavefield to the file given on the command line as raw float64 values,
 * component after component, so the server can seal it as a reference and
 * compute l2_rel of later runs against it.
 *
 * Determinism rules: Kokkos views are NOT zero-initialised, so every array is
 * explicitly filled here; the only excitation is a single-node initial pulse
 * (source weights are zero), so no dependence on source-term indexing.
 */

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>

#include "cartesian_struct_builder.h"
#include "data_type.h"
#include "model.h"
#include "sem_enums.h"
#include "solver.h"
#include "solver_factory.h"

namespace solver {
namespace fe {
namespace proxy {

/// Fixed problem parameters, shared by both harnesses.
constexpr int kEx = 24;
constexpr int kEy = 24;
constexpr int kEz = 24;
constexpr float kLx = 2000.0f;
constexpr float kLy = 2000.0f;
constexpr float kLz = 2000.0f;
constexpr int kSteps = 100;
constexpr float kDt = 0.001f;

/**
 * @brief Builds the structured Cartesian mesh of the proxy problem.
 *
 * The builder signature takes DOMAIN sizes (lx, ly, lz), not element sizes.
 * getModel(false) selects all-damping boundaries (no free surface).
 */
template <int Order>
std::shared_ptr<model::ModelApi<float, int>> buildStructModel(bool is_elastic) {
  model::CartesianStructBuilder<float, int, Order> builder(kEx, kLx, kEy, kLy, kEz, kLz, false, is_elastic);
  return builder.getModel(false);
}

/// @copydoc buildStructModel(bool) with a runtime order (1 or 2).
inline std::shared_ptr<model::ModelApi<float, int>> buildStructModel(int order, bool is_elastic) {
  switch (order) {
    case 1:
      return buildStructModel<1>(is_elastic);
    case 2:
      return buildStructModel<2>(is_elastic);
    default:
      throw std::runtime_error("proxy harness supports orders 1 and 2 only");
  }
}

/**
 * @brief Creates a structured Makutu SEM solver through the public factory,
 * the same entry point every production caller uses.
 */
inline std::unique_ptr<Solver> makeSemSolver(utils::enums::physicType physic, int order) {
  auto solver = solver_factory::createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu,
                                             utils::enums::meshType::kStruct,
                                             utils::enums::modelLocationType::kOnElements, physic, order);
  solver->setAnisotropyType(model::AnisotropyType::kIso);
  return solver;
}

/// @brief Writes one integer into a rank-1 int view (device-space safe).
inline void setInt(vectorInt& v, int i, int value) {
  auto mirror = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(v, mirror);  // mirror content does not matter: value is set before copy-back.
  mirror(i) = value;
  Kokkos::deep_copy(v, mirror);
}

/**
 * @brief Copies the fields to host mirrors, scans them for NaN/Inf and the
 * maximum absolute value, and writes them as raw float64 to @p path.
 *
 * @return Number of values written (0 if the file could not be opened).
 */
inline long writeAndScanFields(const std::string& path, const std::vector<vectorReal>& fields, bool& finite,
                               double& max_abs) {
  finite = true;
  max_abs = 0.0;
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    return 0;
  }
  long count = 0;
  for (const auto& field : fields) {
    auto mirror = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(mirror, field);
    for (size_t i = 0; i < mirror.extent(0); ++i) {
      const double value = static_cast<double>(mirror(i));
      if (!std::isfinite(value)) {
        finite = false;
      }
      max_abs = std::max(max_abs, std::abs(value));
      out.write(reinterpret_cast<const char*>(&value), sizeof(double));
      ++count;
    }
  }
  out.close();
  return count;
}

/// @brief Prints the metrics block parsed by the proxy server.
inline void printMetrics(double time_s_per_step, bool finite, double max_abs, long n_values) {
  std::printf("time_s=%.9f\n", time_s_per_step);
  std::printf("finite=%d\n", finite ? 1 : 0);
  std::printf("max_abs=%.6e\n", max_abs);
  std::printf("n_values=%ld\n", n_values);
}

}  // namespace proxy
}  // namespace fe
}  // namespace solver

#endif  // TESTS_BENCHMARKS_CPP_INCLUDE_PROXY_SEM_COMMON_H_
