/**
 * @file proxy_sem_acoustic.cc
 * @brief Proxy harness: single-physics acoustic SEM solver, fixed problem.
 *
 * Driven by the MIMIR proxy server (see proxy_sem_common.h). Usage:
 *   proxy_sem_acoustic [output-file]
 * Prints one metrics block (time_s, finite, max_abs, n_values) and writes the
 * final pressure field as raw float64 to the output file.
 */

#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include "proxy_sem_common.h"
#include "rhs_acoustic.h"
#include "sem_solver_data.h"
#include "wavefield_acoustic.h"

namespace {

using namespace solver::fe::proxy;
using solver::fe::SEMsolverDataAcoustic;
using solver::fe::RhsAcoustic;
using solver::fe::WavefieldAcoustic;

constexpr int kOrder = 2;  // 24^3 elements, order 2 -> 49^3 = 117,649 dofs.

}  // namespace

int main(int argc, char** argv) {
  const std::string output = argc > 1 ? argv[1] : "proxy_field_acoustic.bin";

  Kokkos::ScopeGuard guard(argc, argv);

  auto mesh = buildStructModel<kOrder>(false);
  auto solver = makeSemSolver(utils::enums::physicType::kAcoustic, kOrder);

  const int numNodes = mesh->getNumberOfNodes();
  const int npp = mesh->getNumberOfPointsPerElement();

  solver->computeFEInit(*mesh, {0.0f, 0.0f, 0.0f}, false, 0.0f);

  // Wavefield: two buffers, everything explicitly initialised.
  vectorReal pPrev = allocateVector<vectorReal>(numNodes, "pPrev");
  vectorReal pCurr = allocateVector<vectorReal>(numNodes, "pCurr");
  Kokkos::deep_copy(pPrev, 0.0f);
  Kokkos::deep_copy(pCurr, 0.0f);
  {
    auto mirror = Kokkos::create_mirror_view(pCurr);
    Kokkos::deep_copy(mirror, pCurr);
    mirror(numNodes / 2) = 1.0f;  // single-node initial pulse: the only excitation.
    Kokkos::deep_copy(pCurr, mirror);
  }

  // Source term: one source, zero weights -> no RHS injection.
  arrayReal rhsTerm = allocateArray2D<arrayReal>(1, kSteps, "rhsTerm");
  vectorInt rhsElem = allocateVector<vectorInt>(1, "rhsElem");
  arrayReal rhsWeights = allocateArray2D<arrayReal>(1, npp, "rhsWeights");
  Kokkos::deep_copy(rhsTerm, 0.0f);
  Kokkos::deep_copy(rhsWeights, 0.0f);
  setInt(rhsElem, 0, 0);

  WavefieldAcoustic wavefield(pPrev, pCurr);
  RhsAcoustic rhs(rhsTerm, rhsElem, rhsWeights);
  SEMsolverDataAcoustic data(wavefield, rhs);

  // Timed loop. The fences make time_s an honest device-side measurement.
  Kokkos::fence();
  const auto t0 = std::chrono::steady_clock::now();
  for (int t = 0; t < kSteps; ++t) {
    solver->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  Kokkos::fence();
  const auto t1 = std::chrono::steady_clock::now();

  const double seconds = std::chrono::duration<double>(t1 - t0).count();
  const double timePerStep = seconds / kSteps;

  bool finite = false;
  double maxAbs = 0.0;
  const long nValues = writeAndScanFields(output, {data.getCurrentField(0)}, finite, maxAbs);

  printMetrics(timePerStep, finite, maxAbs, nValues);
  std::printf("nodes=%d order=%d steps=%d\n", numNodes, kOrder, kSteps);
  return 0;
}
