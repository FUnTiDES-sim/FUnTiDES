/**
 * @file proxy_sem_elastic_tti.cc
 * @brief Proxy harness: single-physics elastic SEM solver with TTI
 *        anisotropy, fixed problem.
 *
 * Same fixed problem as proxy_sem_elastic.cc (24^3, order 1, single-node ux
 * pulse, 100 steps) but the model carries tilted transverse isotropy:
 * initElasticityTensors(kTTI) before computeFEInit and
 * setAnisotropyType(kTTI) on the solver — the exact recipe of the green
 * gtest SemSolverElasticTtiTest. This exercises the TTI path
 * (computeElementContributions_Tti, computeTtiCompact) that the isotropic
 * oracle cannot see: the 2026-10-15 TTI port showed a fix in these kernels
 * is invisible to the existing sealed references.
 *
 * Driven by the MIMIR proxy server (see proxy_sem_common.h). Usage:
 *   proxy_sem_elastic_tti [output-file]
 * Prints one metrics block (time_s, finite, max_abs, n_values) and writes
 * the final (ux, uy, uz) fields as raw float64, component after component,
 * to the output file.
 */

#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include "proxy_sem_common.h"
#include "rhs_elastic.h"
#include "sem_solver_data.h"
#include "wavefield_elastic.h"

namespace {

using namespace solver::proxy;
using solver::SEMsolverDataElastic;
using solver::RhsElastic;
using solver::WavefieldElastic;

// Order 1 mirrors proxy_sem_elastic.cc so the TTI and Iso oracles differ
// only in the constitutive law: 24^3 elements -> 25^3 = 15,625 nodes x 3 dofs.
constexpr int kOrder = 1;

}  // namespace

int main(int argc, char** argv) {
  const std::string output = argc > 1 ? argv[1] : "proxy_field_elastic_tti.bin";

  Kokkos::ScopeGuard guard(argc, argv);

  auto mesh = buildStructModel<kOrder>(true);
  // TTI tensors must be initialised before computeFEInit reads them
  // (recipe of SemSolverElasticTtiTest::SetUp).
  mesh->initElasticityTensors(model::AnisotropyType::kTTI);
  auto solver = makeSemSolver(utils::enums::physicType::kElastic, kOrder, model::AnisotropyType::kTTI);

  const int numNodes = mesh->getNumberOfNodes();
  const int npp = mesh->getNumberOfPointsPerElement();

  solver->computeFEInit(*mesh, {0.0f, 0.0f, 0.0f}, false, 0.0f);

  // Wavefield: two buffers per component, everything explicitly initialised.
  vectorReal uxPrev = allocateVector<vectorReal>(numNodes, "uxPrev");
  vectorReal uyPrev = allocateVector<vectorReal>(numNodes, "uyPrev");
  vectorReal uzPrev = allocateVector<vectorReal>(numNodes, "uzPrev");
  vectorReal uxCurr = allocateVector<vectorReal>(numNodes, "uxCurr");
  vectorReal uyCurr = allocateVector<vectorReal>(numNodes, "uyCurr");
  vectorReal uzCurr = allocateVector<vectorReal>(numNodes, "uzCurr");
  Kokkos::deep_copy(uxPrev, 0.0f);
  Kokkos::deep_copy(uyPrev, 0.0f);
  Kokkos::deep_copy(uzPrev, 0.0f);
  Kokkos::deep_copy(uxCurr, 0.0f);
  Kokkos::deep_copy(uyCurr, 0.0f);
  Kokkos::deep_copy(uzCurr, 0.0f);
  {
    auto mirror = Kokkos::create_mirror_view(uxCurr);
    Kokkos::deep_copy(mirror, uxCurr);
    mirror(numNodes / 2) = 1.0f;  // single-node initial pulse on ux: the only excitation.
    Kokkos::deep_copy(uxCurr, mirror);
  }

  // Source term: one source, zero weights -> no RHS injection.
  arrayReal rhsTermx = allocateArray2D<arrayReal>(1, kSteps, "rhsTermx");
  arrayReal rhsTermy = allocateArray2D<arrayReal>(1, kSteps, "rhsTermy");
  arrayReal rhsTermz = allocateArray2D<arrayReal>(1, kSteps, "rhsTermz");
  vectorInt rhsElem = allocateVector<vectorInt>(1, "rhsElem");
  arrayReal rhsWeights = allocateArray2D<arrayReal>(1, npp, "rhsWeights");
  Kokkos::deep_copy(rhsTermx, 0.0f);
  Kokkos::deep_copy(rhsTermy, 0.0f);
  Kokkos::deep_copy(rhsTermz, 0.0f);
  Kokkos::deep_copy(rhsWeights, 0.0f);
  setInt(rhsElem, 0, 0);

  WavefieldElastic wavefield(uxPrev, uyPrev, uzPrev, uxCurr, uyCurr, uzCurr);
  RhsElastic rhs(rhsTermx, rhsTermy, rhsTermz, rhsElem, rhsWeights);
  SEMsolverDataElastic data(wavefield, rhs);

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
  const long nValues = writeAndScanFields(
      output, {data.getCurrentField(0), data.getCurrentField(1), data.getCurrentField(2)}, finite, maxAbs);

  printMetrics(timePerStep, finite, maxAbs, nValues);
  std::printf("nodes=%d order=%d steps=%d anisotropy=tti\n", numNodes, kOrder, kSteps);
  return 0;
}
