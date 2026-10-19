// =============================================================================
// Elastic SEM element-contribution proxy harness (standalone, no Google
// Benchmark).
//
// Purpose: exercise the production Flat (makutu) elastic element-contribution
// kernel on a topologically Cartesian unstructured mesh and report
//   * the device wall time of the kernel (median over `reps` blocks),
//   * the process peak footprint (VmHWM) and the harness-owned allocation size,
//   * an independent cross-check against the Tensorial-GEMM kernel
//     (relDiffVsFlat), and
//   * a deterministic stride sample of the assembled force field, written as
//     raw little-endian float64 for comparison against a sealed reference.
//
// The kernel under measurement is reached only through the real solver
// (solver::SEMsolver<...>::computeElementContributions), never
// re-implemented here.
//
// Usage:
//   bench_elastic_kernel_proxy --order=1 --on-nodes=0 --aniso=0 --ex=32 \
//       --iters=20 --reps=5 --out=/path/field.bin [--oracle=1] [--sample=2000000]
// =============================================================================

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "Integrals.h"
#include "Qk_Hexahedron_Tensorial.h"
#include "cartesian_unstruct_builder.h"
#include "data_type.h"
#include "model.h"
#include "model_unstruct.h"
#include "rhs_elastic.h"
#include "sem_enums.h"
#include "sem_macros.h"
#include "sem_solver.h"
#include "sem_solver_impl.h"
#include "wavefield_elastic.h"

using namespace solver;
using namespace utils::enums;

namespace {

struct Args {
  int order = 1;
  int on_nodes = 0;
  int aniso = 0;
  int ex = 32;
  int iters = 20;
  int reps = 5;
  int oracle = 1;
  long long sample = 2000000;
  std::string out;
};

// Accepts "--key=value" forms only; unrecognised tokens are ignored (Kokkos
// consumes its own "--kokkos-*" flags from argv).
Args parseArgs(int argc, char** argv) {
  Args a;
  auto value = [](const std::string& s, const char* key, std::string& dst) -> bool {
    const std::string k = std::string("--") + key + "=";
    if (s.rfind(k, 0) != 0) return false;
    dst = s.substr(k.size());
    return true;
  };
  for (int i = 1; i < argc; ++i) {
    const std::string s(argv[i]);
    std::string v;
    if (value(s, "order", v)) a.order = std::atoi(v.c_str());
    else if (value(s, "on-nodes", v)) a.on_nodes = std::atoi(v.c_str());
    else if (value(s, "aniso", v)) a.aniso = std::atoi(v.c_str());
    else if (value(s, "ex", v)) a.ex = std::atoi(v.c_str());
    else if (value(s, "iters", v)) a.iters = std::atoi(v.c_str());
    else if (value(s, "reps", v)) a.reps = std::atoi(v.c_str());
    else if (value(s, "oracle", v)) a.oracle = std::atoi(v.c_str());
    else if (value(s, "sample", v)) a.sample = std::atoll(v.c_str());
    else if (value(s, "out", v)) a.out = v;
  }
  if (a.iters < 1) a.iters = 1;
  if (a.reps < 1) a.reps = 1;
  return a;
}

// Peak resident set size of this process: host pages plus any host-mapped
// device allocation. External to the code under measurement.
long long peakRssBytes() {
  FILE* f = std::fopen("/proc/self/status", "r");
  if (f == nullptr) return -1;
  char line[256];
  long long kb = -1;
  while (std::fgets(line, sizeof(line), f) != nullptr) {
    if (std::strncmp(line, "VmHWM:", 6) == 0) {
      std::sscanf(line + 6, "%lld", &kb);
      break;
    }
  }
  std::fclose(f);
  return kb < 0 ? -1 : kb * 1024;
}

template <typename ViewA, typename ViewB>
void accumulateDiff(const ViewA& ha, const ViewB& hb, long long stride, double& dmax, double& rmax) {
  const long long n = static_cast<long long>(ha.extent(0));
  for (long long i = 0; i < n; i += stride) {
    const double x = static_cast<double>(ha(static_cast<int>(i)));
    const double y = static_cast<double>(hb(static_cast<int>(i)));
    dmax = std::max(dmax, std::abs(x - y));
    rmax = std::max(rmax, std::abs(x));
  }
}

// Max relative difference between the assembled force vectors of two solvers,
// evaluated on the same stride sample the field dump uses.
template <typename SolverA, typename SolverB>
double maxRelDiffForceSampled(SolverA& a, SolverB& b, int n_fields, long long stride) {
  double dmax = 0.0, rmax = 0.0;
  for (int f = 0; f < n_fields; ++f) {
    auto ha = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, a.getForceVector(f));
    auto hb = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, b.getForceVector(f));
    accumulateDiff(ha, hb, stride, dmax, rmax);
  }
  return rmax > 0.0 ? dmax / rmax : dmax;
}

// Deterministic, well-conditioned initial state.
//
// Kokkos::View default-initialises to zero, so a wavefield left untouched makes
// the elastic operator act on zero: K*u == 0, the force vector is identically
// zero, and every correctness signal becomes vacuous -- relDiffVsFlat compares
// 0 against 0, and an l2_rel requirement is satisfied by any kernel, including
// one that computes nothing. The kernel performs the same work either way, so
// this changes the accuracy gate and not the timing.
//
// The pattern is an integer hash of the node index: exactly reproducible (no
// libm, no RNG, no platform dependence), so a sealed reference stays valid.
template <class ViewT>
void fillDeterministic(const char* tag, ViewT view, unsigned int salt) {
  const int n = static_cast<int>(view.extent(0));
  Kokkos::parallel_for(
      tag, Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
        const unsigned int h = (static_cast<unsigned int>(i) + salt) * 2654435761u;
        view(i) = static_cast<float>(static_cast<double>(h % 1000003u) / 1000003.0 - 0.5) * 1.0e-3f;
      });
}

template <int ORDER, bool ON_NODES>
int runCase(const Args& a) {
  using MeshT = model::ModelUnstruct<float, int>;
  using BuilderT = model::CartesianUnstructBuilder<float, int>;
  using ParamsT = model::CartesianParams<float, int>;
  using FlatIntegral = typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type;
  using GemmIntegral = typename IntegralTypeSelector<ORDER, IntegralType::TENSORIAL_GEMM>::type;
  using FlatSolver = solver::SEMsolver<ORDER, FlatIntegral, MeshT, ON_NODES, utils::enums::physicType::kElastic>;
  using GemmSolver = solver::SEMsolver<ORDER, GemmIntegral, MeshT, ON_NODES, utils::enums::physicType::kElastic>;

  constexpr int kNumFields = 3;  // elastic particle velocity
  constexpr int n_rhs = 2;
  constexpr int n_time_steps = 1500;
  const float l = 2000.0f;
  const std::array<float, 3> sponge_size = {200.0f, 200.0f, 200.0f};
  const bool surface_sponge = false;
  const float taper_delta = 100.0f;

  ParamsT params(ORDER, a.ex, a.ex, a.ex, l, l, l, ON_NODES != 0, true);
  BuilderT builder(params);
  std::shared_ptr<model::ModelApi<float, int>> mesh = builder.getModel(true);

  const int n_dof = static_cast<int>(mesh->getNumberOfNodes());
  const int n_elem = static_cast<int>(mesh->getNumberOfElements());
  const int n_points = mesh->getNumberOfPointsPerElement();

  // State arrays: same shapes the production time loop and the repo benchmark
  // use, so the kernel runs on a realistic wavefield.
  auto uxPrev = allocateVector<vectorReal>(n_dof, "uxnGlobalPrev");
  auto uyPrev = allocateVector<vectorReal>(n_dof, "uynGlobalPrev");
  auto uzPrev = allocateVector<vectorReal>(n_dof, "uznGlobalPrev");
  auto uxCurr = allocateVector<vectorReal>(n_dof, "uxnGlobalCurr");
  auto uyCurr = allocateVector<vectorReal>(n_dof, "uynGlobalCurr");
  auto uzCurr = allocateVector<vectorReal>(n_dof, "uznGlobalCurr");
  auto rhsTermx = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTermx");
  auto rhsTermy = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTermy");
  auto rhsTermz = allocateArray2D<arrayReal>(n_rhs, n_time_steps, "rhsTermz");
  auto rhsElement = allocateVector<vectorInt>(n_rhs, "rhsElement");
  auto rhsWeights = allocateArray2D<arrayReal>(n_rhs, n_points, "rhsWeights");
  FENCE

  // Distinct salt per component and buffer: the six arrays stay distinguishable,
  // so a kernel that reads the wrong one cannot pass the oracle by symmetry.
  fillDeterministic("init_ux_prev", uxPrev, 11u);
  fillDeterministic("init_uy_prev", uyPrev, 29u);
  fillDeterministic("init_uz_prev", uzPrev, 47u);
  fillDeterministic("init_ux_curr", uxCurr, 71u);
  fillDeterministic("init_uy_curr", uyCurr, 97u);
  fillDeterministic("init_uz_curr", uzCurr, 131u);
  FENCE

  const long long alloc_bytes = static_cast<long long>(n_dof) * 6 * static_cast<long long>(sizeof(float)) +
                                 static_cast<long long>(n_rhs) * n_time_steps * 3 *
                                    static_cast<long long>(sizeof(float)) +
                                static_cast<long long>(n_rhs) * static_cast<long long>(sizeof(int)) +
                                static_cast<long long>(n_rhs) * n_points * static_cast<long long>(sizeof(float));

  auto wavefield = WavefieldElastic(uxPrev, uyPrev, uzPrev, uxCurr, uyCurr, uzCurr);
  auto rhs = RhsElastic(rhsTermx, rhsTermy, rhsTermz, rhsElement, rhsWeights);
  SEMsolverDataElastic data(wavefield, rhs);

  const auto aniso = static_cast<model::AnisotropyType>(a.aniso);

  // --- Flat (makutu) solver: the kernel under measurement --------------------
  FlatSolver flat;
  flat.setAnisotropyType(aniso);
  if (aniso == model::AnisotropyType::kTTI && !ON_NODES) mesh->initElasticityTensors(aniso);
  flat.computeFEInit(*mesh, sponge_size, surface_sponge, taper_delta);

  // Deterministic sample stride: a pure function of the mesh size, so the
  // sealed reference and every later run sample identical entries.
  const long long stride = std::max<long long>(1, (static_cast<long long>(n_dof) * kNumFields) / a.sample);

  for (int w = 0; w < 3; ++w) flat.computeElementContributions(data);
  Kokkos::fence();

  std::vector<double> times;
  times.reserve(static_cast<size_t>(a.reps));
  for (int r = 0; r < a.reps; ++r) {
    Kokkos::fence();
    Kokkos::Timer timer;
    for (int it = 0; it < a.iters; ++it) flat.computeElementContributions(data);
    Kokkos::fence();
    times.push_back(timer.seconds() / static_cast<double>(a.iters));
  }
  std::sort(times.begin(), times.end());
  const double t_med = times[times.size() / 2];
  const double t_min = times.front();
  const double t_max = times.back();

  // --- single-application state for the oracle and the exported field --------
  // workVectorsGlobal_ is ACCUMULATED by the kernel, so the timed loop above has
  // left `flat` holding 3 + iters*reps = 103 applications, while a freshly
  // constructed GEMM solver holds one. Comparing the two therefore measures the
  // accumulation count and not the kernels: the discrepancy is then exactly
  // 102/103 = 0.990291, which is what this harness reported before this fix --
  // a 99% "disagreement" that was pure bookkeeping. Both sides are zeroed and
  // advanced by exactly one application before any comparison or export.
  auto zeroForces = [](auto& s) {
    for (int f = 0; f < kNumFields; ++f) Kokkos::deep_copy(s.getForceVector(f), 0.0f);
  };

  zeroForces(flat);
  flat.computeElementContributions(data);
  Kokkos::fence();

  // --- independent cross-check: Tensorial GEMM ------------------------------
  double rel_vs_flat = -1.0;
  if (a.oracle != 0) {
    GemmSolver gemm;
    gemm.setAnisotropyType(aniso);
    gemm.computeFEInit(*mesh, sponge_size, surface_sponge, taper_delta);
    zeroForces(gemm);
    gemm.computeElementContributions_Elastic_Gemm(data);
    Kokkos::fence();
    rel_vs_flat = maxRelDiffForceSampled(flat, gemm, kNumFields, stride);
  }

  // --- field dump (raw float64, deterministic stride sample) ----------------
  if (!a.out.empty()) {
    std::vector<double> buf;
    for (int f = 0; f < kNumFields; ++f) {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, flat.getForceVector(f));
      const long long n = static_cast<long long>(h.extent(0));
      for (long long i = 0; i < n; i += stride) buf.push_back(static_cast<double>(h(static_cast<int>(i))));
    }
    FILE* fo = std::fopen(a.out.c_str(), "wb");
    if (fo == nullptr) {
      std::fprintf(stderr, "cannot open output file '%s'\n", a.out.c_str());
      return 3;
    }
    const size_t written = std::fwrite(buf.data(), sizeof(double), buf.size(), fo);
    std::fclose(fo);
    if (written != buf.size()) {
      std::fprintf(stderr, "short write to '%s'\n", a.out.c_str());
      return 3;
    }
  }

  std::printf("PROXY_METRICS_BEGIN\n");
  std::printf("time_s=%.9g\n", t_med);
  std::printf("time_min_s=%.9g\n", t_min);
  std::printf("time_max_s=%.9g\n", t_max);
  std::printf("mem_bytes=%lld\n", peakRssBytes());
  std::printf("alloc_bytes=%lld\n", alloc_bytes);
  std::printf("relDiffVsFlat=%.6g\n", rel_vs_flat);
  std::printf("n_elem=%d\n", n_elem);
  std::printf("n_node=%d\n", n_dof);
  std::printf("n_points_per_element=%d\n", n_points);
  std::printf("order=%d\n", ORDER);
  std::printf("on_nodes=%d\n", ON_NODES ? 1 : 0);
  std::printf("aniso=%d\n", a.aniso);
  std::printf("ex=%d\n", a.ex);
  std::printf("iters=%d\n", a.iters);
  std::printf("reps=%d\n", a.reps);
  std::printf("sample_stride=%lld\n", stride);
  std::printf("backend=%s\n", Kokkos::DefaultExecutionSpace::name());
  std::printf("output_file=%s\n", a.out.empty() ? "-" : a.out.c_str());
  std::printf("PROXY_METRICS_END\n");
  std::fflush(stdout);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  Args a = parseArgs(argc, argv);
  Kokkos::initialize(argc, argv);
  int rc = 0;
  {
    switch (a.order * 10 + (a.on_nodes != 0 ? 1 : 0)) {
      case 10: rc = runCase<1, false>(a); break;
      case 11: rc = runCase<1, true>(a); break;
      case 20: rc = runCase<2, false>(a); break;
      case 21: rc = runCase<2, true>(a); break;
      case 40: rc = runCase<4, false>(a); break;
      case 41: rc = runCase<4, true>(a); break;
      default:
        std::fprintf(stderr, "unsupported --order/--on-nodes combination\n");
        rc = 2;
    }
  }
  Kokkos::finalize();
  return rc;
}
