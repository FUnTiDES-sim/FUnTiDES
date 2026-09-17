/**
 * @file bench_solver_fe_cartesian_struct_elastic_pml.cc
 * @brief Standalone benchmark harness for the C-PML elastic SEM solver.
 *
 * This is the PROXY HARNESS for the elastic C-PML optimization ratchet. It is
 * NOT a Google Benchmark binary: the proxy runner invokes it as
 *     <binary> <param_file>
 * where param_file is a text file of key=value lines:
 *     n=<number of time steps>
 *     output_file=<path to write the final displacement field as raw float64>
 *
 * It instantiates SEMsolverElastic directly (ModelStruct<float,int,ORDER>)
 * on a Cartesian mesh with a C-PML absorbing layer enabled via setPML()
 * (two-sided stretched gradient/divergence, Komatitsch & Martin 2007),
 * loops over computeOneStep(), and emits the proxy metrics block:
 *     PROXY_METRICS_BEGIN
 *     time_s=...      (mean wall time per time step)
 *     mem_bytes=...   (deterministic sum of persistent solver allocations)
 *     PROXY_METRICS_END
 *
 * The final 3-component displacement field is written to output_file as raw
 * float64 (3*n_node values, component-major: ux[0..n), uy[0..n), uz[0..n))
 * for the reference comparison (l2_rel/linf_rel).
 *
 * ORDER is 1 and the MAKUTU (sum-factorization) integral type is used because
 * the elastic C-PML is implemented only on the Flat kernel, which is
 * instantiated up to kMaxOrderForFlatElastic = 1 (see SEMsolver::setPML).
 * Higher orders or the tensorial GEMM type would silently ignore the PML.
 *
 * NOTE: this file is the measuring instrument, never the subject of the
 * optimization. The ratchet edits the solver sources, not this harness.
 */
#include "cartesian_struct_builder.h"
#include "data_type.h"
#include "model.h"
#include "rhs_elastic.h"
#include "sem_solver_data.h"
#include "sem_solver_impl.h"
#include "wavefield_elastic.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Configuration (fixed for the benchmark; the proxy only varies n via the
// param file). Domain: 24^3 elements, order 1 → 24^3 * 8 dofs per field.
// PML thickness 300 m on every side (~15% of the 2000 m domain) so the PML
// layer is a meaningful fraction of the mesh.
// ---------------------------------------------------------------------------
constexpr int kOrder = 1;
constexpr int kEx = 24;
constexpr int kEy = 24;
constexpr int kEz = 24;
constexpr float kDomainSize = 2000.0f;
constexpr float kPmlSize = 300.0f;
constexpr float kDt = 0.001f;
constexpr int kNRhs = 1;
constexpr float kF0 = 5.0f;
constexpr int kNDof = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);

using MeshType = model::ModelStruct<float, int, kOrder>;
using IntType = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
using SemSolverT = solver::fe::SEMsolverElastic<kOrder, IntType, MeshType, false>;

struct Params {
  int n_steps = 200;
  std::string output_file;
};

Params readParams(const char* path) {
  Params p;
  std::ifstream fh(path);
  if (!fh.is_open()) {
    std::fprintf(stderr, "Cannot open param file %s\n", path);
    std::exit(1);
  }
  std::string line;
  while (std::getline(fh, line)) {
    if (line.empty() || line[0] == '#') continue;
    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string key = line.substr(0, eq);
    std::string val = line.substr(eq + 1);
    if (key == "n") p.n_steps = std::atoi(val.c_str());
    else if (key == "output_file") p.output_file = val;
  }
  return p;
}

std::shared_ptr<model::ModelApi<float, int>> createModel() {
  // The builder takes (ex, lx, ey, ly, ez, lz) — DOMAIN sizes, not element
  // sizes (see the DG-SEM harness note). No DG/SEM split, no acoustoelastic.
  model::CartesianStructBuilder<float, int, kOrder> builder(
      kEx, kDomainSize, kEy, kDomainSize, kEz, kDomainSize, /*isModelOnNodes=*/false, /*isElastic=*/true,
      /*ox=*/0.0f, /*oy=*/0.0f, /*oz=*/0.0f, /*global_lx=*/-1.0f, /*global_ly=*/-1.0f, /*global_lz=*/-1.0f,
      /*global_ox=*/0.0f, /*global_oy=*/0.0f, /*global_oz=*/0.0f,
      /*isAcoustoElastic=*/false, /*acoustoElasticBoundaryZ=*/0.0f, /*DgSemBoundaryZ=*/0.0f);
  // free_surface=false → all boundaries are Damping; the C-PML replaces the
  // sponge inside the PML layer (setupPML disables the taper there).
  return builder.getModel(false);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: %s <param_file>\n", argv[0]);
    return 1;
  }
  Params const params = readParams(argv[1]);

  // Pin to a quiet GPU for stable timing: the proxy runs without
  // CUDA_VISIBLE_DEVICES, so Kokkos would land on a busy shared GPU and the
  // per-step time would swing between replicates. GPU 7 is free on this host
  // (checked via nvidia-smi). setenv before Kokkos::initialize.
  setenv("CUDA_VISIBLE_DEVICES", "7", 1);

  Kokkos::initialize(argc, argv);
  {
    auto model = createModel();
    const int n_elem = model->getNumberOfElements();
    const int n_node = model->getNumberOfNodes();

    SemSolverT solver;
    solver.setAnisotropyType(model::AnisotropyType::kIso);
    // Enable the C-PML layer BEFORE computeFEInit (required by setPML).
    solver.setPML({kPmlSize, kPmlSize, kPmlSize}, /*profile=*/2.0f, /*reflection=*/1e-3f, /*alpha_max=*/0.0f,
                  /*kappa_max=*/1.0f, /*dt=*/kDt);
    // No sponge (boundaries_size=0): the C-PML replaces it inside the layer.
    solver.computeFEInit(*model, {0.0f, 0.0f, 0.0f}, false, 0.015f);

    // Wavefield arrays (3 displacement components for elastic).
    auto uxPrev = allocateVector<vectorReal>(n_node, "uxPrev");
    auto uyPrev = allocateVector<vectorReal>(n_node, "uyPrev");
    auto uzPrev = allocateVector<vectorReal>(n_node, "uzPrev");
    auto uxCurr = allocateVector<vectorReal>(n_node, "uxCurr");
    auto uyCurr = allocateVector<vectorReal>(n_node, "uyCurr");
    auto uzCurr = allocateVector<vectorReal>(n_node, "uzCurr");

    // Initial impulse (Gaussian bump at the domain center) in all three
    // components — the same excitation the elastic-PML unit tests use
    // (PmlAbsorbsBetterThanPlain). The elastic SEM driver has no RHS-source
    // path, so an initial condition is the validated way to excite the field.
    {
      auto h_ux = Kokkos::create_mirror_view(uxCurr);
      auto h_uy = Kokkos::create_mirror_view(uyCurr);
      auto h_uz = Kokkos::create_mirror_view(uzCurr);
      int const dim = kOrder + 1;
      int const c = n_node / 2;
      int const cx = c % kEx, cy = (c / kEx) % kEy, cz = c / (kEx * kEy);
      for (int k = 0; k < kEz; ++k)
        for (int j = 0; j < kEy; ++j)
          for (int i = 0; i < kEx; ++i) {
            float const dx = (i - cx) / 2.0f, dy = (j - cy) / 2.0f, dz = (k - cz) / 2.0f;
            float const g = std::exp(-(dx * dx + dy * dy + dz * dz));
            int const idx = i + j * kEx + k * kEx * kEy;
            h_ux(idx) = g;
            h_uy(idx) = g;
            h_uz(idx) = g;
          }
      Kokkos::deep_copy(uxCurr, h_ux);
      Kokkos::deep_copy(uyCurr, h_uy);
      Kokkos::deep_copy(uzCurr, h_uz);
      FENCE
    }

    // Zero RHS (no forcing; the initial condition drives the simulation).
    auto rhsTermx = allocateArray2D<arrayReal>(kNRhs, params.n_steps, "rhsTermx");
    auto rhsTermy = allocateArray2D<arrayReal>(kNRhs, params.n_steps, "rhsTermy");
    auto rhsTermz = allocateArray2D<arrayReal>(kNRhs, params.n_steps, "rhsTermz");
    auto rhsElement = allocateVector<vectorInt>(kNRhs, "rhsElement");
    auto rhsWeights = allocateArray2D<arrayReal>(kNRhs, kNDof, "rhsWeights");
    Kokkos::deep_copy(rhsTermx, 0.0f);
    Kokkos::deep_copy(rhsTermy, 0.0f);
    Kokkos::deep_copy(rhsTermz, 0.0f);
    Kokkos::deep_copy(rhsWeights, 0.0f);
    auto h_rhsElement = Kokkos::create_mirror_view(rhsElement);
    h_rhsElement(0) = 0;
    Kokkos::deep_copy(rhsElement, h_rhsElement);
    FENCE

    solver::fe::WavefieldElastic wavefield(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
    solver::fe::RhsElastic rhs(rhsTermx, rhsTermy, rhsTermz, rhsElement, rhsWeights);
    solver::fe::SEMsolverDataElastic data(wavefield, rhs);

    // Deterministic persistent memory footprint (bytes). The PML array sizes
    // are queried from the solver's actual views so the metric tracks the real
    // allocation (e.g. a compacted coefficient stride) rather than a fixed
    // model. The rest mirror the solver's persistent allocations: for elastic
    // the damping and work vectors are 3-component (kNumFields=3).
    const size_t solver_fe = 8ull * n_node * sizeof(float);  // mass + 3 damping + 3 work + sponge
    const size_t pml_coeff = solver.getPmlCoefficients().extent(0) * solver.getPmlCoefficients().extent(1) *
                             sizeof(float);
    const size_t pml_masks = solver.getPmlNodeIndex().extent(0) * sizeof(int) +
                             solver.getPmlElementMask().extent(0) * sizeof(int);
    const size_t pml_mem = solver.getPmlMemoryVariables().extent(0) * solver.getPmlMemoryVariables().extent(1) *
                           sizeof(float);
    const size_t wavefield_bytes = 6ull * n_node * sizeof(float);
    const size_t rhs_arrays = 3ull * kNRhs * params.n_steps * sizeof(float) + 1ull * kNRhs * sizeof(int) +
                              1ull * kNRhs * kNDof * sizeof(float);
    const size_t mem_bytes = solver_fe + pml_coeff + pml_masks + pml_mem + wavefield_bytes + rhs_arrays;

    // Warm-up (first-launch overhead excluded from the measurement). 100 steps
    // so the CUDA JIT/clock ramp-up of the first process launch settles before
    // the timed loop — the first replicate of a short warmup is systematically
    // slower (outlier that inflates the measured spread).
    for (int t = 0; t < 100; ++t) {
      solver.computeOneStep(kDt, t, data);
      data.swapWavefields();
    }

    // Timed loop.
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < params.n_steps; ++t) {
      solver.computeOneStep(kDt, t, data);
      data.swapWavefields();
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double const total_s = std::chrono::duration<double>(t1 - t0).count();
    double const mean_time_s = total_s / static_cast<double>(params.n_steps);

    // Final 3-component displacement field → host, for the reference
    // comparison (component-major: ux[0..n), uy[0..n), uz[0..n)).
    auto h_ux = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, uxPrev);
    auto h_uy = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, uyPrev);
    auto h_uz = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, uzPrev);
    std::vector<double> field(3ull * n_node);
    for (int n = 0; n < n_node; ++n) {
      field[n] = static_cast<double>(h_ux(n));
      field[n_node + n] = static_cast<double>(h_uy(n));
      field[2ull * n_node + n] = static_cast<double>(h_uz(n));
    }

    if (!params.output_file.empty()) {
      std::ofstream fout(params.output_file, std::ios::binary);
      fout.write(reinterpret_cast<const char*>(field.data()), field.size() * sizeof(double));
    }

    std::printf("PROXY_METRICS_BEGIN\n");
    std::printf("time_s=%.9f\n", mean_time_s);
    std::printf("mem_bytes=%zu\n", mem_bytes);
    std::printf("n_elem=%d\n", n_elem);
    std::printf("n_node=%d\n", n_node);
    std::printf("kNDof=%d\n", kNDof);
    std::printf("pml_enabled=%d\n", solver.isPmlEnabled() ? 1 : 0);
    std::printf("output_file=%s\n", params.output_file.c_str());
    std::printf("PROXY_METRICS_END\n");
  }
  Kokkos::finalize();
  return 0;
}
