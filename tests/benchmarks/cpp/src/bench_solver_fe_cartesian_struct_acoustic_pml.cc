/**
 * @file bench_solver_fe_cartesian_struct_acoustic_pml.cc
 * @brief Standalone benchmark harness for the C-PML acoustic SEM solver.
 *
 * This is the PROXY HARNESS for the C-PML optimization ratchet. It is NOT a
 * Google Benchmark binary: the proxy runner invokes it as
 *     <binary> <param_file>
 * where param_file is a text file of key=value lines:
 *     n=<number of time steps>
 *     output_file=<path to write the final pressure field as raw float64>
 *
 * It instantiates SEMsolverAcoustic directly (ModelStruct<float,int,ORDER>)
 * on a Cartesian mesh with a C-PML absorbing layer enabled via setPML()
 * (two-sided stretched gradient/divergence, Komatitsch & Martin 2007),
 * loops over computeOneStep(), and emits the proxy metrics block:
 *     PROXY_METRICS_BEGIN
 *     time_s=...      (mean wall time per time step)
 *     mem_bytes=...   (deterministic sum of persistent solver allocations)
 *     PROXY_METRICS_END
 *
 * The final pressure field is written to output_file as raw float64
 * (n_node values, row-major) for the reference comparison (l2_rel/linf_rel).
 *
 * The tensorial GEMM integral type is used so the GEMM-vs-Flat dispatch in
 * computeElementContributions_Acoustic is live: with PML enabled the solver
 * is forced onto the sum-factorization (Flat) kernel, which is the code path
 * this benchmark measures.
 *
 * NOTE: this file is the measuring instrument, never the subject of the
 * optimization. The ratchet edits the solver sources, not this harness.
 */
#include "cartesian_struct_builder.h"
#include "data_type.h"
#include "model.h"
#include "sem_solver_data.h"
#include "sem_solver_impl.h"
#include "source_and_receiver_utils.h"
#include "utils.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Configuration (fixed for the benchmark; the proxy only varies n via the
// param file). Domain: 24^3 elements, order 2 → 24^3 * 27 dofs per field.
// PML thickness 300 m on every side (~15% of the 2000 m domain) so the PML
// layer is a meaningful fraction of the mesh.
// ---------------------------------------------------------------------------
constexpr int kOrder = 2;
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
using IntType = typename IntegralTypeSelector<kOrder, IntegralType::TENSORIAL_GEMM>::type;
using SemSolverT = solver::fe::SEMsolverAcoustic<kOrder, IntType, MeshType, false>;

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
      kEx, kDomainSize, kEy, kDomainSize, kEz, kDomainSize, /*isModelOnNodes=*/false, /*isElastic=*/false,
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
    // Enable the C-PML layer BEFORE computeFEInit (required by setPML).
    solver.setPML({kPmlSize, kPmlSize, kPmlSize}, /*profile=*/2.0f, /*reflection=*/1e-3f, /*alpha_max=*/0.0f,
                  /*kappa_max=*/1.0f, /*dt=*/kDt);
    // No sponge (boundaries_size=0): the C-PML replaces it inside the layer.
    solver.computeFEInit(*model, {0.0f, 0.0f, 0.0f}, false, 0.015f);

    // Wavefield + RHS arrays (mirror the SEM driver in sem_proxy.cc).
    auto pnGlobalPrev = allocateVector<vectorReal>(n_node, "pnGlobalPrev");
    auto pnGlobalCurr = allocateVector<vectorReal>(n_node, "pnGlobalCurr");
    auto rhsTerm = allocateArray2D<arrayReal>(kNRhs, params.n_steps, "rhsTerm");
    auto rhsElement = allocateVector<vectorInt>(kNRhs, "rhsElement");
    auto rhsWeights = allocateArray2D<arrayReal>(kNRhs, kNDof, "rhsWeights");

    // NOTE: the device views are CudaSpace; host code must fill host mirrors
    // and deep_copy, never write to the device views directly.
    auto h_rhsElement = Kokkos::create_mirror_view(rhsElement);
    auto h_rhsTerm = Kokkos::create_mirror_view(rhsTerm);
    auto h_rhsWeights = Kokkos::create_mirror_view(rhsWeights);

    // Source in the interior (mirrors sem_proxy.cc InitSource):
    //   floor((rel_x*ex)/lx) + floor((rel_y*ey)/ly)*ex + floor((rel_z*ez)/lz)*ey*ex
    constexpr float kSrcX = 1000.0f, kSrcY = 1000.0f, kSrcZ = 1000.0f;
    const int src_elem = static_cast<int>(std::floor((kSrcX * kEx) / kDomainSize)) +
                         static_cast<int>(std::floor((kSrcY * kEy) / kDomainSize)) * kEx +
                         static_cast<int>(std::floor((kSrcZ * kEz) / kDomainSize)) * kEy * kEx;
    h_rhsElement(0) = src_elem;

    SolverUtils utils;
    float const tpeak = 1.0f / kF0;
    std::vector<float> sourceTerm = utils.computeSourceTerm(params.n_steps, kDt, kF0, 2, tpeak);
    for (int j = 0; j < params.n_steps; ++j) {
      h_rhsTerm(0, j) = sourceTerm[j];
    }

    // Lagrange interpolation weights at the source point.
    auto* typed_model = dynamic_cast<MeshType*>(model.get());
    if (!typed_model) {
      std::fprintf(stderr, "Model is not ModelStruct\n");
      return 1;
    }
    for (int s = 0; s < kNRhs; ++s) {
      const int e = h_rhsElement(s);
      float cornerCoords[8][3];
      int corner_iter = 0;
      for (int kv = 0; kv < 2; ++kv)
        for (int jv = 0; jv < 2; ++jv)
          for (int iv = 0; iv < 2; ++iv) {
            const int gn = typed_model->globalNodeIndex(e, iv, jv, kv);
            cornerCoords[corner_iter][0] = typed_model->nodeCoord(gn, 0);
            cornerCoords[corner_iter][1] = typed_model->nodeCoord(gn, 1);
            cornerCoords[corner_iter][2] = typed_model->nodeCoord(gn, 2);
            corner_iter++;
          }
      std::array<float, 3> srcCoord{kSrcX, kSrcY, kSrcZ};
      SourceAndReceiverUtils::ComputeRHSWeights<kOrder>(cornerCoords, srcCoord, h_rhsWeights);
    }
    Kokkos::deep_copy(rhsElement, h_rhsElement);
    Kokkos::deep_copy(rhsTerm, h_rhsTerm);
    Kokkos::deep_copy(rhsWeights, h_rhsWeights);
    FENCE

    solver::fe::WavefieldAcoustic wavefield(pnGlobalPrev, pnGlobalCurr);
    solver::fe::RhsAcoustic rhs(rhsTerm, rhsElement, rhsWeights);
    solver::fe::SEMsolverDataAcoustic data(wavefield, rhs);

    // Deterministic persistent memory footprint (bytes). The PML array sizes
    // are queried from the solver's actual views so the metric tracks the real
    // allocation (e.g. a compacted coefficient stride) rather than a fixed
    // model. The rest mirror the solver's persistent allocations.
    const size_t solver_fe = 4ull * n_node * sizeof(float);  // mass + damping + work + sponge
    const size_t pml_coeff = solver.getPmlCoefficients().extent(0) * solver.getPmlCoefficients().extent(1) *
                             sizeof(float);
    const size_t pml_masks = solver.getPmlNodeIndex().extent(0) * sizeof(int) +
                             solver.getPmlElementMask().extent(0) * sizeof(int);
    const size_t pml_mem = solver.getPmlMemoryVariables().extent(0) * solver.getPmlMemoryVariables().extent(1) *
                           sizeof(float);
    const size_t wavefield_bytes = 2ull * n_node * sizeof(float);
    const size_t rhs_arrays = 1ull * kNRhs * params.n_steps * sizeof(float) + 1ull * kNRhs * sizeof(int) +
                              1ull * kNRhs * kNDof * sizeof(float);
    const size_t mem_bytes = solver_fe + pml_coeff + pml_masks + pml_mem + wavefield_bytes + rhs_arrays;

    // Warm-up (first-launch overhead excluded from the measurement).
    for (int t = 0; t < 5; ++t) {
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

    // Final pressure field → host, for the reference comparison.
    auto h_pn = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, pnGlobalPrev);
    std::vector<double> field(n_node);
    for (int n = 0; n < n_node; ++n) field[n] = static_cast<double>(h_pn(n));

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
