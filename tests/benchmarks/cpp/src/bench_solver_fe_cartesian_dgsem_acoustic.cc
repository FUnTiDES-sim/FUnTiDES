/**
 * @file bench_solver_fe_cartesian_dgsem_acoustic.cc
 * @brief Standalone benchmark harness for the DG-SEM coupled acoustic solver.
 *
 * This is the PROXY HARNESS for the DG-SEM optimization ratchet. It is NOT a
 * Google Benchmark binary: the proxy runner invokes it as
 *     <binary> <param_file>
 * where param_file is a text file of key=value lines:
 *     n=<number of time steps>
 *     output_file=<path to write the final DG field as raw float64>
 *
 * It instantiates DGSEMsolver directly (ModelStruct<float,int,ORDER>) on a
 * Cartesian mesh split into a DG half (z < lz/2) and an SEM half (z >= lz/2)
 * via setElementTags(), loops over computeOneStep(), and emits the proxy
 * metrics block:
 *     PROXY_METRICS_BEGIN
 *     time_s=...      (mean wall time per time step)
 *     mem_bytes=...   (deterministic sum of persistent solver allocations)
 *     PROXY_METRICS_END
 *
 * The final DG pressure field is written to output_file as raw float64
 * (n_elem x n_dof, row-major) for the reference comparison (l2_rel/linf_rel).
 *
 * NOTE: this file is the measuring instrument, never the subject of the
 * optimization. The ratchet edits the solver sources, not this harness.
 */
#include "cartesian_struct_builder.h"
#include "data_type.h"
#include "dg-sem_solver_data.h"
#include "dg-sem_solver_impl.h"
#include "model.h"
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
// param file). Domain: 24^3 elements, order 2 → 24^3 * 27 dofs per DG field.
// ---------------------------------------------------------------------------
constexpr int kOrder = 2;
constexpr int kEx = 24;
constexpr int kEy = 24;
constexpr int kEz = 24;
constexpr float kDomainSize = 2000.0f;
constexpr float kDt = 0.001f;
constexpr int kNRhs = 1;
constexpr float kF0 = 5.0f;
constexpr int kNDof = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);

using MeshType = model::ModelStruct<float, int, kOrder>;
using IntType = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
using DgSemSolverT = solver::fe::DGSEMsolver<kOrder, IntType, MeshType, false, utils::enums::physicType::kAcoustic>;

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
  // NOTE: the builder takes (ex, lx, ey, ly, ez, lz) — DOMAIN sizes, not element
  // sizes. Passing h = kDomainSize/kEx here made the mesh 83.33 units tall
  // (one element per axis) and broke the DG/SEM split. Use kDomainSize directly.
  model::CartesianStructBuilder<float, int, kOrder> builder(
      kEx, kDomainSize, kEy, kDomainSize, kEz, kDomainSize, /*isModelOnNodes=*/false, /*isElastic=*/false,
      /*ox=*/0.0f, /*oy=*/0.0f, /*oz=*/0.0f, /*global_lx=*/-1.0f, /*global_ly=*/-1.0f, /*global_lz=*/-1.0f,
      /*global_ox=*/0.0f, /*global_oy=*/0.0f, /*global_oz=*/0.0f,
      /*isAcoustoElastic=*/false, /*acoustoElasticBoundaryZ=*/0.0f,
      /*DgSemBoundaryZ=*/kDomainSize / 2.0f);
  // Match the reference driver (sem_proxy.cc): free_surface=false → all
  // boundaries are Damping, exactly like the reference V100 state.
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
  // per-step time would swing 4-5x between replicates. GPU 7 is the free one
  // on this host (checked via nvidia-smi). setenv before Kokkos::initialize.
  setenv("CUDA_VISIBLE_DEVICES", "7", 1);

  Kokkos::initialize(argc, argv);
  {
    auto model = createModel();
    const int n_elem = model->getNumberOfElements();
    const int n_node = model->getNumberOfNodes();

    DgSemSolverT solver;
    // Match the reference driver: no setElementTags (Z-threshold heuristic),
    // no sponge (boundaries_size=0), taper_delta=0.015.
    solver.computeFEInit(*model, {0.0f, 0.0f, 0.0f}, false, 0.015f);

    // Wavefield + RHS arrays (mirror the DG-SEM driver in sem_proxy.cc).
    auto pnDGPrev = allocateArray2D<arrayReal>(n_elem, kNDof, "pnDGPrev");
    auto pnDGCurr = allocateArray2D<arrayReal>(n_elem, kNDof, "pnDGCurr");
    auto pnSEMPrev = allocateVector<vectorReal>(n_node, "pnSEMPrev");
    auto pnSEMCurr = allocateVector<vectorReal>(n_node, "pnSEMCurr");
    auto rhsTermDG = allocateArray2D<arrayReal>(kNRhs, params.n_steps, "rhsTermDG");
    auto rhsTermSEM = allocateArray2D<arrayReal>(kNRhs, params.n_steps, "rhsTermSEM");
    auto rhsElement = allocateVector<vectorInt>(kNRhs, "rhsElement");
    auto rhsWeights = allocateArray2D<arrayReal>(kNRhs, kNDof, "rhsWeights");

    // Source in the SEM half (mirrors sem_proxy.cc: source goes to the domain
    // containing src_coord; the reference driver places it in SEM). The DG
    // source term stays zero so only one domain is injected. Weights are the
    // Lagrange interpolation weights at the source point (ComputeRHSWeights),
    // not uniform 1/kNDof — uniform weights excite non-physical modes and
    // diverge at the CFL-realistic dt.
    //
    // NOTE: the device views are CudaSpace; host code must fill host mirrors
    // and deep_copy, never write to the device views directly.
    auto h_rhsElement = Kokkos::create_mirror_view(rhsElement);
    auto h_rhsTermDG = Kokkos::create_mirror_view(rhsTermDG);
    auto h_rhsTermSEM = Kokkos::create_mirror_view(rhsTermSEM);
    auto h_rhsWeights = Kokkos::create_mirror_view(rhsWeights);

    // Source at the reference driver's default (srcx,srcy,srcz)=(1010,1010,1010).
    // Element index computed exactly like sem_proxy.cc InitSource():
    //   floor((rel_x*ex)/lx) + floor((rel_y*ey)/ly)*ex + floor((rel_z*ez)/lz)*ey*ex
    constexpr float kSrcX = 1010.0f, kSrcY = 1010.0f, kSrcZ = 1010.0f;
    const int src_elem = static_cast<int>(std::floor((kSrcX * kEx) / kDomainSize)) +
                         static_cast<int>(std::floor((kSrcY * kEy) / kDomainSize)) * kEx +
                         static_cast<int>(std::floor((kSrcZ * kEz) / kDomainSize)) * kEy * kEx;
    h_rhsElement(0) = src_elem;

    SolverUtils utils;
    float const tpeak = 1.0f / kF0;
    std::vector<float> sourceTerm = utils.computeSourceTerm(params.n_steps, kDt, kF0, 2, tpeak);
    for (int j = 0; j < params.n_steps; ++j) {
      h_rhsTermDG(0, j) = 0.0f;
      h_rhsTermSEM(0, j) = sourceTerm[j];
    }

    // Lagrange interpolation weights at the source point for each source.
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
      // Source point = the reference driver's src_coord_ (1010,1010,1010).
      std::array<float, 3> srcCoord{kSrcX, kSrcY, kSrcZ};
      SourceAndReceiverUtils::ComputeRHSWeights<kOrder>(cornerCoords, srcCoord, h_rhsWeights);
    }
    Kokkos::deep_copy(rhsElement, h_rhsElement);
    Kokkos::deep_copy(rhsTermDG, h_rhsTermDG);
    Kokkos::deep_copy(rhsTermSEM, h_rhsTermSEM);
    Kokkos::deep_copy(rhsWeights, h_rhsWeights);
    FENCE

    solver::fe::DGSEMWavefieldAcoustic wavefield(pnDGPrev, pnDGCurr, pnSEMPrev, pnSEMCurr);
    solver::fe::DGSEMRhsAcoustic rhs(rhsTermDG, rhsTermSEM, rhsElement, rhsWeights);
    solver::fe::DGSEMsolverData data(wavefield, rhs);

    // Deterministic persistent memory footprint (bytes).
    const size_t dg_wavefield = 2ull * n_elem * kNDof * sizeof(float);
    const size_t sem_wavefield = 2ull * n_node * sizeof(float);
    const size_t dg_scratch = 4ull * n_elem * kNDof * sizeof(float);  // rhs/mass/stiff/damp local
    const size_t sem_work = 1ull * n_node * sizeof(float);            // workVectorsGlobal_[0]
    const size_t mem_bytes = dg_wavefield + sem_wavefield + dg_scratch + sem_work;

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

    // Final DG field → host, for the reference comparison.
    auto h_pn = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, pnDGPrev);
    std::vector<double> field(n_elem * kNDof);
    for (int e = 0; e < n_elem; ++e)
      for (int d = 0; d < kNDof; ++d) field[e * kNDof + d] = static_cast<double>(h_pn(e, d));

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
    std::printf("output_file=%s\n", params.output_file.c_str());
    std::printf("PROXY_METRICS_END\n");
  }
  Kokkos::finalize();
  return 0;
}
