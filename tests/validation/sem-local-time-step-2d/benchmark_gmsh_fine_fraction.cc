// Fine-fraction study of the two-level local time step (plan step 2) on gmsh meshes with a
// geometrically refined patch (generate_gmsh_fine_fraction.py). For each mesh:
//  - global step dt_g = margin * 2 / sqrt(lambda_max(M^-1 K)) (power iteration);
//  - LTS macro step 2 dt_g; partition by the geometric estimate scaled so that the smallest
//    element sits exactly at the fine step: fine elements are those whose estimate is below
//    twice the smallest one;
//  - stability (amplitude, LTS discrete energy) and accuracy (final field against the global
//    scheme at dt_g / 4) checks, then repeated, interleaved timings of the propagation loops.
// Usage: sem2d_gmsh_fine_fraction output.csv repetitions duration_s case_dir...
#include <iostream>
#include <map>
#include <sstream>

#include "gmsh_2d_reader.h"
#include "sem2d_benchmark_common.h"

namespace {
using namespace solver::validation::benchmark;
using Mesh = model::ModelUnstruct<float, int>;

Mesh readCase(const std::filesystem::path& root) {
  std::ifstream in(root / "materials.csv");
  std::string line;
  std::getline(in, line);
  std::map<int, model::GmshAcousticMaterial> materials;
  while (std::getline(in, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int tag;
    float vp, rho;
    if (!(row >> tag >> vp >> rho)) throw std::runtime_error("Invalid materials in " + root.string());
    materials.emplace(tag, model::GmshAcousticMaterial{vp, rho});
  }
  return Mesh(model::readGmsh2D(root / "mesh.msh", 2, materials,
                              {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}}));
}

vectorReal initialField(Mesh& mesh) {
  vectorReal result("fine fraction initial field", mesh.getNumberOfNodes());
  auto h = Kokkos::create_mirror_view(result);
  double top = 0;
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) top = std::max(top, static_cast<double>(mesh.nodeCoord(i, 2)));
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double x = mesh.nodeCoord(i, 0), z = mesh.nodeCoord(i, 2);
    // Gaussian pulse at 1200 m depth, left of the patch, and its image above the free surface
    double zs = top - 1200;
    h(i) = mesh.isFreeSurface(i) ? 0
                                 : std::exp(-.5 * (std::pow((x - 2400) / 120, 2) + std::pow((z - zs) / 120, 2))) -
                                       std::exp(-.5 * (std::pow((x - 2400) / 120, 2) + std::pow((z - 2 * top + zs) / 120, 2)));
  }
  Kokkos::deep_copy(result, h);
  return result;
}

double maxAbs(const std::vector<float>& field) {
  double m = 0;
  for (float v : field) {
    if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
    m = std::max(m, static_cast<double>(std::abs(v)));
  }
  return m;
}

double relativeError(const std::vector<float>& field, const std::vector<float>& reference,
                     const std::vector<float>& mass) {
  double square = 0, norm = 0;
  for (size_t i = 0; i < field.size(); ++i) {
    double d = static_cast<double>(field[i]) - reference[i];
    square += mass[i] * d * d;
    norm += static_cast<double>(mass[i]) * reference[i] * reference[i];
  }
  return std::sqrt(square / norm);
}

struct Timed { double loop_ms; std::vector<float> field; double energy_start, energy_end; };

Timed propagate(Mesh& mesh, const vectorReal& initial, bool local, float dt, int steps, const Parameters& p,
                bool precomputed_metric = true) {
  Run run(mesh, initial, local, dt, p);
  if (local) static_cast<LTS&>(*run.integrator).setPrecomputedSupportMetric(precomputed_metric);
  double e0 = local ? static_cast<LTS&>(*run.integrator).discreteEnergy(run.data) : 0;
  Kokkos::fence();
  auto start = Clock::now();
  for (int i = 0; i < steps; ++i) run.advance(i);
  Kokkos::fence();
  double loop = milliseconds(start);
  double e1 = local ? static_cast<LTS&>(*run.integrator).discreteEnergy(run.data) : 0;
  return {loop, hostField(run.data.getCurrentField(0)), e0, e1};
}

void study(const std::filesystem::path& csv_path, int repeats, double duration,
           const std::vector<std::filesystem::path>& cases) {
  std::ofstream csv(csv_path);
  csv.exceptions(std::ios::failbit | std::ios::badbit);
  csv << std::setprecision(10)
      << "case,elements,nodes,fine_elements,support_elements,fine_fraction,support_fraction,global_limit_s,"
         "dt_global_s,macro_dt_s,global_steps,global_ms,local_ms,gain,global_ms_q25,global_ms_q75,local_ms_q25,"
         "local_ms_q75,error_global,error_local,energy_ratio_local,amplitude_ratio_local,local_same_kernel_ms,"
         "gain_same_kernel,stable\n";
  for (auto const& root : cases) {
    auto mesh = readCase(root);
    auto initial = initialField(mesh);
    double limit = spectralEstimate(mesh, root);
    int steps = 2 * static_cast<int>(std::ceil(duration / (0.9 * limit) / 2));   // even: whole macro steps
    float dt = static_cast<float>(duration / steps), macro = 2 * dt;
    auto probe = solver::buildSEM2DLocalTimeStepPlan(mesh, 1.f, 1e-9f);
    double min_estimate = std::numeric_limits<double>::infinity();
    for (auto const& s : probe.element_steps) min_estimate = std::min(min_estimate, static_cast<double>(s.estimate));
    Parameters p{};
    p.macro_dt = macro;
    p.cfl = static_cast<float>(dt / min_estimate * 1.0001);   // smallest element exactly at the fine step
    auto plan = solver::buildSEM2DLocalTimeStepPlan(mesh, p.cfl, macro);
    int const n_elem = mesh.getNumberOfElements();
    std::cout << root.filename().string() << ": " << n_elem << " elements, " << plan.fine_elements.size() << " fine ("
              << 100.0 * plan.fine_elements.size() / n_elem << " %), support " << plan.fine_support_elements.size()
              << ", global limit " << limit * 1e3 << " ms, dt " << dt * 1e3 << " ms, " << steps << " steps"
              << std::endl;

    // Checks (outside the timings)
    Global mass_solver;
    mass_solver.computeFEInit(mesh, {0, 0, 0}, false, 0);
    auto mass = hostField(mass_solver.getMassMatrixAcoustic());
    auto reference = propagate(mesh, initial, false, dt / 4, 4 * steps, p).field;
    auto global = propagate(mesh, initial, false, dt, steps, p);
    auto local = propagate(mesh, initial, true, macro, steps / 2, p);
    double amplitude = maxAbs(local.field) / maxAbs(hostField(initial));
    double energy_ratio = local.energy_end / local.energy_start;
    if (!(amplitude < 10) || !(local.energy_start > 0) || !(energy_ratio <= 1.0001)) {
      // Recorded, not fatal: the geometric partition does not certify the spectrum (plan step 3)
      std::cout << "  LTS unstable (amplitude ratio " << amplitude << ", energy ratio " << energy_ratio << ")"
                << std::endl;
      csv << root.filename().string() << ',' << n_elem << ',' << mesh.getNumberOfNodes() << ','
          << plan.fine_elements.size() << ',' << plan.fine_support_elements.size() << ','
          << static_cast<double>(plan.fine_elements.size()) / n_elem << ','
          << static_cast<double>(plan.fine_support_elements.size()) / n_elem << ',' << limit << ',' << dt << ','
          << macro << ',' << steps << ",,,,,,,,,," << energy_ratio << ',' << amplitude << ",,,0\n";
      csv.flush();
      continue;
    }
    double error_global = relativeError(global.field, reference, mass);
    double error_local = relativeError(local.field, reference, mass);

    std::vector<double> tg, tl, ts;
    propagate(mesh, initial, true, macro, steps / 2, p);   // warm-up
    for (int r = 0; r < repeats; ++r) {
      // rotated order: global, LTS, LTS with the sem2d element kernel
      for (int k = 0; k < 3; ++k) switch ((k + r) % 3) {
          case 0: tg.push_back(propagate(mesh, initial, false, dt, steps, p).loop_ms); break;
          case 1: tl.push_back(propagate(mesh, initial, true, macro, steps / 2, p).loop_ms); break;
          default: ts.push_back(propagate(mesh, initial, true, macro, steps / 2, p, false).loop_ms);
        }
    }
    double g = quantile(tg, .5), l = quantile(tl, .5), same = quantile(ts, .5);
    csv << root.filename().string() << ',' << n_elem << ',' << mesh.getNumberOfNodes() << ','
        << plan.fine_elements.size() << ',' << plan.fine_support_elements.size() << ','
        << static_cast<double>(plan.fine_elements.size()) / n_elem << ','
        << static_cast<double>(plan.fine_support_elements.size()) / n_elem << ',' << limit << ',' << dt << ','
        << macro << ',' << steps << ',' << g << ',' << l << ',' << g / l << ',' << quantile(tg, .25) << ','
        << quantile(tg, .75) << ',' << quantile(tl, .25) << ',' << quantile(tl, .75) << ',' << error_global << ','
        << error_local << ',' << energy_ratio << ',' << amplitude << ',' << same << ',' << g / same << ",1\n";
    csv.flush();
    std::cout << "  global " << g << " ms, local " << l << " ms, gain " << g / l << " (same kernel " << g / same
              << "), errors " << error_global
              << " / " << error_local << std::endl;
  }
}
}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int status = 0;
  try {
    if (argc < 5)
      throw std::invalid_argument("Usage: sem2d_gmsh_fine_fraction output.csv repetitions duration_s case_dir...");
    std::vector<std::filesystem::path> cases(argv + 4, argv + argc);
    study(argv[1], std::stoi(argv[2]), std::stod(argv[3]), cases);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    status = 1;
  }
  Kokkos::finalize();
  return status;
}
