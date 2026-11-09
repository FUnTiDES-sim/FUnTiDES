// Repeated propagation-only timings of the already validated four-run case.
#include <iostream>
#include <map>
#include <sstream>

#include "gmsh_2d_reader.h"
#include "sem2d_benchmark_common.h"

namespace {
using namespace solver::validation::benchmark;
using Mesh = model::ModelUnstruct<float, int>;

Mesh readMesh(const std::filesystem::path& root) {
  std::ifstream in(root / "materials.csv");
  std::string line;
  std::getline(in, line);
  std::map<int, model::GmshAcousticMaterial> materials;
  while (std::getline(in, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int tag;
    float vp, rho;
    if (!(row >> tag >> vp >> rho)) throw std::runtime_error("Invalid materials");
    materials.emplace(tag, model::GmshAcousticMaterial{vp, rho});
  }
  return Mesh(model::readGmsh2D(root / "bilayer.msh", 2, materials,
                              {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}}));
}

vectorReal initialField(Mesh& mesh) {
  vectorReal result("timing common initial field", mesh.getNumberOfNodes());
  auto h = Kokkos::create_mirror_view(result);
  double pi = std::acos(-1.);
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double x = mesh.nodeCoord(i, 0), z = mesh.nodeCoord(i, 2);
    double gaussian = std::exp(-.5 * (std::pow((x - 2200) / 160, 2) + std::pow((z - 2800) / 160, 2))) -
                      std::exp(-.5 * (std::pow((x - 2200) / 160, 2) + std::pow((z - 3200) / 160, 2)));
    double noise = 0;
    if (x > 2990 && x < 3010 && z > 890 && z < 910)
      noise = 1e-6 * std::pow(std::sin(pi * (x - 2990) / 20), 2) *
              std::pow(std::sin(pi * (z - 890) / 20), 2) * std::cos(2 * pi * (x - 3000)) *
              std::cos(2 * pi * (z - 900));
    h(i) = mesh.isFreeSurface(i) ? 0 : gaussian + noise;
  }
  Kokkos::deep_copy(result, h);
  return result;
}

double maxAmplitude(const std::vector<float>& field) {
  double maximum = 0;
  for (float v : field) {
    if (!std::isfinite(v)) throw std::runtime_error("Nonfinite timed field");
    maximum = std::max(maximum, static_cast<double>(std::abs(v)));
  }
  return maximum;
}

double error(const std::vector<float>& field, const std::vector<float>& reference, const std::vector<float>& mass) {
  double square = 0, norm = 0;
  for (size_t i = 0; i < field.size(); ++i) {
    double d = static_cast<double>(field[i]) - reference[i];
    square += mass[i] * d * d;
    norm += static_cast<double>(mass[i]) * reference[i] * reference[i];
  }
  if (!(norm > 0) || !std::isfinite(square)) throw std::runtime_error("Invalid timed error norm");
  return std::sqrt(square / norm);
}

void execute(const std::filesystem::path& baseline, const std::filesystem::path& root,
             const std::filesystem::path& parameters_file, int repeats) {
  // Values are exported from the passed test_result/calibration, never reselected.
  std::ifstream params(parameters_file);
  float dt, cfl;
  int steps, divergent_steps;
  if (!(params >> dt >> cfl >> steps >> divergent_steps) || dt <= 0 || cfl <= 0 || steps < 1 || divergent_steps < 1)
    throw std::runtime_error("Expected timing parameters: macro_dt cfl macro_steps divergent_steps");
  auto coarse = readMesh(baseline), refined = readMesh(root);
  auto coarse_initial = initialField(coarse), refined_initial = initialField(refined);
  Parameters p{};
  p.duration = dt * steps;
  p.macro_dt = dt;
  p.cfl = cfl;
  p.absorbing = p.free_top = true;
  Global mass_solver;
  mass_solver.computeFEInit(refined, {0, 0, 0}, false, 0);
  auto mass = hostField(mass_solver.getMassMatrixAcoustic());
  std::vector<float> reference(refined.getNumberOfNodes());
  std::vector<bool> seen(reference.size());
  std::ifstream snapshots(root / "snapshots.csv");
  std::string line;
  std::getline(snapshots, line);
  while (std::getline(snapshots, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    double time, x, z;
    float local, ref;
    int node;
    if (!(row >> time >> node >> x >> z >> local >> ref) || node < 0 || node >= static_cast<int>(reference.size()))
      throw std::runtime_error("Invalid validation snapshot");
    if (std::abs(time - p.duration) < 1e-5) reference[node] = ref, seen[node] = true;
  }
  if (std::find(seen.begin(), seen.end(), false) != seen.end()) throw std::runtime_error("Missing final reference");
  struct Scenario { const char* name; bool local, control, divergent; float dt; int steps; };
  std::vector<Scenario> scenarios{{"Control_global_macro", false, true, false, dt, steps},
                                 {"Refined_global_macro_unstable", false, false, true, dt, divergent_steps},
                                 {"Refined_global_fine", false, false, false, dt / 2, 2 * steps},
                                 {"Refined_local", true, false, false, dt, steps}};
  struct Measurement { double setup, loop, amplitude, relative_error; };
  auto sample = [&](const Scenario& s) {
    Mesh& mesh = s.control ? coarse : refined;
    auto const& field = s.control ? coarse_initial : refined_initial;
    Kokkos::fence();
    auto start = Clock::now();
    Run simulation(mesh, field, s.local, s.dt, p);
    double setup = milliseconds(start);
    Kokkos::fence();
    start = Clock::now();
    for (int i = 0; i < s.steps; ++i) simulation.advance(i);
    Kokkos::fence();
    double loop = milliseconds(start);
    // Copies and numerical checks occur after the propagation timer stops.
    auto final = hostField(simulation.data.getCurrentField(0));
    double amplitude = maxAmplitude(final) / maxAmplitude(hostField(field));
    if (s.divergent ? amplitude <= 1e5 : amplitude >= 10) throw std::runtime_error("Timed stability witness changed");
    double relative_error = s.control || s.divergent ? -1 : error(final, reference, mass);
    if (relative_error >= .01) throw std::runtime_error("Timed numerical error exceeds validation budget");
    return Measurement{setup, loop, amplitude, relative_error};
  };
  for (int warmup = 0; warmup < 2; ++warmup)
    for (auto const& s : scenarios) sample(s);
  std::ofstream raw(root / "timing_samples.csv"), summary(root / "timing_summary.csv");
  for (auto* f : {&raw, &summary}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  raw << "case,trial,dt_s,steps,physical_duration_s,setup_ms,loop_ms,amplitude_ratio,relative_l2_error\n";
  summary << "case,dt_s,steps,physical_duration_s,setup_median_ms,loop_median_ms,loop_q25_ms,loop_q75_ms,solver_total_median_ms,relative_l2_error\n";
  std::vector<std::vector<double>> setup(4), loop(4), total(4), errors(4);
  std::vector<int> order{0, 1, 2, 3};
  std::mt19937 rng(20261108);
  for (int trial = 0; trial < repeats; ++trial) {
    std::shuffle(order.begin(), order.end(), rng);
    for (int j : order) {
      auto const& s = scenarios[j];
      auto m = sample(s);
      raw << s.name << ',' << trial << ',' << s.dt << ',' << s.steps << ',' << s.dt * s.steps << ','
          << m.setup << ',' << m.loop << ',' << m.amplitude << ',';
      if (m.relative_error >= 0) raw << m.relative_error;
      raw << '\n';
      setup[j].push_back(m.setup);
      loop[j].push_back(m.loop);
      total[j].push_back(m.setup + m.loop);
      errors[j].push_back(m.relative_error);
    }
    raw.flush();
    std::cout << "Timing repeat " << trial + 1 << '/' << repeats << std::endl;
  }
  for (int j = 0; j < 4; ++j) {
    auto const& s = scenarios[j];
    summary << s.name << ',' << s.dt << ',' << s.steps << ',' << s.dt * s.steps << ',' << quantile(setup[j], .5)
            << ',' << quantile(loop[j], .5) << ',' << quantile(loop[j], .25) << ',' << quantile(loop[j], .75)
            << ',' << quantile(total[j], .5) << ',';
    if (errors[j][0] >= 0) summary << quantile(errors[j], .5);
    summary << '\n';
    std::cout << s.name << ": " << quantile(loop[j], .5) << " ms, physical duration " << s.dt * s.steps << " s" << std::endl;
  }
}
}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int status = 1;
  try {
    if (argc < 4) throw std::invalid_argument("Usage: sem2d_gmsh_inclusion_timing baseline_dir result_dir parameters_file [repetitions=7]");
    int repeats = argc > 4 ? std::stoi(argv[4]) : 7;
    if (repeats < 3) throw std::invalid_argument("At least 3 repetitions required");
    execute(argv[1], argv[2], argv[3], repeats);
    status = 0;
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
  }
  Kokkos::finalize();
  return status;
}
