#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include "sem2d_benchmark_common.h"
#include "three_layer_case.h"

namespace {
using C = solver::validation::ThreeLayerResult;
using LTS = solver::SEMsolver2DLocalTimeStep<2, false>;
using Global = solver::SEMsolver2D<2, false>;
using Clock = std::chrono::steady_clock;
using namespace solver::validation::benchmark;

void exportDiagnostics(model::ModelUnstruct<float, int>& mesh, const vectorReal& initial, const Parameters& p,
                       const std::filesystem::path& output, const std::vector<float>& mass) {
  int steps = std::lround(p.duration / p.macro_dt);
  if (std::abs(steps * p.macro_dt - p.duration) > 1e-6) throw std::runtime_error("Unaligned diagnostic duration");
  Run local(mesh, initial, true, p.macro_dt, p),
      reference(mesh, initial, false, p.macro_dt / 32, p, p.double_reference),
      unsafe_global(mesh, initial, false, p.macro_dt, p);
  std::ofstream snapshots(output / "snapshots.csv"), receivers(output / "receivers.csv"),
      history(output / "history.csv");
  for (auto* f : {&snapshots, &receivers, &history}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  snapshots << "time,node,x,z,local,reference\n";
  receivers << "time,receiver,x,z,local,reference\n";
  history << "time,relative_l2_error,local_energy,global_macro_amplitude,reference_amplitude,cumulative_dissipation,"
             "energy_balance_error\n";
  double initial_energy = static_cast<LTS&>(*local.integrator).discreteEnergy(local.data), max_drift = 0;
  double max_error = 0, max_global_amplitude = 0, lost = 0, max_balance_error = 0, final_energy = initial_energy;
  double max_surface_pressure = 0;
  auto damping = solver::validation::layerHost(local.integrator->getDampingMatrix(0));
  if (!(initial_energy > 0)) throw std::runtime_error("Nonpositive initial diagnostic energy");
  for (int step = 0; step <= steps; ++step) {
    std::vector<float> previous;
    if (step) {
      previous = solver::validation::layerHost(local.data.getPreviousField(0));
      local.advance(step - 1);
      unsafe_global.advance(step - 1);
      for (int sub = 0; sub < 32; ++sub) reference.advance((step - 1) * 32 + sub);
    }
    auto a = solver::validation::layerHost(local.data.getCurrentField(0));
    auto b = solver::validation::layerHost(reference.data.getCurrentField(0));
    auto g = solver::validation::layerHost(unsafe_global.data.getCurrentField(0));
    double error = 0, norm = 0, ag = 0, ar = 0;
    for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
      if (!std::isfinite(a[i]) || !std::isfinite(b[i])) throw std::runtime_error("Diagnostic local/reference diverged");
      error += static_cast<double>(mass[i]) * (a[i] - b[i]) * (a[i] - b[i]);
      norm += static_cast<double>(mass[i]) * b[i] * b[i];
      ag = !std::isfinite(g[i]) ? std::numeric_limits<double>::infinity()
                                : std::max(ag, static_cast<double>(std::abs(g[i])));
      ar = std::max(ar, static_cast<double>(std::abs(b[i])));
      if (p.free_top && i >= mesh.getNumberOfNodes() - C::nx) {
        max_surface_pressure =
            std::max(max_surface_pressure, static_cast<double>(std::max(std::abs(a[i]), std::abs(b[i]))));
        if (a[i] != 0 || b[i] != 0) throw std::runtime_error("Nonzero pressure on the free surface");
      }
    }
    double energy = static_cast<LTS&>(*local.integrator).discreteEnergy(local.data);
    if (step)
      for (size_t i = 0; i < a.size(); ++i) {
        double velocity = (static_cast<double>(a[i]) - previous[i]) / (2 * p.macro_dt);
        lost += p.macro_dt * damping[i] * velocity * velocity;
      }
    double balance_error = std::abs(energy + lost - initial_energy) / initial_energy;
    max_balance_error = std::max(max_balance_error, balance_error);
    if (!std::isfinite(energy) || energy < -1e-6 * initial_energy)
      throw std::runtime_error("Invalid diagnostic energy");
    final_energy = energy;
    max_drift = std::max(max_drift, std::abs((energy - initial_energy) / initial_energy));
    max_error = std::max(max_error, std::sqrt(error / norm));
    max_global_amplitude = std::max(max_global_amplitude, ag);
    history << step * p.macro_dt << ',' << std::sqrt(error / norm) << ',' << energy << ',' << ag << ',' << ar << ','
            << lost << ',' << balance_error << '\n';
    for (int receiver = 0; receiver < 51; ++receiver) {
      int const ix = 10 + 2 * receiver;
      int const row = C::nz - 1 - std::lround(p.receiver_depth / (C::h / 2));
      int const node = ix + row * C::nx;
      receivers << step * p.macro_dt << ',' << receiver << ',' << mesh.nodeCoord(node, 0) << ','
                << mesh.nodeCoord(node, 2) << ',' << a[node] << ',' << b[node] << '\n';
    }
    if (step == 0 || step == steps / 6 || step == steps / 3 || step == steps / 2 || step == 2 * steps / 3 ||
        step == steps)
      for (int i = 0; i < mesh.getNumberOfNodes(); ++i)
        snapshots << step * p.macro_dt << ',' << i << ',' << mesh.nodeCoord(i, 0) << ',' << mesh.nodeCoord(i, 2) << ','
                  << a[i] << ',' << b[i] << '\n';
  }
  std::ofstream diagnostics(output / "diagnostics.json");
  diagnostics.exceptions(std::ios::failbit | std::ios::badbit);
  diagnostics << std::setprecision(12) << "{\"max_relative_energy_drift\": " << max_drift
              << ", \"max_relative_energy_balance_error\": " << max_balance_error
              << ", \"final_energy_fraction\": " << final_energy / initial_energy
              << ", \"dissipated_energy_fraction\": " << lost / initial_energy
              << ", \"free_surface_max_pressure\": " << max_surface_pressure
              << ", \"max_trajectory_relative_l2_error\": " << max_error
              << ", \"global_macro_exceeded_1e5\": " << (max_global_amplitude > 1e5 ? "true" : "false") << "}\n";
  if (max_balance_error > 1e-4) throw std::runtime_error("Excessive LTS energy balance error");
}

void benchmark(const std::filesystem::path& output, int repeats) {
  Parameters p;
  std::ifstream parameters(output / "parameters.txt");
  if (!(parameters >> p.x >> p.z >> p.sigma >> p.duration >> p.macro_dt >> p.cfl) ||
      !(p.sigma > 0 && p.duration > 0 && p.macro_dt > 0 && p.cfl > 0))
    throw std::runtime_error("Invalid parameters");
  auto d = solver::validation::threeLayerMesh();
  int boundary_mode = 0;
  if (parameters >> boundary_mode) {
    if (boundary_mode < 0 || boundary_mode > 3) throw std::runtime_error("Invalid boundary mode");
    p.absorbing = boundary_mode == 1 || boundary_mode == 2;
    p.free_top = boundary_mode == 2 || boundary_mode == 3;
  } else if (!parameters.eof())
    throw std::runtime_error("Invalid boundary mode");
  int reference_mode = 0;
  if (parameters >> reference_mode) {
    if (reference_mode != 0 && reference_mode != 1) throw std::runtime_error("Invalid reference precision mode");
    p.double_reference = reference_mode == 1;
  } else if (!parameters.eof())
    throw std::runtime_error("Invalid reference precision mode");
  float depth = 0;
  if (parameters >> depth) {
    if (!std::isfinite(depth) || depth < 0 || depth > C::ez * C::h ||
        std::abs(depth / (C::h / 2) - std::round(depth / (C::h / 2))) > 1e-5)
      throw std::runtime_error("Receiver depth must coincide with a GLL row");
    p.receiver_depth = depth;
  } else if (!parameters.eof())
    throw std::runtime_error("Invalid receiver depth");
  if (p.free_top && (p.receiver_depth == 0 || p.z >= C::ez * C::h))
    throw std::runtime_error("Pressure source and receivers must be below the pressure-free top");
  if (p.absorbing || p.free_top) {
    d.boundaries_t_ = allocateVector<dualVectorInt>(d.n_node_, "inclusion absorbing boundary");
    for (int j = 0; j < C::nz; ++j)
      for (int i = 0; i < C::nx; ++i)
        d.boundaries_t_.view_host()(i + j * C::nx) =
            static_cast<int>(p.free_top && j == C::nz - 1 ? model::BoundaryFlag::Surface
                             : p.absorbing && (i == 0 || j == 0 || i == C::nx - 1 || j == C::nz - 1)
                                 ? model::BoundaryFlag::Damping
                                 : model::BoundaryFlag::InteriorNode);
    d.boundaries_t_.modify_host();
  }
  std::ifstream input(output / "model.csv");
  std::string line;
  std::getline(input, line);
  int count = 0;
  while (std::getline(input, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int element, region;
    float vp, rho, target_dt;
    if (!(row >> element >> vp >> rho >> region >> target_dt) || element != count || element >= d.n_element_ ||
        !(vp > 0 && rho > 0 && target_dt > 0) || !std::isfinite(vp + rho + target_dt))
      throw std::runtime_error("Invalid element row");
    d.model_vp_element_.view_host()(element) = vp;
    d.model_rho_element_.view_host()(element) = rho;
    ++count;
  }
  if (count != d.n_element_) throw std::runtime_error("Incomplete element model");
  d.model_vp_element_.modify_host();
  d.model_rho_element_.modify_host();
  model::ModelUnstruct<float, int> mesh(d);
  auto plan = solver::buildSEM2DLocalTimeStepPlan(mesh, p.cfl, p.macro_dt);
  vectorReal initial("inclusion Gaussian", mesh.getNumberOfNodes());
  auto h = Kokkos::create_mirror_view(initial);
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double x = (mesh.nodeCoord(i, 0) - p.x) / p.sigma, z = (mesh.nodeCoord(i, 2) - p.z) / p.sigma;
    h(i) = std::exp(-0.5 * (x * x + z * z));
    if (p.free_top) {
      double image_z = (mesh.nodeCoord(i, 2) - (2 * C::ez * C::h - p.z)) / p.sigma;
      h(i) -= std::exp(-.5 * (x * x + image_z * image_z));
      if (mesh.isFreeSurface(i)) h(i) = 0;
    }
  }
  Kokkos::deep_copy(initial, h);
  Global mass_solver;
  mass_solver.computeFEInit(mesh, {0, 0, 0}, false, 0);
  auto mass = solver::validation::layerHost(mass_solver.getMassMatrixAcoustic());
  int reference_steps = std::lround(p.duration / (p.macro_dt / 32));
  auto reference = run({"reference", false, reference_steps}, mesh, initial, p, true).field;
  double norm = 0;
  for (size_t i = 0; i < mass.size(); ++i) norm += static_cast<double>(mass[i]) * reference[i] * reference[i];
  auto error = [&](const std::vector<float>& values) {
    double square = 0;
    for (size_t i = 0; i < mass.size(); ++i) {
      if (!std::isfinite(values[i])) return std::numeric_limits<double>::infinity();
      double diff = static_cast<double>(values[i]) - reference[i];
      square += mass[i] * diff * diff;
    }
    return std::sqrt(square / norm);
  };
  double reference_delta = error(run({"reference check", false, reference_steps / 2}, mesh, initial, p, true).field);
  if (!(reference_delta < .001))
    throw std::runtime_error("Reference refinement discrepancy exceeds 0.1%: " + std::to_string(reference_delta));
  double limit = spectralEstimate(mesh, output);
  std::vector<Case> cases;
  for (bool local : {false, true})
    for (int steps : {100, 120, 150, 180, 200, 225, 240, 270, 300, 360, 450, 600, 1200})
      cases.push_back({(local ? "LTS_" : "Global_") + std::to_string(std::lround(steps * p.duration / 1.8f)), local,
                       static_cast<int>(std::lround(steps * p.duration / 1.8f))});
  for (int w = 0; w < 2; ++w)
    for (auto const& c : cases) run(c, mesh, initial, p);
  std::ofstream raw(output / "timing_samples.csv"), summary(output / "timing_summary.csv");
  for (auto* f : {&raw, &summary}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  raw << "case,local,dt_ms,steps,trial,setup_ms,loop_ms,relative_l2_error\n";
  summary << "case,local,dt_ms,steps,setup_median_ms,loop_median_ms,loop_q25_ms,loop_q75_ms,total_median_ms,relative_"
             "l2_error\n";
  std::vector<std::vector<double>> setup(cases.size()), loop(cases.size()), total(cases.size()), errors(cases.size());
  std::vector<int> order(cases.size());
  std::iota(order.begin(), order.end(), 0);
  std::mt19937 rng(20261109);
  for (int trial = 0; trial < repeats; ++trial) {
    std::shuffle(order.begin(), order.end(), rng);
    for (int j : order) {
      auto const& c = cases[j];
      auto s = run(c, mesh, initial, p);
      double e = error(s.field);
      raw << c.name << ',' << c.local << ',' << p.duration * 1000 / c.steps << ',' << c.steps << ',' << trial << ','
          << s.setup_ms << ',' << s.loop_ms << ',' << e << '\n';
      setup[j].push_back(s.setup_ms);
      loop[j].push_back(s.loop_ms);
      total[j].push_back(s.setup_ms + s.loop_ms);
      errors[j].push_back(e);
    }
  }
  for (size_t j = 0; j < cases.size(); ++j) {
    auto const& c = cases[j];
    summary << c.name << ',' << c.local << ',' << p.duration * 1000 / c.steps << ',' << c.steps << ','
            << quantile(setup[j], .5) << ',' << quantile(loop[j], .5) << ',' << quantile(loop[j], .25) << ','
            << quantile(loop[j], .75) << ',' << quantile(total[j], .5) << ',' << quantile(errors[j], .5) << '\n';
  }
  exportDiagnostics(mesh, initial, p, output, mass);
  std::ofstream protocol(output / "protocol.json");
  protocol.exceptions(std::ios::failbit | std::ios::badbit);
  protocol << std::setprecision(12) << "{\"backend\": \"" << Kokkos::DefaultExecutionSpace::name()
           << "\", \"repetitions\": " << repeats << ", \"warmups\": 2, \"seed\": 20261109, "
           << "\"fine_elements\": " << plan.fine_elements.size()
           << ", \"support_elements\": " << plan.fine_support_elements.size()
           << ", \"reference_delta_l2\": " << reference_delta << ", \"global_dt_limit_estimate_s\": " << limit
           << ", \"absorbing_boundaries\": " << (p.absorbing ? "true" : "false")
           << ", \"free_surface_top\": " << (p.free_top ? "true" : "false")
           << ", \"receiver_depth_m\": " << p.receiver_depth << ", \"source_kind\": \""
           << (p.free_top ? "initial_gaussian_dirichlet_image" : "initial_gaussian")
           << "\", \"reference_dt_s\": " << p.duration / reference_steps
           << ", \"reference_double_time_accumulation\": " << (p.double_reference ? "true" : "false")
           << ", \"target_steps_are_executed\": false, \"loop_excludes\": \"setup, diagnostics, host copies, IO\"}\n";
  std::cout << output.filename() << ": " << plan.fine_elements.size() << " fine elements, global spectral dt estimate "
            << limit * 1000 << " ms, reference delta " << reference_delta << std::endl;
}
}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int status = 0;
  try {
    if (argc < 2) throw std::invalid_argument("Usage: inclusion_benchmark prepared_case_directory [repetitions]");
    int repeats = argc > 2 ? std::stoi(argv[2]) : 7;
    if (repeats < 3) throw std::invalid_argument("Use at least three repetitions");
    benchmark(argv[1], repeats);
  } catch (std::exception const& e) {
    std::cerr << e.what() << '\n';
    status = 1;
  }
  Kokkos::finalize();
  return status;
}
