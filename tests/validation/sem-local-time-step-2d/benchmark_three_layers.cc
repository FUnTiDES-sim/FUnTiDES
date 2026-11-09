#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>

#include "three_layer_case.h"

namespace {

using Clock = std::chrono::steady_clock;
using C = solver::validation::ThreeLayerResult;

struct Case {
  const char* name;
  bool local;
  float dt;
  int steps;
};

struct Measurement {
  double setup_ms, loop_ms, relative_error;
  std::vector<float> final;
};

double milliseconds(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

Measurement run(const Case& config, model::ModelUnstruct<float, int>& mesh, const vectorReal& initial) {
  Kokkos::fence();
  auto start = Clock::now();
  std::unique_ptr<solver::Solver> integrator;
  if (config.local)
    integrator = std::make_unique<solver::SEMsolver2DLocalTimeStep<2, false>>();
  else
    integrator = std::make_unique<solver::SEMsolver2D<2, false>>();
  integrator->computeFEInit(mesh, {0, 0, 0}, false, 0);
  if (config.local)
    static_cast<solver::SEMsolver2DLocalTimeStep<2, false>&>(*integrator).configureLocalTimeStep(config.dt, C::cfl);

  int const n = mesh.getNumberOfNodes();
  vectorReal current("benchmark current", n), previous("benchmark previous", n);
  Kokkos::deep_copy(current, initial);
  solver::RhsAcoustic empty(arrayReal("benchmark no source", 0, 0), vectorInt("benchmark no element", 0),
                            arrayReal("benchmark no weights", 0, 0));
  solver::SEMsolverDataAcoustic data(solver::WavefieldAcoustic(previous, current), empty);
  if (config.local) {
    static_cast<solver::SEMsolver2DLocalTimeStep<2, false>&>(*integrator).initializePreviousField(data);
  } else {
    integrator->computeForces(config.dt, 0, data);
    auto mass = integrator->getMassMatrixAcoustic();
    auto force = integrator->getForceVector(0);
    float const dt = config.dt;
    Kokkos::parallel_for(
        "benchmark initial previous", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int const i) { previous(i) = current(i) - 0.5f * dt * dt * force(i) / mass(i); });
  }
  Kokkos::fence();
  double const setup = milliseconds(start);

  // Time stepping only: no receiver sampling, energy, CSV or host copies here.
  start = Clock::now();
  for (int step = 0; step < config.steps; ++step) {
    integrator->computeOneStep(config.dt, step, data);
    data.swapWavefields();
  }
  Kokkos::fence();
  double const loop = milliseconds(start);
  auto final = solver::validation::layerHost(data.getCurrentField(0));
  for (float p : final)
    if (!std::isfinite(p)) throw std::runtime_error(std::string(config.name) + ": nonfinite field");
  return {setup, loop, 0., std::move(final)};
}

double quantile(std::vector<double> values, double p) {
  std::sort(values.begin(), values.end());
  double const index = p * (values.size() - 1);
  auto const lo = static_cast<size_t>(index), hi = std::min(lo + 1, values.size() - 1);
  return values[lo] + (index - lo) * (values[hi] - values[lo]);
}

void benchmark(const std::filesystem::path& output, int repeats) {
  auto d = solver::validation::threeLayerMesh();
  model::ModelUnstruct<float, int> mesh(d);
  vectorReal initial("benchmark Gaussian", mesh.getNumberOfNodes());
  auto h_initial = Kokkos::create_mirror_view(initial);
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double const x = (mesh.nodeCoord(i, 0) - 3000) / 160.;
    double const z = (mesh.nodeCoord(i, 2) - 2500) / 160.;
    h_initial(i) = static_cast<float>(std::exp(-0.5 * (x * x + z * z)));
  }
  Kokkos::deep_copy(initial, h_initial);
  std::vector<Case> const cases{{"LTS_12_6ms", true, C::coarse_dt, 150},
                                {"Global_9ms", false, C::coarse_dt * 0.75f, 200},
                                {"Global_6ms", false, C::coarse_dt / 2, 300},
                                {"Global_3ms", false, C::coarse_dt / 4, 600},
                                {"Global_1.5ms", false, C::coarse_dt / 8, 1200}};
  Case const reference{"Reference_0.75ms", false, C::coarse_dt / 16, 2400};
  auto fine_reference = run(reference, mesh, initial).final;
  solver::SEMsolver2D<2, false> mass_solver;
  mass_solver.computeFEInit(mesh, {0, 0, 0}, false, 0);
  auto mass = solver::validation::layerHost(mass_solver.getMassMatrixAcoustic());
  double norm = 0;
  for (size_t i = 0; i < mass.size(); ++i) norm += static_cast<double>(mass[i]) * fine_reference[i] * fine_reference[i];
  auto error = [&](const std::vector<float>& values) {
    double square = 0;
    for (size_t i = 0; i < values.size(); ++i) {
      double const diff = static_cast<double>(values[i]) - fine_reference[i];
      square += mass[i] * diff * diff;
    }
    return std::sqrt(square / norm);
  };

  // Warm up each implementation and alternate trial order to limit time drift.
  for (int warmup = 0; warmup < 2; ++warmup)
    for (auto const& config : cases) run(config, mesh, initial);
  std::filesystem::create_directories(output);
  std::ofstream raw(output / "timing_samples.csv");
  raw.exceptions(std::ios::failbit | std::ios::badbit);
  raw << std::setprecision(12) << "case,trial,dt_ms,steps,setup_ms,loop_ms,relative_l2_error\n";
  std::vector<std::vector<Measurement>> measurements(cases.size());
  std::vector<int> order(cases.size());
  std::iota(order.begin(), order.end(), 0);
  std::mt19937 rng(20261106);
  for (int trial = 0; trial < repeats; ++trial) {
    std::shuffle(order.begin(), order.end(), rng);
    for (int index : order) {
      auto const& config = cases[index];
      auto sample = run(config, mesh, initial);
      sample.relative_error = error(sample.final);
      sample.final.clear();
      raw << config.name << ',' << trial << ',' << config.dt * 1000 << ',' << config.steps << ',' << sample.setup_ms
          << ',' << sample.loop_ms << ',' << sample.relative_error << '\n';
      measurements[index].push_back(std::move(sample));
    }
  }
  std::ofstream summary(output / "timing_summary.csv");
  summary.exceptions(std::ios::failbit | std::ios::badbit);
  summary
      << std::setprecision(12)
      << "case,dt_ms,steps,setup_median_ms,loop_median_ms,loop_q25_ms,loop_q75_ms,total_median_ms,relative_l2_error\n";
  for (size_t j = 0; j < cases.size(); ++j) {
    std::vector<double> setup, loop, total;
    for (auto const& m : measurements[j]) {
      setup.push_back(m.setup_ms);
      loop.push_back(m.loop_ms);
      total.push_back(m.setup_ms + m.loop_ms);
    }
    summary << cases[j].name << ',' << cases[j].dt * 1000 << ',' << cases[j].steps << ',' << quantile(setup, .5) << ','
            << quantile(loop, .5) << ',' << quantile(loop, .25) << ',' << quantile(loop, .75) << ','
            << quantile(total, .5) << ',' << measurements[j].front().relative_error << '\n';
    std::cout << cases[j].name << ": loop " << quantile(loop, .5) << " ms; setup " << quantile(setup, .5)
              << " ms; error " << 100 * measurements[j].front().relative_error << " %\n";
  }
  std::ofstream metadata(output / "timing_protocol.json");
  metadata.exceptions(std::ios::failbit | std::ios::badbit);
  metadata << "{\n  \"backend\": \"" << Kokkos::DefaultExecutionSpace::name() << "\",\n"
           << "  \"concurrency\": " << Kokkos::DefaultExecutionSpace().concurrency() << ",\n"
           << "  \"repetitions\": " << repeats << ", \"warmups_per_case\": 2,\n"
           << "  \"duration_s\": 1.8, \"reference_for_error_dt_ms\": 0.75,\n"
           << "  \"elements\": 1800, \"nodes\": 7381, \"order\": 2,\n"
           << "  \"loop_excludes\": \"setup, energy, sampling, file IO, host copies\"\n}\n";
}

}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int result = 0;
  try {
    auto const out = argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("validation-timings");
    int const repeats = argc > 2 ? std::stoi(argv[2]) : 11;
    if (repeats < 3) throw std::invalid_argument("Use at least three repetitions");
    benchmark(out, repeats);
  } catch (std::exception const& e) {
    std::cerr << e.what() << '\n';
    result = 1;
  }
  Kokkos::finalize();
  return result;
}
