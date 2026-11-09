// Performance/accuracy study: fixed geometry, variable thickness of the fast layer.
// Deliberately keeps the production integrators and original validation unchanged.
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

#include "three_layer_case.h"

namespace {
using Clock = std::chrono::steady_clock;
using C = solver::validation::ThreeLayerResult;
struct Case {
  std::string name;
  bool local;
  int steps;
  float dt() const { return 1.8f / steps; }
};
struct Sample {
  double setup_ms, loop_ms;
  std::vector<float> field;
};
double ms(Clock::time_point start) { return std::chrono::duration<double, std::milli>(Clock::now() - start).count(); }
Sample run(const Case& config, model::ModelUnstruct<float, int>& mesh, const vectorReal& initial) {
  Kokkos::fence();
  auto start = Clock::now();
  std::unique_ptr<solver::Solver> integrator;
  if (config.local)
    integrator = std::make_unique<solver::SEMsolver2DLocalTimeStep<2, false>>();
  else
    integrator = std::make_unique<solver::SEMsolver2D<2, false>>();
  integrator->computeFEInit(mesh, {0, 0, 0}, false, 0);
  // Scaling cfl with dt preserves exactly the same partition in the LTS dt sweep.
  if (config.local)
    static_cast<solver::SEMsolver2DLocalTimeStep<2, false>&>(*integrator)
        .configureLocalTimeStep(config.dt(), C::cfl * config.dt() / C::coarse_dt);
  int const n = mesh.getNumberOfNodes();
  vectorReal current("fraction current", n), previous("fraction previous", n);
  Kokkos::deep_copy(current, initial);
  solver::RhsAcoustic empty(arrayReal("no source", 0, 0), vectorInt("no element", 0), arrayReal("no weights", 0, 0));
  solver::SEMsolverDataAcoustic data(solver::WavefieldAcoustic(previous, current), empty);
  if (config.local) {
    static_cast<solver::SEMsolver2DLocalTimeStep<2, false>&>(*integrator).initializePreviousField(data);
  } else {
    integrator->computeForces(config.dt(), 0, data);
    auto mass = integrator->getMassMatrixAcoustic();
    auto force = integrator->getForceVector(0);
    float const dt = config.dt();
    Kokkos::parallel_for(
        "fraction initial previous", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int i) { previous(i) = current(i) - 0.5f * dt * dt * force(i) / mass(i); });
  }
  Kokkos::fence();
  double const setup = ms(start);
  start = Clock::now();
  for (int step = 0; step < config.steps; ++step) {
    integrator->computeOneStep(config.dt(), step, data);
    data.swapWavefields();
  }
  Kokkos::fence();
  double const loop = ms(start);
  return {setup, loop, solver::validation::layerHost(data.getCurrentField(0))};
}
double quantile(std::vector<double> v, double p) {
  std::sort(v.begin(), v.end());
  double const x = p * (v.size() - 1);
  size_t const a = static_cast<size_t>(x), b = std::min(a + 1, v.size() - 1);
  if (a == b || v[a] == v[b]) return v[a];
  return v[a] + (x - a) * (v[b] - v[a]);
}

void benchmark(const std::filesystem::path& output, int repeats) {
  std::filesystem::create_directories(output);
  std::ofstream raw(output / "fraction_samples.csv"), summary(output / "fraction_summary.csv"),
      scenarios(output / "fraction_scenarios.csv");
  for (auto* f : {&raw, &summary, &scenarios}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  std::string const common = "scenario,case,local,dt_ms,steps";
  raw << common << ",trial,setup_ms,loop_ms,relative_l2_error,finite\n";
  summary << common
          << ",setup_median_ms,loop_median_ms,loop_q25_ms,loop_q75_ms,total_median_ms,relative_l2_error,finite\n";
  scenarios << "scenario,elements,fine_elements,support_elements,nodes,fine_nodes,interface_nodes,reference_delta_l2\n";
  std::vector<Case> const cases{
      {"LTS_12ms", true, 150},      {"LTS_9ms", true, 200},       {"LTS_6ms", true, 300},
      {"LTS_3ms", true, 600},       {"Global_12ms", false, 150},  {"Global_10ms", false, 180},
      {"Global_9ms", false, 200},   {"Global_7.5ms", false, 240}, {"Global_6ms", false, 300},
      {"Global_4.5ms", false, 400}, {"Global_3ms", false, 600},   {"Global_1.5ms", false, 1200}};
  std::mt19937 rng(20261107);
  // -1 reproduces the original 1500/2250/3000 m/s three-layer baseline.
  for (int rows : {0, 1, 3, 6, 10, 20, 30, -1}) {
    std::string const name = rows < 0 ? "original_three_layers" : "fast_rows_" + std::to_string(rows);
    auto d = solver::validation::threeLayerMesh();
    if (rows >= 0) {
      for (int e = 0; e < d.n_element_; ++e) d.model_vp_element_.view_host()(e) = e / C::ex < rows ? 3000.f : 1500.f;
      d.model_vp_element_.modify_host();
    }
    model::ModelUnstruct<float, int> mesh(d);
    auto plan = solver::buildSEM2DLocalTimeStepPlan(mesh, C::cfl, C::coarse_dt);
    vectorReal initial("fraction Gaussian", mesh.getNumberOfNodes());
    auto h = Kokkos::create_mirror_view(initial);
    for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
      double const x = (mesh.nodeCoord(i, 0) - 3000) / 160.;
      double const z = (mesh.nodeCoord(i, 2) - 2500) / 160.;
      h(i) = static_cast<float>(std::exp(-0.5 * (x * x + z * z)));
    }
    Kokkos::deep_copy(initial, h);
    solver::SEMsolver2D<2, false> mass_solver;
    mass_solver.computeFEInit(mesh, {0, 0, 0}, false, 0);
    auto mass = solver::validation::layerHost(mass_solver.getMassMatrixAcoustic());
    auto reference = run({"reference_0.375ms", false, 4800}, mesh, initial).field;
    double norm = 0;
    for (size_t i = 0; i < mass.size(); ++i) norm += static_cast<double>(mass[i]) * reference[i] * reference[i];
    if (!(norm > 0) || !std::isfinite(norm)) throw std::runtime_error("Invalid reference norm");
    auto error = [&](const std::vector<float>& values) {
      double square = 0;
      for (size_t i = 0; i < mass.size(); ++i) {
        if (!std::isfinite(values[i])) return std::numeric_limits<double>::infinity();
        double const diff = static_cast<double>(values[i]) - reference[i];
        square += mass[i] * diff * diff;
      }
      return std::sqrt(square / norm);
    };
    double const reference_delta = error(run({"reference_0.75ms", false, 2400}, mesh, initial).field);
    scenarios << name << ',' << d.n_element_ << ',' << plan.fine_elements.size() << ','
              << plan.fine_support_elements.size() << ',' << d.n_node_ << ',' << plan.fine_nodes.size() << ','
              << plan.interface_nodes.size() << ',' << reference_delta << '\n';
    // Check that dt refinement is not silently changing the LTS zone.
    for (auto const& c : cases)
      if (c.local) {
        auto p = solver::buildSEM2DLocalTimeStepPlan(mesh, C::cfl * c.dt() / C::coarse_dt, c.dt());
        if (p.fine_elements != plan.fine_elements || p.fine_nodes != plan.fine_nodes)
          throw std::runtime_error("The LTS partition changed during the dt sweep");
      }
    for (int w = 0; w < 2; ++w)
      for (auto const& c : cases) run(c, mesh, initial);
    std::vector<std::vector<double>> setup(cases.size()), loop(cases.size()), total(cases.size()), errors(cases.size());
    std::vector<int> order(cases.size());
    std::iota(order.begin(), order.end(), 0);
    for (int trial = 0; trial < repeats; ++trial) {
      std::shuffle(order.begin(), order.end(), rng);
      for (int index : order) {
        auto const& c = cases[index];
        auto s = run(c, mesh, initial);
        double const e = error(s.field);
        raw << name << ',' << c.name << ',' << c.local << ',' << c.dt() * 1000 << ',' << c.steps << ',' << trial << ','
            << s.setup_ms << ',' << s.loop_ms << ',' << e << ',' << std::isfinite(e) << '\n';
        setup[index].push_back(s.setup_ms);
        loop[index].push_back(s.loop_ms);
        total[index].push_back(s.setup_ms + s.loop_ms);
        errors[index].push_back(e);
      }
    }
    for (size_t j = 0; j < cases.size(); ++j) {
      auto const& c = cases[j];
      summary << name << ',' << c.name << ',' << c.local << ',' << c.dt() * 1000 << ',' << c.steps << ','
              << quantile(setup[j], .5) << ',' << quantile(loop[j], .5) << ',' << quantile(loop[j], .25) << ','
              << quantile(loop[j], .75) << ',' << quantile(total[j], .5) << ',' << quantile(errors[j], .5) << ','
              << std::isfinite(errors[j].front()) << '\n';
    }
    raw.flush();
    summary.flush();
    scenarios.flush();
    std::cout << name << ": fine " << plan.fine_elements.size() << '/' << d.n_element_
              << "; reference refinement delta " << reference_delta << std::endl;
  }
  std::ofstream protocol(output / "fraction_protocol.json");
  protocol.exceptions(std::ios::failbit | std::ios::badbit);
  protocol << "{\n  \"backend\": \"" << Kokkos::DefaultExecutionSpace::name() << "\",\n"
           << "  \"concurrency\": " << Kokkos::DefaultExecutionSpace().concurrency() << ",\n"
           << "  \"repetitions\": " << repeats << ", \"warmups_per_case\": 2,\n"
           << "  \"seed\": 20261107, \"duration_s\": 1.8, \"order\": 2,\n"
           << "  \"reference_dt_ms\": 0.375, \"reference_check_dt_ms\": 0.75,\n"
           << "  \"lts_partition\": \"fixed per scenario; cfl scales with dt\",\n"
           << "  \"loop_excludes\": \"setup, diagnostics, file IO, host copies\"\n}\n";
}
// Refine the 2% budget boundary on integer step counts, then measure both selected
// configurations in interleaved order. This removes the coarse dt grid's threshold
// jumps. Bisection locates a crossing, not a proof of global error monotonicity.
void refineBudget(const std::filesystem::path& output, int repeats) {
  std::filesystem::create_directories(output);
  std::ofstream trace(output / "budget_search.csv"), raw(output / "budget_samples.csv"),
      summary(output / "budget_summary.csv");
  for (auto* f : {&trace, &raw, &summary}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  trace << "scenario,local,steps,dt_ms,relative_l2_error\n";
  raw << "scenario,local,steps,dt_ms,trial,setup_ms,loop_ms,relative_l2_error\n";
  summary << "scenario,local,steps,dt_ms,setup_median_ms,loop_median_ms,loop_q25_ms,loop_q75_ms,total_median_ms,"
             "relative_l2_error\n";
  std::mt19937 rng(20261108);
  for (int rows : {0, 1, 3, 6, 10, 20, 30, -1}) {
    std::string const name = rows < 0 ? "original_three_layers" : "fast_rows_" + std::to_string(rows);
    auto d = solver::validation::threeLayerMesh();
    if (rows >= 0) {
      for (int e = 0; e < d.n_element_; ++e) d.model_vp_element_.view_host()(e) = e / C::ex < rows ? 3000.f : 1500.f;
      d.model_vp_element_.modify_host();
    }
    model::ModelUnstruct<float, int> mesh(d);
    vectorReal initial("budget Gaussian", mesh.getNumberOfNodes());
    auto h = Kokkos::create_mirror_view(initial);
    for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
      double const x = (mesh.nodeCoord(i, 0) - 3000) / 160., z = (mesh.nodeCoord(i, 2) - 2500) / 160.;
      h(i) = static_cast<float>(std::exp(-0.5 * (x * x + z * z)));
    }
    Kokkos::deep_copy(initial, h);
    solver::SEMsolver2D<2, false> mass_solver;
    mass_solver.computeFEInit(mesh, {0, 0, 0}, false, 0);
    auto mass = solver::validation::layerHost(mass_solver.getMassMatrixAcoustic());
    auto reference = run({"reference", false, 4800}, mesh, initial).field;
    double norm = 0;
    for (size_t i = 0; i < mass.size(); ++i) norm += static_cast<double>(mass[i]) * reference[i] * reference[i];
    if (!(norm > 0) || !std::isfinite(norm)) throw std::runtime_error("Invalid budget reference");
    auto error = [&](const std::vector<float>& values) {
      double square = 0;
      for (size_t i = 0; i < mass.size(); ++i) {
        if (!std::isfinite(values[i])) return std::numeric_limits<double>::infinity();
        double const diff = static_cast<double>(values[i]) - reference[i];
        square += mass[i] * diff * diff;
      }
      return std::sqrt(square / norm);
    };
    std::vector<Case> selected;
    for (bool local : {false, true}) {
      std::map<int, double> cache;
      auto evaluate = [&](int steps) {
        if (!cache.count(steps)) {
          Case c{local ? "local" : "global", local, steps};
          double const e = error(run(c, mesh, initial).field);
          cache[steps] = e;
          trace << name << ',' << local << ',' << steps << ',' << c.dt() * 1000 << ',' << e << '\n';
        }
        return cache[steps];
      };
      int lo = 150, hi = 1200;
      if (evaluate(hi) > .02) throw std::runtime_error("No passing step count in the budget bracket");
      if (evaluate(lo) <= .02)
        hi = lo;
      else {
        while (hi - lo > 1) {
          int mid = (lo + hi) / 2;
          if (evaluate(mid) <= .02)
            hi = mid;
          else
            lo = mid;
        }
        if (!(evaluate(hi - 1) > .02)) throw std::runtime_error("Invalid budget boundary");
      }
      if (!(evaluate(hi) <= .02)) throw std::runtime_error("Budget candidate fails");
      if (local) {
        auto fixed = solver::buildSEM2DLocalTimeStepPlan(mesh, C::cfl, C::coarse_dt);
        float const dt = 1.8f / hi;
        auto refined = solver::buildSEM2DLocalTimeStepPlan(mesh, C::cfl * dt / C::coarse_dt, dt);
        if (fixed.fine_nodes != refined.fine_nodes || fixed.fine_elements != refined.fine_elements)
          throw std::runtime_error("Budget refinement changed the partition");
      }
      selected.push_back({local ? "local" : "global", local, hi});
    }
    for (int w = 0; w < 2; ++w)
      for (auto const& c : selected) run(c, mesh, initial);
    std::vector<std::vector<double>> setups(2), loops(2), totals(2), errors(2);
    std::vector<int> order{0, 1};
    for (int trial = 0; trial < repeats; ++trial) {
      std::shuffle(order.begin(), order.end(), rng);
      for (int j : order) {
        auto const& c = selected[j];
        auto s = run(c, mesh, initial);
        double const e = error(s.field);
        if (!(e <= .02)) throw std::runtime_error("Measured candidate exceeds the budget");
        raw << name << ',' << c.local << ',' << c.steps << ',' << c.dt() * 1000 << ',' << trial << ',' << s.setup_ms
            << ',' << s.loop_ms << ',' << e << '\n';
        setups[j].push_back(s.setup_ms);
        loops[j].push_back(s.loop_ms);
        totals[j].push_back(s.setup_ms + s.loop_ms);
        errors[j].push_back(e);
      }
    }
    for (int j = 0; j < 2; ++j) {
      auto const& c = selected[j];
      summary << name << ',' << c.local << ',' << c.steps << ',' << c.dt() * 1000 << ',' << quantile(setups[j], .5)
              << ',' << quantile(loops[j], .5) << ',' << quantile(loops[j], .25) << ',' << quantile(loops[j], .75)
              << ',' << quantile(totals[j], .5) << ',' << quantile(errors[j], .5) << '\n';
    }
    trace.flush();
    raw.flush();
    summary.flush();
    std::cout << name << ": global " << selected[0].steps << " steps, local " << selected[1].steps
              << "; global/local = " << quantile(loops[0], .5) / quantile(loops[1], .5) << std::endl;
  }
  std::ofstream protocol(output / "budget_protocol.json");
  protocol.exceptions(std::ios::failbit | std::ios::badbit);
  protocol << "{\n  \"budget\": 0.02, \"min_steps\": 150, \"max_steps\": 1200,\n"
           << "  \"search\": \"integer bisection; adjacent failing count checked unless min_steps passes\",\n"
           << "  \"reference_dt_ms\": 0.375, \"duration_s\": 1.8,\n"
           << "  \"repetitions\": " << repeats << ", \"warmups_per_case\": 2, \"seed\": 20261108\n}\n";
}
}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int status = 0;
  try {
    auto output = argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("fraction-results");
    int const repeats = argc > 2 ? std::stoi(argv[2]) : 7;
    if (repeats < 3) throw std::invalid_argument("Use at least three repetitions");
    if (argc > 3 && std::string(argv[3]) == "refine")
      refineBudget(output, repeats);
    else
      benchmark(output, repeats);
  } catch (std::exception const& e) {
    std::cerr << e.what() << '\n';
    status = 1;
  }
  Kokkos::finalize();
  return status;
}
