#ifndef FUNTIDES_SEM2D_BENCHMARK_COMMON_H_
#define FUNTIDES_SEM2D_BENCHMARK_COMMON_H_
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "sem_solver_2d_local_time_step.h"
namespace solver::validation::benchmark {
using LTS = solver::SEMsolver2DLocalTimeStep<2, false>;
using Global = solver::SEMsolver2D<2, false>;
using Clock = std::chrono::steady_clock;
inline std::vector<float> hostField(const vectorReal& field) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field);
  return {h.data(), h.data() + h.extent(0)};
}

struct Parameters {
  float x, z, sigma, duration, macro_dt, cfl;
  bool absorbing{false};
  bool double_reference{false};
  bool free_top{false};
  float receiver_depth{0};
};
struct Case {
  std::string name;
  bool local;
  int steps;
};
struct Sample {
  double setup_ms, loop_ms;
  std::vector<float> field;
};
inline double milliseconds(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
inline double quantile(std::vector<double> v, double p) {
  std::sort(v.begin(), v.end());
  double x = p * (v.size() - 1);
  size_t a = static_cast<size_t>(x), b = std::min(a + 1, v.size() - 1);
  return a == b || v[a] == v[b] ? v[a] : v[a] + (x - a) * (v[b] - v[a]);
}

// All preparation is outside the integration timing. The generic target steps in
// model.csv are intentionally not consumed by this fixed-ratio integrator.
struct Run {
  using Precise = Kokkos::View<double*, vectorReal::array_layout, vectorReal::device_type>;
  std::unique_ptr<solver::Solver> integrator;
  solver::SEMsolverDataAcoustic data;
  Precise precise_current, precise_previous;
  vectorInt surface_mask;
  float dt;
  Run(model::ModelUnstruct<float, int>& mesh, const vectorReal& initial, bool local, float step, const Parameters& p,
      bool double_update = false)
      : data(solver::WavefieldAcoustic(vectorReal("previous", mesh.getNumberOfNodes()),
                                       vectorReal("current", mesh.getNumberOfNodes())),
             solver::RhsAcoustic(arrayReal("no source", 0, 0), vectorInt("no element", 0),
                                 arrayReal("no weights", 0, 0))),
        dt(step) {
    if (local)
      integrator = std::make_unique<LTS>();
    else
      integrator = std::make_unique<Global>();
    integrator->computeFEInit(mesh, {0, 0, 0}, false, 0);
    Kokkos::deep_copy(data.getCurrentField(0), initial);
    if (local) {
      auto& lts = static_cast<LTS&>(*integrator);
      lts.configureLocalTimeStep(dt, p.cfl * dt / p.macro_dt);
      lts.initializePreviousField(data);
    } else {
      integrator->computeForces(dt, 0, data);
      auto cur = data.getCurrentField(0), prev = data.getPreviousField(0);
      auto mass = integrator->getMassMatrixAcoustic(), force = integrator->getForceVector(0);
      Kokkos::parallel_for(
          "inclusion initial previous", Kokkos::RangePolicy<>(0, mesh.getNumberOfNodes()),
          KOKKOS_LAMBDA(int i) { prev(i) = cur(i) - 0.5f * step * step * force(i) / mass(i); });
      if (double_update) {
        int n = mesh.getNumberOfNodes();
        precise_current = Precise("reference current f64", n);
        precise_previous = Precise("reference previous f64", n);
        surface_mask = vectorInt("reference Dirichlet mask", n);
        auto mask = surface_mask;
        auto geometry = mesh;
        auto c = precise_current, old = precise_previous;
        double h = step;
        Kokkos::parallel_for(
            "reference Taylor f64", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(int i) {
              c(i) = cur(i);
              old(i) = c(i) - .5 * h * h * static_cast<double>(force(i)) / mass(i);
              mask(i) = geometry.isFreeSurface(i);
              if (mask(i)) c(i) = old(i) = cur(i) = 0;
              prev(i) = old(i);
            });
      }
    }
    Kokkos::fence();
  }
  void advance(int step) {
    if (!precise_current.extent(0))
      integrator->computeOneStep(dt, step, data);
    else {
      // Spatial forces remain float32; only the cancellation-prone temporal
      // recurrence is accumulated in float64, exclusively for reference runs.
      integrator->computeForces(dt, step, data);
      auto c = precise_current, old = precise_previous;
      auto mask = surface_mask;
      auto mass = integrator->getMassMatrixAcoustic(), damping = integrator->getDampingMatrix(0);
      auto force = integrator->getForceVector(0), next = data.getPreviousField(0);
      double h = dt;
      Kokkos::parallel_for(
          "reference centered update f64", Kokkos::RangePolicy<>(0, c.extent(0)), KOKKOS_LAMBDA(int i) {
            double m = mass(i), d = damping(i);
            if (mask(i)) {
              c(i) = old(i) = next(i) = 0;
              return;
            }
            old(i) = (2 * m * c(i) - (m - .5 * h * d) * old(i) - h * h * force(i)) / (m + .5 * h * d);
            next(i) = old(i);
          });
      Kokkos::fence();
      std::swap(precise_current, precise_previous);
    }
    data.swapWavefields();
  }
};

inline Sample run(const Case& c, model::ModelUnstruct<float, int>& mesh, const vectorReal& initial, const Parameters& p,
                  bool reference = false) {
  Kokkos::fence();
  auto start = Clock::now();
  Run simulation(mesh, initial, c.local, p.duration / c.steps, p, reference && p.double_reference);
  double setup = milliseconds(start);
  start = Clock::now();
  for (int step = 0; step < c.steps; ++step) simulation.advance(step);
  Kokkos::fence();
  double loop = milliseconds(start);
  return {setup, loop, hostField(simulation.data.getCurrentField(0))};
}

// Power iteration on the symmetric mass-normalized operator A=M^-1/2 K M^-1/2.
// Records the residual; this is a numerical spectral estimate, not a certified bound.
inline double spectralEstimate(model::ModelUnstruct<float, int>& mesh, const std::filesystem::path& output) {
  LTS op;
  op.computeFEInit(mesh, {0, 0, 0}, false, 0);
  auto mass = op.getMassMatrixAcoustic();
  int const n = mesh.getNumberOfNodes();
  vectorReal x("spectral x", n), u("spectral u", n), force("spectral K u", n), ax("spectral A x", n);
  auto h = Kokkos::create_mirror_view(x);
  std::mt19937 rng(4711);
  std::uniform_real_distribution<float> uniform(-1, 1);
  double square = 0;
  for (int i = 0; i < n; ++i) {
    h(i) = uniform(rng);
    square += h(i) * h(i);
  }
  for (int i = 0; i < n; ++i) h(i) /= std::sqrt(square);
  Kokkos::deep_copy(x, h);
  double lambda = 0, residual = 0;
  int iteration = 0;
  for (; iteration < 2500; ++iteration) {
    Kokkos::parallel_for(
        "spectral pressure", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int i) { u(i) = x(i) / Kokkos::sqrt(mass(i)); });
    op.applyStiffness(u, force);
    Kokkos::parallel_for(
        "spectral mass normalize", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int i) { ax(i) = force(i) / Kokkos::sqrt(mass(i)); });
    Kokkos::parallel_reduce(
        "spectral rayleigh", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int i, double& sum) { sum += static_cast<double>(x(i)) * ax(i); }, lambda);
    double norm = 0, remainder = 0;
    Kokkos::parallel_reduce(
        "spectral norm", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int i, double& sum) { sum += static_cast<double>(ax(i)) * ax(i); }, norm);
    Kokkos::parallel_reduce(
        "spectral residual", Kokkos::RangePolicy<>(0, n),
        KOKKOS_LAMBDA(int i, double& sum) {
          double r = ax(i) - lambda * x(i);
          sum += r * r;
        },
        remainder);
    residual = std::sqrt(remainder) / lambda;
    if (residual < 1e-5) break;
    double inverse = 1 / std::sqrt(norm);
    Kokkos::parallel_for(
        "spectral normalize", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(int i) { x(i) = ax(i) * inverse; });
  }
  if (!(lambda > 0) || !std::isfinite(lambda)) throw std::runtime_error("Invalid spectral estimate");
  double limit = 2 / std::sqrt(lambda);
  std::ofstream out(output / "spectral.json");
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << std::setprecision(12) << "{\"lambda_max_estimate\": " << lambda
      << ", \"global_dt_limit_estimate_s\": " << limit << ", \"relative_residual\": " << residual
      << ", \"iterations\": " << std::min(iteration + 1, 2500)
      << ", \"iteration_limit_reached\": " << (iteration == 2500 ? "true" : "false") << "}\n";
  return limit;
}

}  // namespace solver::validation::benchmark
#endif
