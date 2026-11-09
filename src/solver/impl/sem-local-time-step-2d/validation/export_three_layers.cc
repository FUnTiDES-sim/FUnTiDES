#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "three_layer_case.h"

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int code = 0;
  try {
    using C = solver::validation::ThreeLayerResult;
    auto const result = solver::validation::runThreeLayerValidation();
    auto out = argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("validation-three-layers");
    std::filesystem::create_directories(out);
    auto open = [&](const char* name) {
      std::ofstream f(out / name);
      f.exceptions(std::ios::failbit | std::ios::badbit);
      f << std::setprecision(12);
      return f;
    };
    auto elements = open("elements.csv");
    elements << "element,x,z,vp,estimated_dt,selected_dt,fine\n";
    std::vector<bool> fine(C::ex * C::ez, false);
    for (int e : result.plan.fine_elements) fine[e] = true;
    for (int e = 0; e < C::ex * C::ez; ++e)
      elements << e << ',' << ((e % C::ex) + 0.5) * C::h << ',' << ((e / C::ex) + 0.5) * C::h << ',' << result.vp[e]
               << ',' << result.plan.element_steps[e].estimate << ','
               << (fine[e] ? result.plan.fine_dt : result.plan.coarse_dt) << ',' << fine[e] << '\n';
    auto fields = open("fields.csv");
    fields << "frame,time,node,x,z,local,reference\n";
    for (size_t frame = 0; frame < result.snapshots.size(); ++frame) {
      auto const& s = result.snapshots[frame];
      for (size_t i = 0; i < s.local.size(); ++i)
        fields << frame << ',' << s.time << ',' << i << ',' << result.x[i] << ',' << result.z[i] << ',' << s.local[i]
               << ',' << s.reference[i] << '\n';
    }
    auto history = open("history.csv");
    history << "time,energy,relative_energy_drift,relative_l2_error\n";
    double const initial_energy = result.history.front().energy;
    double max_drift = 0;
    for (auto const& s : result.history) {
      double const drift = s.energy / initial_energy - 1;
      max_drift = std::max(max_drift, std::abs(drift));
      history << s.time << ',' << s.energy << ',' << drift << ',' << s.relative_error << '\n';
    }
    auto metadata = open("metadata.json");
    metadata << "{\n  \"nx\": " << C::nx << ", \"nz\": " << C::nz << ",\n"
             << "  \"coarse_dt\": " << C::coarse_dt << ", \"fine_dt\": " << result.plan.fine_dt << ",\n"
             << "  \"reference_dt\": " << C::reference_dt << ", \"cfl\": " << C::cfl << ",\n"
             << "  \"relative_l2_error\": " << result.history.back().relative_error << ",\n"
             << "  \"max_relative_energy_drift\": " << max_drift << "\n}\n";
    std::cout << "Wrote three-layer validation data to " << out
              << "\nFinal relative L2 error: " << result.history.back().relative_error
              << "\nMax relative energy drift: " << max_drift << '\n';
  } catch (std::exception const& e) {
    std::cerr << e.what() << '\n';
    code = 1;
  }
  Kokkos::finalize();
  return code;
}
