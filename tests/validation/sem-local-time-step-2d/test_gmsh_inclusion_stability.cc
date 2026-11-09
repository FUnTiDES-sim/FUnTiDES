#include <iostream>
#include <map>
#include <sstream>

#include "gmsh_2d_reader.h"
#include "quad_observation.h"
#include "sem2d_benchmark_common.h"

namespace {
using namespace solver::validation::benchmark;
using Mesh = model::ModelUnstruct<float, int>;

Mesh read(const std::filesystem::path& root) {
  std::ifstream in(root / "materials.csv");
  std::string line;
  std::getline(in, line);
  std::map<int, model::GmshAcousticMaterial> materials;
  while (std::getline(in, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int tag;
    float vp, rho;
    if (!(row >> tag >> vp >> rho)) throw std::runtime_error("Invalid materials.csv");
    materials.emplace(tag, model::GmshAcousticMaterial{vp, rho});
  }
  return Mesh(model::readGmsh2D(root / "bilayer.msh", 2, materials,
                              {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}}));
}

vectorReal initial(Mesh& mesh) {
  vectorReal result("common physical initial field", mesh.getNumberOfNodes());
  auto h = Kokkos::create_mirror_view(result);
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double x = mesh.nodeCoord(i, 0), z = mesh.nodeCoord(i, 2);
    double gaussian = std::exp(-.5 * (std::pow((x - 2200) / 160, 2) + std::pow((z - 2800) / 160, 2))) -
                      std::exp(-.5 * (std::pow((x - 2200) / 160, 2) + std::pow((z - 3200) / 160, 2)));
    double noise = 0;
    if (x > 2990 && x < 3010 && z > 890 && z < 910) {
      double pi = std::acos(-1.);
      noise = 1e-6 * std::pow(std::sin(pi * (x - 2990) / 20), 2) *
              std::pow(std::sin(pi * (z - 890) / 20), 2) * std::cos(2 * pi * (x - 3000)) *
              std::cos(2 * pi * (z - 900));
    }
    h(i) = mesh.isFreeSurface(i) ? 0 : gaussian + noise;
  }
  Kokkos::deep_copy(result, h);
  return result;
}

struct Peak { double value = 0; int node = 0; bool finite = true; };
Peak peak(const std::vector<float>& field) {
  Peak p;
  for (size_t i = 0; i < field.size(); ++i) {
    if (!std::isfinite(field[i])) p.finite = false;
    if (std::abs(field[i]) > p.value) p.value = std::abs(field[i]), p.node = i;
  }
  return p;
}
double relative(const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& mass) {
  double error = 0, norm = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    double d = static_cast<double>(a[i]) - b[i];
    error += mass[i] * d * d;
    norm += static_cast<double>(mass[i]) * b[i] * b[i];
  }
  if (!(norm > 0) || !std::isfinite(error)) throw std::runtime_error("Invalid error norm");
  return std::sqrt(error / norm);
}
bool inInclusion(Mesh& mesh, int node) {
  double x = mesh.nodeCoord(node, 0), z = mesh.nodeCoord(node, 2);
  return x >= 2800 && x <= 3200 && z >= 700 && z <= 1100;
}

bool execute(const std::filesystem::path& baseline, const std::filesystem::path& root, double factor) {
  auto coarse = read(baseline), mesh = read(root);
  auto initial_coarse = initial(coarse), initial_refined = initial(mesh);
  std::filesystem::create_directories(root / "coarse_spectrum");
  std::filesystem::create_directories(root / "refined_spectrum");
  double coarse_limit = spectralEstimate(coarse, root / "coarse_spectrum");
  double refined_limit = spectralEstimate(mesh, root / "refined_spectrum");
  Parameters p{};
  p.duration = 3.6f;
  p.double_reference = p.free_top = p.absorbing = true;
  int steps = 6 * static_cast<int>(std::ceil(p.duration / (factor * refined_limit) / 6));
  p.macro_dt = p.duration / steps;
  if (!(p.macro_dt > refined_limit && p.macro_dt < .9 * coarse_limit))
    throw std::runtime_error("No separated instability window: refine or recalibrate the patch");
  auto probe = solver::buildSEM2DLocalTimeStepPlan(mesh, 1.f, 1e-9f);
  double min_ratio = 1e30, outside_ratio = 1e30;
  for (int e = 0; e < mesh.getNumberOfElements(); ++e) {
    min_ratio = std::min(min_ratio, static_cast<double>(probe.element_steps[e].estimate));
    if (mesh.getModelVpOnElement(e) != 600)
      outside_ratio = std::min(outside_ratio, static_cast<double>(probe.element_steps[e].estimate));
  }
  p.cfl = 1.02 * std::max(.5 * p.macro_dt / min_ratio, p.macro_dt / outside_ratio);
  auto plan = solver::buildSEM2DLocalTimeStepPlan(mesh, p.cfl, p.macro_dt);
  if (plan.fine_elements.empty()) throw std::runtime_error("No fine elements in inclusion");
  for (int e : plan.fine_elements)
    if (mesh.getModelVpOnElement(e) != 600) throw std::runtime_error("Fine element outside inclusion");
  std::cout << "Calibrated macro " << p.macro_dt * 1000 << " ms, fine " << p.macro_dt * 500
            << " ms; global limits refined/control " << refined_limit * 1000 << '/' << coarse_limit * 1000
            << " ms; " << plan.fine_elements.size() << " fine elements, all at 600 m/s" << std::endl;
  {
    std::ofstream calibration(root / "calibration.json");
    calibration << std::setprecision(12) << "{\"macro_dt_s\":" << p.macro_dt << ",\"fine_dt_s\":"
                << p.macro_dt / 2 << ",\"coarse_global_limit_estimate_s\":" << coarse_limit
                << ",\"refined_global_limit_estimate_s\":" << refined_limit << ",\"cfl_partition\":" << p.cfl
                << ",\"fine_elements\":" << plan.fine_elements.size() << ",\"fine_outside_inclusion\":0}\n";
  }
  Run local(mesh, initial_refined, true, p.macro_dt, p), fine(mesh, initial_refined, false, p.macro_dt / 2, p),
      reference(mesh, initial_refined, false, p.macro_dt / 8, p, true),
      control(coarse, initial_coarse, false, p.macro_dt, p);
  auto unsafe = std::make_unique<Run>(mesh, initial_refined, false, p.macro_dt, p);
  auto mass = hostField(local.integrator->getMassMatrixAcoustic());
  auto damping = hostField(local.integrator->getDampingMatrix(0));
  double e0 = static_cast<LTS&>(*local.integrator).discreteEnergy(local.data);
  if (!(e0 > 0)) throw std::runtime_error("Nonpositive initial LTS energy");
  double initial_amplitude = peak(hostField(initial_refined)).value;
  double control_initial_amplitude = peak(hostField(initial_coarse)).value;
  std::vector<model::QuadObservation> observations;
  for (int j = 0; j < 51; ++j) observations.push_back(model::locateQuadObservation(mesh, 500 + 100 * j, 2950));
  std::ofstream history(root / "stability_history.csv"), receivers(root / "receivers.csv"),
      snapshots(root / "snapshots.csv");
  for (auto* f : {&history, &receivers, &snapshots}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  history << "time,control_peak,global_macro_peak,global_fine_peak,local_peak,reference_peak,energy,dissipation,balance_error,local_l2_error\n";
  receivers << "time,receiver,x,z,local,global_fine,reference\n";
  snapshots << "time,node,x,z,local,reference\n";
  double lost = 0, max_balance = 0, max_error = 0, local_max = 0, fine_max = 0, control_max = 0;
  double receiver_error = 0, receiver_norm = 0, first_growth_time = -1, divergence_time = -1;
  double first_growth_x = 0, first_growth_z = 0, divergence_amplitude = 0;
  bool first_growth_in_inclusion = false;
  std::vector<float> a, b, f;
  for (int step = 0; step <= steps; ++step) {
    if (step) {
      auto old = hostField(local.data.getPreviousField(0));
      local.advance(step - 1);
      control.advance(step - 1);
      for (int sub = 0; sub < 2; ++sub) fine.advance(2 * (step - 1) + sub);
      for (int sub = 0; sub < 8; ++sub) reference.advance(8 * (step - 1) + sub);
      if (unsafe) unsafe->advance(step - 1);
      a = hostField(local.data.getCurrentField(0));
      for (size_t i = 0; i < a.size(); ++i) {
        double v = (static_cast<double>(a[i]) - old[i]) / (2 * p.macro_dt);
        lost += p.macro_dt * damping[i] * v * v;
      }
    } else a = hostField(local.data.getCurrentField(0));
    b = hostField(reference.data.getCurrentField(0));
    f = hostField(fine.data.getCurrentField(0));
    auto c = hostField(control.data.getCurrentField(0));
    auto pa = peak(a), pb = peak(b), pf = peak(f), pc = peak(c);
    if (!(pa.finite && pb.finite && pf.finite && pc.finite)) throw std::runtime_error("A required stable run diverged");
    local_max = std::max(local_max, pa.value / initial_amplitude);
    fine_max = std::max(fine_max, pf.value / initial_amplitude);
    control_max = std::max(control_max, pc.value / control_initial_amplitude);
    if (std::max({local_max, fine_max, control_max}) >= 10) throw std::runtime_error("Required stable run grew tenfold");
    for (int i = 0; i < mesh.getNumberOfNodes(); ++i)
      if (mesh.isFreeSurface(i) && (a[i] != 0 || b[i] != 0 || f[i] != 0))
        throw std::runtime_error("Nonzero free-surface pressure");
    double unsafe_peak = -1;
    if (unsafe) {
      auto g = hostField(unsafe->data.getCurrentField(0));
      auto pg = peak(g);
      if (!pg.finite) throw std::runtime_error("Unsafe witness reached NaN before recording divergence");
      unsafe_peak = pg.value;
      if (first_growth_time < 0 && pg.value > 10 * initial_amplitude) {
        first_growth_time = step * static_cast<double>(p.macro_dt);
        first_growth_x = mesh.nodeCoord(pg.node, 0);
        first_growth_z = mesh.nodeCoord(pg.node, 2);
        first_growth_in_inclusion = inInclusion(mesh, pg.node);
      }
      if (pg.value > 1e5 * initial_amplitude) {
        divergence_time = step * static_cast<double>(p.macro_dt);
        divergence_amplitude = pg.value / initial_amplitude;
        unsafe.reset();
        std::cout << "Global macro divergence at " << divergence_time << " s; first growth at ("
                  << first_growth_x << ',' << first_growth_z << ")" << std::endl;
      }
    }
    double energy = static_cast<LTS&>(*local.integrator).discreteEnergy(local.data);
    double balance = std::abs(energy + lost - e0) / e0;
    if (!std::isfinite(energy) || energy < 0 || balance > 1e-4) throw std::runtime_error("LTS energy balance failed");
    double error = relative(a, b, mass);
    max_balance = std::max(max_balance, balance);
    max_error = std::max(max_error, error);
    history << step * p.macro_dt << ',' << pc.value << ',';
    if (unsafe_peak >= 0) history << unsafe_peak;
    history << ',' << pf.value << ',' << pa.value << ',' << pb.value << ',' << energy << ',' << lost << ','
            << balance << ',' << error << '\n';
    for (size_t j = 0; j < observations.size(); ++j) {
      auto const& o = observations[j];
      double va = model::observeQuad(o, a), vb = model::observeQuad(o, b);
      receiver_error += (va - vb) * (va - vb);
      receiver_norm += vb * vb;
      receivers << step * p.macro_dt << ',' << j << ',' << o.x << ',' << o.z << ',' << va << ','
                << model::observeQuad(o, f) << ',' << vb << '\n';
    }
    if (step == 0 || step == steps / 6 || step == steps / 3 || step == steps / 2 || step == 2 * steps / 3 || step == steps)
      for (size_t i = 0; i < a.size(); ++i)
        snapshots << step * p.macro_dt << ',' << i << ',' << mesh.nodeCoord(i, 0) << ',' << mesh.nodeCoord(i, 2)
                  << ',' << a[i] << ',' << b[i] << '\n';
    if (step && step % (steps / 6) == 0)
      std::cout << "Physical time " << step * p.macro_dt << '/' << p.duration << " s" << std::endl;
  }
  double local_error = relative(a, b, mass), fine_error = relative(f, b, mass);
  auto reference_check = run({"Reference check", false, 4 * steps}, mesh, initial_refined, p, true).field;
  double reference_delta = relative(reference_check, b, mass);
  double traces_error = std::sqrt(receiver_error / receiver_norm);
  bool passed = divergence_time >= 0 && first_growth_in_inclusion && local_error < .01 && fine_error < .01 &&
                traces_error < .01 && reference_delta < .001;
  std::ofstream result(root / "test_result.json");
  result << std::boolalpha << std::setprecision(12)
         << "{\"passed\":" << passed << ",\"physical_duration_s\":" << p.duration
         << ",\"macro_steps\":" << steps << ",\"macro_dt_s\":" << p.macro_dt
         << ",\"fine_dt_s\":" << p.macro_dt / 2 << ",\"reference_dt_s\":" << p.macro_dt / 8
         << ",\"fine_elements\":" << plan.fine_elements.size() << ",\"fine_outside_inclusion\":0"
         << ",\"first_tenfold_growth_time_s\":" << first_growth_time
         << ",\"first_growth_x_m\":" << first_growth_x << ",\"first_growth_z_up_m\":" << first_growth_z
         << ",\"first_growth_in_inclusion\":" << first_growth_in_inclusion
         << ",\"global_macro_divergence_time_s\":" << divergence_time
         << ",\"global_macro_amplification_at_stop\":" << divergence_amplitude
         << ",\"control_max_amplification\":" << control_max << ",\"global_fine_max_amplification\":" << fine_max
         << ",\"local_max_amplification\":" << local_max << ",\"local_final_mass_l2_error\":" << local_error
         << ",\"global_fine_final_mass_l2_error\":" << fine_error
         << ",\"local_receiver_relative_l2_error\":" << traces_error << ",\"reference_delta_l2\":" << reference_delta
         << ",\"max_local_energy_balance_error\":" << max_balance << ",\"max_local_trajectory_mass_l2_error\":" << max_error << "}\n";
  std::cout << (passed ? "PASS" : "FAIL") << ": local error " << local_error << ", fine error " << fine_error
            << ", trace error " << traces_error << ", reference delta " << reference_delta << std::endl;
  return passed;
}
}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int status = 1;
  try {
    if (argc < 3) throw std::invalid_argument("Usage: sem2d_gmsh_inclusion_stability baseline_dir refined_dir [macro_factor=1.2]");
    double factor = argc > 3 ? std::stod(argv[3]) : 1.2;
    if (!(factor > 1 && factor < 2)) throw std::invalid_argument("Macro factor must be between 1 and 2");
    status = execute(argv[1], argv[2], factor) ? 0 : 1;
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
  }
  Kokkos::finalize();
  return status;
}
