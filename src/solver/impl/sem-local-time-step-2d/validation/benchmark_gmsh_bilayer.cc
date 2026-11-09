#include <iostream>
#include <numeric>
#include <sstream>

#include "gmsh_2d_reader.h"
#include "quad_observation.h"
#include "sem2d_benchmark_common.h"

namespace {
using namespace solver::validation::benchmark;

std::map<int, model::GmshAcousticMaterial> readMaterials(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  std::map<int, model::GmshAcousticMaterial> materials;
  while (std::getline(in, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int tag;
    float vp, rho;
    if (!(row >> tag >> vp >> rho) || !materials.emplace(tag, model::GmshAcousticMaterial{vp, rho}).second)
      throw std::runtime_error("Invalid material table");
  }
  if (materials.empty()) throw std::runtime_error("Empty material table");
  return materials;
}

double relativeError(const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& mass) {
  double error = 0, norm = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (!std::isfinite(a[i]) || !std::isfinite(b[i])) return std::numeric_limits<double>::infinity();
    double d = static_cast<double>(a[i]) - b[i];
    error += mass[i] * d * d;
    norm += static_cast<double>(mass[i]) * b[i] * b[i];
  }
  if (!(norm > 0)) throw std::runtime_error("Invalid reference norm");
  return std::sqrt(error / norm);
}

void exportMesh(model::ModelUnstruct<float, int>& mesh, const solver::SEM2DLocalTimeStepPlan& plan,
                const std::filesystem::path& out) {
  std::ofstream nodes(out / "nodes.csv"), elements(out / "elements.csv"), triangles(out / "triangles.csv");
  for (auto* f : {&nodes, &elements, &triangles}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  nodes << "node,x,z,boundary\n";
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i)
    nodes << i << ',' << mesh.nodeCoord(i, 0) << ',' << mesh.nodeCoord(i, 2) << ','
          << static_cast<int>(mesh.boundaryType(i)) << '\n';
  elements << "element,n0,n1,n2,n3,vp,rho,min_spacing,dt_estimate,fine,support\n";
  triangles << "n0,n1,n2\n";
  std::vector<bool> fine(mesh.getNumberOfElements()), support(fine.size());
  for (int e : plan.fine_elements) fine[e] = true;
  for (int e : plan.fine_support_elements) support[e] = true;
  int q = mesh.getOrder();
  for (int e = 0; e < mesh.getNumberOfElements(); ++e) {
    elements << e << ',' << mesh.globalNodeIndex(e, 0, 0, 0) << ',' << mesh.globalNodeIndex(e, q, 0, 0) << ','
             << mesh.globalNodeIndex(e, q, 0, q) << ',' << mesh.globalNodeIndex(e, 0, 0, q) << ','
             << mesh.getModelVpOnElement(e) << ',' << mesh.getModelRhoOnElement(e) << ','
             << plan.element_steps[e].min_node_spacing << ',' << plan.element_steps[e].estimate << ',' << fine[e] << ','
             << support[e] << '\n';
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < q; ++i) {
        int a = mesh.globalNodeIndex(e, i, 0, j), b = mesh.globalNodeIndex(e, i + 1, 0, j),
            c = mesh.globalNodeIndex(e, i + 1, 0, j + 1), d = mesh.globalNodeIndex(e, i, 0, j + 1);
        triangles << a << ',' << b << ',' << c << '\n' << a << ',' << c << ',' << d << '\n';
      }
  }
}

void diagnostics(model::ModelUnstruct<float, int>& mesh, const vectorReal& initial, const Parameters& p, int steps,
                 int refinement, const std::vector<float>& mass, const std::filesystem::path& out) {
  Run local(mesh, initial, true, p.macro_dt, p), reference(mesh, initial, false, p.macro_dt / refinement, p, true);
  std::vector<model::QuadObservation> observations;
  double top = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) top = std::max(top, static_cast<double>(mesh.nodeCoord(i, 2)));
  for (int j = 0; j < 51; ++j)
    observations.push_back(model::locateQuadObservation(mesh, 500 + 100 * j, top - p.receiver_depth));
  std::ofstream snaps(out / "snapshots.csv"), traces(out / "receivers.csv"), history(out / "history.csv");
  for (auto* f : {&snaps, &traces, &history}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  snaps << "time,node,local,reference\n";
  traces << "time,receiver,x,z,local,reference\n";
  history << "time,energy,cumulative_dissipation,balance_error,relative_l2_error\n";
  auto damping = hostField(local.integrator->getDampingMatrix(0));
  double initial_energy = static_cast<LTS&>(*local.integrator).discreteEnergy(local.data), lost = 0, max_balance = 0,
         max_error = 0;
  if (!(initial_energy > 0)) throw std::runtime_error("Nonpositive LTS initial energy");
  for (int step = 0; step <= steps; ++step) {
    auto old = hostField(local.data.getPreviousField(0));
    if (step) {
      local.advance(step - 1);
      for (int sub = 0; sub < refinement; ++sub) reference.advance((step - 1) * refinement + sub);
    }
    auto a = hostField(local.data.getCurrentField(0)), b = hostField(reference.data.getCurrentField(0));
    for (size_t i = 0; i < a.size(); ++i) {
      if (!std::isfinite(a[i]) || !std::isfinite(b[i])) throw std::runtime_error("Divergent Gmsh diagnostic");
      if (mesh.isFreeSurface(i) && (a[i] != 0 || b[i] != 0))
        throw std::runtime_error("Nonzero pressure on free surface");
      if (step) {
        double v = (static_cast<double>(a[i]) - old[i]) / (2 * p.macro_dt);
        lost += p.macro_dt * damping[i] * v * v;
      }
    }
    double energy = static_cast<LTS&>(*local.integrator).discreteEnergy(local.data);
    double balance = std::abs(energy + lost - initial_energy) / initial_energy;
    max_balance = std::max(max_balance, balance);
    double error = relativeError(a, b, mass);
    max_error = std::max(max_error, error);
    if (!std::isfinite(energy) || energy < -1e-6 * initial_energy || balance > 1e-4)
      throw std::runtime_error("LTS energy balance failed");
    history << step * p.macro_dt << ',' << energy << ',' << lost << ',' << balance << ',' << error << '\n';
    for (size_t j = 0; j < observations.size(); ++j) {
      auto const& o = observations[j];
      traces << step * p.macro_dt << ',' << j << ',' << o.x << ',' << o.z << ',' << model::observeQuad(o, a) << ','
             << model::observeQuad(o, b) << '\n';
    }
    if (step == 0 || step == steps / 6 || step == steps / 3 || step == steps / 2 || step == 2 * steps / 3 ||
        step == steps)
      for (size_t i = 0; i < a.size(); ++i)
        snaps << step * p.macro_dt << ',' << i << ',' << a[i] << ',' << b[i] << '\n';
  }
  std::ofstream json(out / "diagnostics.json");
  json << std::setprecision(12) << "{\"max_relative_energy_balance_error\":" << max_balance
       << ",\"max_trajectory_relative_l2_error\":" << max_error << "}\n";
}

void execute(const std::filesystem::path& input, const std::filesystem::path& materials_file,
             const std::filesystem::path& output, int repeats) {
  std::filesystem::create_directories(output);
  auto data = model::readGmsh2D(input, 2, readMaterials(materials_file),
                                {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}});
  model::ModelUnstruct<float, int> mesh(data);
  Parameters p;
  std::ifstream simulation(input.parent_path() / "simulation.txt");
  if (!(simulation >> p.duration >> p.x >> p.z >> p.sigma >> p.receiver_depth >> p.cfl) ||
      !(p.duration > 0 && p.sigma > 0 && p.cfl > 0))
    throw std::runtime_error("Invalid simulation.txt");
  p.absorbing = p.free_top = p.double_reference = true;
  auto probe = solver::buildSEM2DLocalTimeStepPlan(mesh, p.cfl, 1e-9f);
  float estimate = std::numeric_limits<float>::infinity();
  for (auto const& e : probe.element_steps) estimate = std::min(estimate, e.estimate);
  int steps = 6 * static_cast<int>(std::ceil(p.duration / (2 * estimate) / 6));
  p.macro_dt = p.duration / steps;
  auto plan = solver::buildSEM2DLocalTimeStepPlan(mesh, p.cfl, p.macro_dt);
  exportMesh(mesh, plan, output);
  vectorReal initial("Gmsh initial pressure", mesh.getNumberOfNodes());
  auto h = Kokkos::create_mirror_view(initial);
  double top = data.oz_ + data.lz_;
  for (int i = 0; i < mesh.getNumberOfNodes(); ++i) {
    double x = (mesh.nodeCoord(i, 0) - p.x) / p.sigma, z = (mesh.nodeCoord(i, 2) - p.z) / p.sigma,
           zi = (mesh.nodeCoord(i, 2) - (2 * top - p.z)) / p.sigma;
    h(i) = mesh.isFreeSurface(i) ? 0 : std::exp(-.5 * (x * x + z * z)) - std::exp(-.5 * (x * x + zi * zi));
  }
  Kokkos::deep_copy(initial, h);
  Global mass_solver;
  mass_solver.computeFEInit(mesh, {0, 0, 0}, false, 0);
  auto mass = hostField(mass_solver.getMassMatrixAcoustic());
  double spectral_limit = spectralEstimate(mesh, output);
  std::cout << "Imported " << mesh.getNumberOfElements() << " quads, " << mesh.getNumberOfNodes()
            << " GLL nodes; local dt " << p.macro_dt << ", global fine dt " << p.macro_dt / 2 << ", fine elements "
            << plan.fine_elements.size() << std::endl;
  int refinement = 8;
  auto reference = run({"Reference", false, steps * refinement}, mesh, initial, p, true).field;
  auto reference_check = run({"Reference check", false, steps * refinement / 2}, mesh, initial, p, true).field;
  double reference_delta = relativeError(reference_check, reference, mass);
  if (reference_delta > .001) {
    refinement = 16;
    reference_check = std::move(reference);
    reference = run({"Reference refined", false, steps * refinement}, mesh, initial, p, true).field;
    reference_delta = relativeError(reference_check, reference, mass);
  }
  if (!(reference_delta < .001)) throw std::runtime_error("Gmsh reference refinement discrepancy exceeds 0.1%");
  std::vector<Case> cases{{"Global_fine", false, 2 * steps}, {"Global_macro", false, steps}, {"Local", true, steps}};
  for (int warmup = 0; warmup < 2; ++warmup)
    for (const auto& c : cases) run(c, mesh, initial, p);
  std::ofstream raw(output / "timing_samples.csv"), summary(output / "timing_summary.csv");
  for (auto* f : {&raw, &summary}) {
    f->exceptions(std::ios::failbit | std::ios::badbit);
    *f << std::setprecision(12);
  }
  raw << "case,trial,dt_s,steps,setup_ms,loop_ms,relative_l2_error\n";
  summary << "case,dt_s,steps,setup_median_ms,loop_median_ms,loop_q25_ms,loop_q75_ms,relative_l2_error\n";
  std::vector<std::vector<double>> setups(3), loops(3), errors(3);
  std::vector<int> order{0, 1, 2};
  std::mt19937 rng(20261110);
  for (int trial = 0; trial < repeats; ++trial) {
    std::shuffle(order.begin(), order.end(), rng);
    for (int j : order) {
      auto s = run(cases[j], mesh, initial, p);
      double error = relativeError(s.field, reference, mass);
      raw << cases[j].name << ',' << trial << ',' << p.duration / cases[j].steps << ',' << cases[j].steps << ','
          << s.setup_ms << ',' << s.loop_ms << ',' << error << '\n';
      setups[j].push_back(s.setup_ms);
      loops[j].push_back(s.loop_ms);
      errors[j].push_back(error);
    }
    std::cout << "Timing repeat " << trial + 1 << '/' << repeats << std::endl;
  }
  for (int j = 0; j < 3; ++j) {
    summary << cases[j].name << ',' << p.duration / cases[j].steps << ',' << cases[j].steps << ','
            << quantile(setups[j], .5) << ',' << quantile(loops[j], .5) << ',' << quantile(loops[j], .25) << ','
            << quantile(loops[j], .75) << ',' << quantile(errors[j], .5) << '\n';
    std::cout << cases[j].name << ": " << quantile(loops[j], .5) << " ms; error " << quantile(errors[j], .5)
              << std::endl;
  }
  diagnostics(mesh, initial, p, steps, refinement, mass, output);
  std::ofstream protocol(output / "protocol.json");
  protocol << std::setprecision(12) << "{\"elements\":" << mesh.getNumberOfElements()
           << ",\"nodes\":" << mesh.getNumberOfNodes() << ",\"sem_order\":2,\"cfl\":" << p.cfl
           << ",\"macro_dt_s\":" << p.macro_dt << ",\"fine_dt_s\":" << p.macro_dt / 2
           << ",\"duration_s\":" << p.duration << ",\"fine_elements\":" << plan.fine_elements.size()
           << ",\"support_elements\":" << plan.fine_support_elements.size()
           << ",\"reference_steps_per_macro\":" << refinement << ",\"reference_delta_l2\":" << reference_delta
           << ",\"global_spectral_limit_estimate_s\":" << spectral_limit << ",\"repetitions\":" << repeats
           << ",\"warmups\":2,\"seed\":20261110,\"backend\":\"" << Kokkos::DefaultExecutionSpace::name()
           << "\",\"timing_scope\":\"propagation loop only; excludes setup, diagnostics, IO and host copies\"}\n";
}
}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int result = 0;
  try {
    if (argc < 4)
      throw std::invalid_argument("Usage: sem2d_gmsh_benchmark mesh.msh materials.csv output_directory [repetitions]");
    int repeats = argc > 4 ? std::stoi(argv[4]) : 7;
    if (repeats < 3) throw std::invalid_argument("Use at least three repetitions");
    execute(argv[1], argv[2], argv[3], repeats);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    result = 1;
  }
  Kokkos::finalize();
  return result;
}
