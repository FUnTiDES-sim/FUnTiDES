#pragma once

#include <cxxopts.hpp>
#include <json.hpp>
#include <stdexcept>
#include <string>

class SemProxyOptions
{
 public:
  // Defaults
  int order = 2;
  int ex = 50, ey = 50, ez = 50;
  float lx = 2000.f, ly = 2000.f, lz = 2000.f;
  float srcx = 1010.f, srcy = 1010.f, srcz = 1010.f;
  float rcvx = 1410.f, rcvy = 1010.f, rcvz = 1010.f;
  std::string implem = "makutu";  // makutu|shiva
  std::string method = "sem";     // sem|dg
  std::string mesh = "cartesian";
  float dt = 0.006;
  float timemax = 0.7;
  bool autodt = false;
  // snapshots
  bool snapshots = false;
  int snap_time_interval = 10;
  // sponge boundaries parameters
  float boundaries_size = 0;
  bool surface_sponge = false;
  float taper_delta = 0.015;
  // Boolean to tell if the model is charged on nodes or on element
  bool isModelOnNodes = false;
  bool isElastic = false;

  void validate() const
  {
    if (order < 1) throw std::runtime_error("order must be >= 1");
    if (ex <= 0 || ey <= 0 || ez <= 0)
      throw std::runtime_error("ex/ey/ez must be > 0");
    if (lx <= 0 || ly <= 0 || lz <= 0)
      throw std::runtime_error("lx/ly/lz must be > 0");
  }

  // Load from JSON file
  void load_from_json(const std::string& json_path)
  {
    std::ifstream file(json_path);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open JSON file: " + json_path);
    }

    nlohmann::json j;
    try {
      file >> j;
    } catch (const nlohmann::json::exception& e) {
      throw std::runtime_error("JSON parsing error: " + std::string(e.what()));
    }

    // Simulation parameters
    if (j.contains("simulation")) {
      auto& sim = j["simulation"];
      if (sim.contains("order")) order = sim["order"];
      if (sim.contains("method")) method = sim["method"];
      if (sim.contains("implementation")) implem = sim["implementation"];
      if (sim.contains("mesh")) mesh = sim["mesh"];
      if (sim.contains("dt")) dt = sim["dt"];
      if (sim.contains("timemax")) timemax = sim["timemax"];
      if (sim.contains("autodt")) autodt = sim["autodt"];
    }

    // Domain parameters
    if (j.contains("domain")) {
      auto& dom = j["domain"];
      if (dom.contains("ex")) ex = dom["ex"];
      if (dom.contains("ey")) ey = dom["ey"];
      if (dom.contains("ez")) ez = dom["ez"];
      if (dom.contains("lx")) lx = dom["lx"];
      if (dom.contains("ly")) ly = dom["ly"];
      if (dom.contains("lz")) lz = dom["lz"];
    }

    // Source position
    if (j.contains("source")) {
      auto& src = j["source"];
      if (src.contains("x")) srcx = src["x"];
      if (src.contains("y")) srcy = src["y"];
      if (src.contains("z")) srcz = src["z"];
    }

    // Receiver position
    if (j.contains("receiver")) {
      auto& rcv = j["receiver"];
      if (rcv.contains("x")) rcvx = rcv["x"];
      if (rcv.contains("y")) rcvy = rcv["y"];
      if (rcv.contains("z")) rcvz = rcv["z"];
    }

    // Snapshots
    if (j.contains("snapshots")) {
      auto& snap = j["snapshots"];
      if (snap.contains("enabled")) snapshots = snap["enabled"];
      if (snap.contains("time_interval")) snap_time_interval = snap["time_interval"];
    }

    // Boundaries
    if (j.contains("boundaries")) {
      auto& bound = j["boundaries"];
      if (bound.contains("surface_sponge")) surface_sponge = bound["surface_sponge"];
      if (bound.contains("size")) boundaries_size = bound["size"];
      if (bound.contains("taper_delta")) taper_delta = bound["taper_delta"];
    }

    // Model parameters
    if (j.contains("model")) {
      auto& model = j["model"];
      if (model.contains("is_on_nodes")) isModelOnNodes = model["is_on_nodes"];
      if (model.contains("is_elastic")) isElastic = model["is_elastic"];
    }
  }

  // Bind CLI flags to this instance (no --help here)
  static void bind_cli(cxxopts::Options& opts, SemProxyOptions& o)
  {
    opts.add_options()("o,order", "Order of approximation",
                       cxxopts::value<int>(o.order))(
        "ex", "Number of elements on X (Cartesian mesh)",
        cxxopts::value<int>(o.ex))("ey",
                                   "Number of elements on Y (Cartesian mesh)",
                                   cxxopts::value<int>(o.ey))(
        "ez", "Number of elements on Z (Cartesian mesh)",
        cxxopts::value<int>(o.ez))("lx", "Domain size X (Cartesian)",
                                   cxxopts::value<float>(o.lx))(
        "ly", "Domain size Y (Cartesian)", cxxopts::value<float>(o.ly))(
        "lz", "Domain size Z (Cartesian)", cxxopts::value<float>(o.lz))(
        "implem", "Implementation: makutu|shiva",
        cxxopts::value<std::string>(o.implem))(
        "method", "Method: sem|dg", cxxopts::value<std::string>(o.method))(
        "mesh", "Mesh: cartesian|ucartesian",
        cxxopts::value<std::string>(o.mesh))(
        "dt", "Time step selection in s (default = 0.001s)",
        cxxopts::value<float>(o.dt))(
        "timemax", "Duration of the simulation in s (default = 1.5s)",
        cxxopts::value<float>(o.timemax))(
        "auto-dt", "Select automatique dt via CFL equation.",
        cxxopts::value<bool>(o.autodt))("s,snapshots", "Enable snapshot.",
                                        cxxopts::value<bool>(o.snapshots))(
        "snap-interval",
        "Interval on iteration between two snapshots. (default=10)",
        cxxopts::value<int>(o.snap_time_interval))(
        "boundaries-size", "Size of absorbing boundaries (meters)",
        cxxopts::value<float>(o.boundaries_size))(
        "sponge-surface", "Considere the surface's nodes as non sponge nodes",
        cxxopts::value<bool>(o.surface_sponge))(
        "taper-delta", "Taper delta for sponge boundaries value",
        cxxopts::value<float>(o.taper_delta))(
        "is-model-on-nodes",
        "Boolean to tell if the model is charged on nodes (true) or on element "
        "(false)",
        cxxopts::value<bool>(o.isModelOnNodes))(
        "is-elastic", "Elastic simulation", cxxopts::value<bool>(o.isElastic));
  }
};
