#pragma once

#include <cxxopts.hpp>
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
