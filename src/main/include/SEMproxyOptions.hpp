#pragma once

#include <cxxopts.hpp>
#include <string>
#include <stdexcept>

class SemProxyOptions {
public:
  // Defaults
  int   order = 2;
  int   ex = 50, ey = 50, ez = 50;
  float lx = 2000.f, ly = 2000.f, lz = 2000.f;
  std::string implem = "optim";   // classic|optim|geos|shiva
  std::string method = "sem";     // sem|dg
  std::string mesh = "cartesian";
  // snapshots
  bool snapshots = false;
  int snap_time_interval = 10;
  std::string snap_folder = "snapshots";

  void validate() const {
    if (order < 1) throw std::runtime_error("order must be >= 1");
    if (ex <= 0 || ey <= 0 || ez <= 0) throw std::runtime_error("ex/ey/ez must be > 0");
    if (lx <= 0 || ly <= 0 || lz <= 0) throw std::runtime_error("lx/ly/lz must be > 0");
  }

  // Bind CLI flags to this instance (no --help here)
  static void bind_cli(cxxopts::Options& opts, SemProxyOptions& o) {
    opts.add_options()
      ("o,order", "Order of approximation",
          cxxopts::value<int>(o.order))
      ("ex", "Number of elements on X (Cartesian mesh)",
          cxxopts::value<int>(o.ex))
      ("ey", "Number of elements on Y (Cartesian mesh)",
          cxxopts::value<int>(o.ey))
      ("ez", "Number of elements on Z (Cartesian mesh)",
          cxxopts::value<int>(o.ez))
      ("lx", "Domain size X (Cartesian)",
          cxxopts::value<float>(o.lx))
      ("ly", "Domain size Y (Cartesian)",
          cxxopts::value<float>(o.ly))
      ("lz", "Domain size Z (Cartesian)",
          cxxopts::value<float>(o.lz))
      ("implem", "Implementation: classic|optim|geos|shiva",
          cxxopts::value<std::string>(o.implem))
      ("method", "Method: sem|dg",
          cxxopts::value<std::string>(o.method))
      ("mesh", "Mesh: cartesian|ucartesian",
          cxxopts::value<std::string>(o.mesh))
      ("s,snapshots", "Enable snapshot.",
          cxxopts::value<bool>(o.snapshots))
      ("snap-folder", "Folder where to save snapshots. (default=snapshots)",
          cxxopts::value<std::string>(o.snap_folder))
      ("snap-interval", "Interval on iteration between two snapshots. (default=10)",
          cxxopts::value<int>(o.snap_time_interval));
  }
};
