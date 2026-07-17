#ifndef FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_OPTIONS_H_
#define FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_OPTIONS_H_

#pragma once

#include <cxxopts.hpp>
#include <stdexcept>
#include <string>
#include <vector>

class SemProxyOptions {
 public:
  // Defaults
  int order = 2;
  int ex = 50, ey = 50, ez = 50;
  float lx = 2000.f, ly = 2000.f, lz = 2000.f;
  float srcx = 1010.f, srcy = 1010.f, srcz = 1010.f;
  float rcvx = 1310.f, rcvy = 1310.f, rcvz = 1310.f;
  float f0 = 5.0f;       // dominant frequency of the source in Hz
  int ricker_order = 2;  // order of the Ricker wavelet source
  float tpeak = 0.2f;    // peak time of the Ricker wavelet source

  std::string implem = "makutu";  // makutu
  std::string method = "sem";     // sem
  std::string mesh = "cartesian";
  std::string anisotropy = "iso";  // iso|vti|tti
  float dt = 0.006;
  float timemax = 0.7;
  bool autodt = false;
  // snapshots
  bool snapshots = false;
  int snap_time_interval = 20;
  // sponge boundaries parameters
  bool use_sponge = false;
  float boundaries_size = 0;
  bool surface_sponge = false;
  float taper_delta = 0.015;
  // Boolean to tell if the model is charged on nodes or on element
  bool isModelOnNodes = false;
  bool isElastic = false;
  bool isAcoustoElastic = false;
  float acoustoElasticBoundaryZ = 0.f;
  bool free_surface = false;
  std::string model_file{""};
  float qp = -1.0f;  // quality factor for P-waves (<0 = not set)
  float qs = -1.0f;  // quality factor for S-waves (<0 = not set)
  std::vector<float> sls_reference_angular_frequencies{};
  std::vector<float> sls_anelasticity_coefficients{};

  // DAS (Distributed Acoustic Sensing) receiver parameters
  std::string das_type = "none";  ///< none | dipole | strain
  float das_dip = 0.f;            ///< Fiber dip angle in degrees
  float das_azimuth = 0.f;        ///< Fiber azimuth angle in degrees
  float das_gauge_length = 1.f;   ///< Gauge length in meters
  int das_samples = 5;            ///< Number of integration points along fiber

  void validate() const {
    if (order < 1) throw std::runtime_error("order must be >= 1");
    if (ex <= 0 || ey <= 0 || ez <= 0) throw std::runtime_error("ex/ey/ez must be > 0");
    if (lx <= 0 || ly <= 0 || lz <= 0) throw std::runtime_error("lx/ly/lz must be > 0");
  }

  // Bind CLI flags to this instance (no --help here)
  static void bind_cli(cxxopts::Options& opts, SemProxyOptions& o) {
    opts.add_options()("o,order", "Order of approximation", cxxopts::value<int>(o.order))(
        "ex", "Number of elements on X (Cartesian mesh)", cxxopts::value<int>(o.ex))(
        "ey", "Number of elements on Y (Cartesian mesh)", cxxopts::value<int>(o.ey))(
        "ez", "Number of elements on Z (Cartesian mesh)", cxxopts::value<int>(o.ez))(
        "lx", "Domain size X (Cartesian)", cxxopts::value<float>(o.lx))("ly", "Domain size Y (Cartesian)",
                                                                        cxxopts::value<float>(o.ly))(
        "lz", "Domain size Z (Cartesian)", cxxopts::value<float>(o.lz))("implem", "Implementation: makutu",
                                                                        cxxopts::value<std::string>(o.implem))(
        "method", "Method: sem|dg", cxxopts::value<std::string>(o.method))("mesh", "Mesh: cartesian|ucartesian",
                                                                           cxxopts::value<std::string>(o.mesh))(
        "dt", "Time step selection in s (default = 0.001s)", cxxopts::value<float>(o.dt))(
        "timemax", "Duration of the simulation in s (default = 1.5s)", cxxopts::value<float>(o.timemax))(
        "auto-dt", "Select automatique dt via CFL equation.", cxxopts::value<bool>(o.autodt))(
        "s,snapshots", "Enable snapshot.", cxxopts::value<bool>(o.snapshots))(
        "snap-interval", "Interval on iteration between two snapshots. (default=10)",
        cxxopts::value<int>(o.snap_time_interval))("use-sponge", "Enable sponge (taper) absorbing boundaries",
                                                   cxxopts::value<bool>(o.use_sponge))(
        "boundaries-size", "Size of absorbing boundaries (meters)", cxxopts::value<float>(o.boundaries_size))(
        "sponge-surface", "Considere the surface's nodes as non sponge nodes", cxxopts::value<bool>(o.surface_sponge))(
        "taper-delta", "Taper delta for sponge boundaries value", cxxopts::value<float>(o.taper_delta))(
        "is-model-on-nodes",
        "Boolean to tell if the model is charged on nodes (true) or on element "
        "(false)",
        cxxopts::value<bool>(o.isModelOnNodes))("is-elastic", "Elastic simulation", cxxopts::value<bool>(o.isElastic))(
        "is-acousto-elastic", "Acousto-elastic coupled simulation", cxxopts::value<bool>(o.isAcoustoElastic))(
        "acousto-elastic-boundary-z", "Z coordinate of the fluid–solid interface (meters)",
        cxxopts::value<float>(o.acoustoElasticBoundaryZ))(
        "free-surface", "Enable free surface on top boundary (Z+). Default: true",
        cxxopts::value<bool>(o.free_surface))("anisotropy", "Anisotropy type for elastic: iso|vti|tti (default=iso)",
                                              cxxopts::value<std::string>(o.anisotropy))(
        "model-file", "Path to .ftmd heterogeneous model file (optional)", cxxopts::value<std::string>(o.model_file))(
        "das-type", "DAS receiver type: none|dipole|strain (default=none)", cxxopts::value<std::string>(o.das_type))(
        "das-dip", "DAS fiber dip angle in degrees (default=0)", cxxopts::value<float>(o.das_dip))(
        "das-azimuth", "DAS fiber azimuth angle in degrees (default=0)", cxxopts::value<float>(o.das_azimuth))(
        "das-gauge-length", "DAS gauge length in meters (default=1)", cxxopts::value<float>(o.das_gauge_length))(
        "das-samples", "Number of integration points along DAS fiber (default=5)", cxxopts::value<int>(o.das_samples))(
        "f0", "Dominant frequency of the source in Hz", cxxopts::value<float>(o.f0))(
        "tpeak", "Peak time of the Ricker wavelet source", cxxopts::value<float>(o.tpeak))(
        "ricker-order", "Order of the Ricker wavelet source", cxxopts::value<int>(o.ricker_order))(
        "srcx", "Source position X (meters)", cxxopts::value<float>(o.srcx))("srcy", "Source position Y (meters)",
                                                                             cxxopts::value<float>(o.srcy))(
        "srcz", "Source position Z (meters)", cxxopts::value<float>(o.srcz))("rcvx", "Receiver position X (meters)",
                                                                             cxxopts::value<float>(o.rcvx))(
        "rcvy", "Receiver position Y (meters)", cxxopts::value<float>(o.rcvy))("rcvz", "Receiver position Z (meters)",
                                                                               cxxopts::value<float>(o.rcvz))(
        "qp", "Quality factor for P-waves (default: no attenuation)", cxxopts::value<float>(o.qp))(
        "qs", "Quality factor for S-waves (default: no attenuation)", cxxopts::value<float>(o.qs))(
        "sls-reference-angular-frequencies", "Comma-separated SLS reference angular frequencies (rad/s)",
        cxxopts::value<std::vector<float>>(o.sls_reference_angular_frequencies))(
        "sls-anelasticity-coefficients",
        "Comma-separated SLS anelasticity coefficients (same size as "
        "frequencies)",
        cxxopts::value<std::vector<float>>(o.sls_anelasticity_coefficients));
  }
};

#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_PROXY_OPTIONS_H_
