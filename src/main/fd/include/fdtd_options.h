#pragma once

#include <cxxopts.hpp>
#include <stdexcept>
#include <string>

class fdtd_options
{
 public:
  // Grid parameters
  struct GridParams
  {
    int nx{200}, ny{200}, nz{200};
    float dx{10.f}, dy{10.f}, dz{10.f};
    std::string mesh{"cartesian"};
  } grid;

  // Stencil configuration
  struct StencilParams
  {
    int lx{4}, ly{4}, lz{4};
    std::string implem{"remez"};  // remez|taylor
  } stencil;

  // Source configuration
  struct SourceParams
  {
    int xs{-1}, ys{-1}, zs{-1};  // -1 means center of grid
    float f0{10.f};
    int sourceOrder{2};
  } source;

  // Velocity model
  struct VelocityParams
  {
    float vmin{1500.f}, vmax{4500.f};
    std::string fileModel{""};
    bool usefilemodel{false};
  } velocity;

  // Time stepping
  struct TimeParams
  {
    float timeStep{0.f};  // 0 means auto-compute from CFL
    float timeMax{1.f};
    std::string method{"FDTD"};
  } time;

  // Boundary conditions
  struct BoundaryParams
  {
    bool usePML{false};
    int pmlSize{20};
    int spongeSize{20};
    float spongeAlpha{0.015f};
  } boundary;

  // Output configuration
  struct OutputParams
  {
    bool saveSnapShots{false};
    int snapShotInterval{10};
  } output;

  void validate() const
  {
    validateGrid();
    validateStencil();
    validateSource();
    validateVelocity();
    validateTime();
    validateBoundary();
    validateOutput();
  }

  static void bind_cli(cxxopts::Options& opts, fdtd_options& o)
  {
    // Grid options
    opts.add_options("Grid")("nx", "Number of grid points on X",
                             cxxopts::value<int>(o.grid.nx))(
        "ny", "Number of grid points on Y", cxxopts::value<int>(o.grid.ny))(
        "nz", "Number of grid points on Z", cxxopts::value<int>(o.grid.nz))(
        "dx", "Grid spacing on X", cxxopts::value<float>(o.grid.dx))(
        "dy", "Grid spacing on Y", cxxopts::value<float>(o.grid.dy))(
        "dz", "Grid spacing on Z", cxxopts::value<float>(o.grid.dz))(
        "mesh", "Mesh type (cartesian)",
        cxxopts::value<std::string>(o.grid.mesh));

    // Stencil options
    opts.add_options("Stencil")("lx", "Half stencil size X",
                                cxxopts::value<int>(o.stencil.lx))(
        "ly", "Half stencil size Y", cxxopts::value<int>(o.stencil.ly))(
        "lz", "Half stencil size Z", cxxopts::value<int>(o.stencil.lz))(
        "implem", "Stencil implementation (remez|taylor)",
        cxxopts::value<std::string>(o.stencil.implem));

    // Source options
    opts.add_options("Source")("xs", "Source X position",
                               cxxopts::value<int>(o.source.xs))(
        "ys", "Source Y position", cxxopts::value<int>(o.source.ys))(
        "zs", "Source Z position", cxxopts::value<int>(o.source.zs))(
        "f0", "Peak frequency", cxxopts::value<float>(o.source.f0))(
        "sourceOrder", "Source time derivative order",
        cxxopts::value<int>(o.source.sourceOrder));

    // Velocity options
    opts.add_options("Velocity")("vmin", "Minimum velocity",
                                 cxxopts::value<float>(o.velocity.vmin))(
        "vmax", "Maximum velocity", cxxopts::value<float>(o.velocity.vmax))(
        "fileModel", "Velocity model file",
        cxxopts::value<std::string>(o.velocity.fileModel));

    // Time stepping options
    opts.add_options("Time")("timeStep", "Time step (0=auto)",
                             cxxopts::value<float>(o.time.timeStep))(
        "timeMax", "Maximum simulation time",
        cxxopts::value<float>(o.time.timeMax))(
        "method", "Time stepping method",
        cxxopts::value<std::string>(o.time.method));

    // Boundary options
    opts.add_options("Boundary")("usePML", "Use PML boundaries",
                                 cxxopts::value<bool>(o.boundary.usePML))(
        "pmlSize", "PML layer thickness",
        cxxopts::value<int>(o.boundary.pmlSize))(
        "spongeSize", "Sponge layer thickness",
        cxxopts::value<int>(o.boundary.spongeSize))(
        "spongeAlpha", "Sponge strength",
        cxxopts::value<float>(o.boundary.spongeAlpha));

    // Output options
    opts.add_options("Output")("saveSnapShots", "Enable snapshot saving",
                               cxxopts::value<bool>(o.output.saveSnapShots))(
        "snapShotInterval", "Steps between snapshots",
        cxxopts::value<int>(o.output.snapShotInterval));

    // Help
    opts.add_options()("h,help", "Print help");
  }

 private:
  void validateGrid() const
  {
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0)
      throw std::runtime_error("Grid dimensions must be positive");
    if (grid.dx <= 0 || grid.dy <= 0 || grid.dz <= 0)
      throw std::runtime_error("Grid spacing must be positive");
  }

  void validateStencil() const
  {
    if (stencil.lx <= 0 || stencil.ly <= 0 || stencil.lz <= 0)
      throw std::runtime_error("Stencil sizes must be positive");
  }

  void validateSource() const
  {
    if (source.f0 <= 0)
      throw std::runtime_error("Source frequency must be positive");
    if (source.sourceOrder < 1)
      throw std::runtime_error("Source order must be >= 1");
  }

  void validateVelocity() const
  {
    if (velocity.vmin <= 0 || velocity.vmax <= 0)
      throw std::runtime_error("Velocities must be positive");
    if (velocity.vmin >= velocity.vmax)
      throw std::runtime_error("vmin must be less than vmax");
  }

  void validateTime() const
  {
    if (time.timeMax <= 0)
      throw std::runtime_error("Maximum time must be positive");
    if (time.timeStep < 0)
      throw std::runtime_error("Time step cannot be negative");
  }

  void validateBoundary() const
  {
    if (boundary.pmlSize < 0)
      throw std::runtime_error("PML size cannot be negative");
    if (boundary.spongeSize < 0)
      throw std::runtime_error("Sponge size cannot be negative");
    if (boundary.spongeAlpha < 0)
      throw std::runtime_error("Sponge alpha cannot be negative");
  }

  void validateOutput() const
  {
    if (output.snapShotInterval <= 0)
      throw std::runtime_error("Snapshot interval must be positive");
  }
};
