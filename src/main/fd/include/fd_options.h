//************************************************************************
// Finite Difference Time Domain (FDTD) Acoustic Simulation
// Version 0.0.1
//
// src/main/fd/include/fdtd_options.hpp
//
// Configuration options for FDTD simulation
//
// This header defines the FdtdOptions class which encapsulates all
// configuration parameters for FDTD simulations including grid geometry,
// stencil configuration, source parameters, velocity models, time stepping,
// boundary conditions, and output settings. It provides command-line
// argument binding and validation functionality.
//************************************************************************

#ifndef SRC_MAIN_FD_INCLUDE_FDTD_OPTIONS_HPP_
#define SRC_MAIN_FD_INCLUDE_FDTD_OPTIONS_HPP_

#include <cxxopts.hpp>
#include <stdexcept>
#include <string>

namespace fdtd
{

namespace options
{
/**
 * @class FdtdOptions
 * @brief Configuration container for FDTD simulation parameters.
 *
 * This class organizes all simulation parameters into logical groups
 * (grid, stencil, source, velocity, time, boundary, output) and provides
 * validation and command-line interface binding functionality.
 *
 * Typical usage:
 * @code
 *   FdtdOptions options;
 *   cxxopts::Options cli("fdtd", "FDTD Acoustic Simulator");
 *   FdtdOptions::BindCli(cli, options);
 *   cli.parse(argc, argv);
 *   options.Validate();
 * @endcode
 */
class FdtdOptions
{
 public:
  /**
   * @struct GridParams
   * @brief Parameters defining the computational grid geometry.
   */
  struct GridParams
  {
    int nx{200};                    ///< Number of grid points in X direction
    int ny{200};                    ///< Number of grid points in Y direction
    int nz{200};                    ///< Number of grid points in Z direction
    float dx{10.f};                 ///< Grid spacing in X direction (meters)
    float dy{10.f};                 ///< Grid spacing in Y direction (meters)
    float dz{10.f};                 ///< Grid spacing in Z direction (meters)
    std::string mesh{"cartesian"};  ///< Mesh type (currently only cartesian)
  } grid;

  /**
   * @struct StencilParams
   * @brief Parameters for finite difference stencil configuration.
   */
  struct StencilParams
  {
    int lx{4};                    ///< Half-width of stencil in X direction
    int ly{4};                    ///< Half-width of stencil in Y direction
    int lz{4};                    ///< Half-width of stencil in Z direction
    std::string implem{"remez"};  ///< Stencil implementation (remez|taylor)
  } stencil;

  /**
   * @struct SourceParams
   * @brief Parameters defining the seismic source configuration.
   */
  struct SourceParams
  {
    int xs{-1};           ///< Source X position (-1 = grid center)
    int ys{-1};           ///< Source Y position (-1 = grid center)
    int zs{-1};           ///< Source Z position (-1 = grid center)
    float f0{10.f};       ///< Peak frequency (Hz)
    int source_order{2};  ///< Time derivative order of source wavelet
  } source;

  /**
   * @struct VelocityParams
   * @brief Parameters for velocity model configuration.
   */
  struct VelocityParams
  {
    float vmin{1500.f};          ///< Minimum velocity in model (m/s)
    float vmax{4500.f};          ///< Maximum velocity in model (m/s)
    std::string file_model{""};  ///< Path to velocity model file
    bool use_file_model{false};  ///< Flag to use file-based model
  } velocity;

  /**
   * @struct TimeParams
   * @brief Parameters for time integration configuration.
   */
  struct TimeParams
  {
    float time_step{0.f};        ///< Time step size (0 = auto-compute from CFL)
    float time_max{1.f};         ///< Maximum simulation time (seconds)
    std::string method{"FDTD"};  ///< Time stepping method
  } time;

  /**
   * @struct BoundaryParams
   * @brief Parameters for absorbing boundary conditions.
   */
  struct BoundaryParams
  {
    bool use_pml{true};          ///< Enable PML boundary conditions
    bool use_sponge{false};      ///< Enable PML boundary conditions
    int pml_size{4};             ///< Thickness of PML layers (grid points)
    int sponge_size{20};         ///< Thickness of sponge layers (grid points)
    float sponge_alpha{0.015f};  ///< Damping coefficient for sponge
  } boundary;

  /**
   * @struct OutputParams
   * @brief Parameters controlling simulation output.
   */
  struct OutputParams
  {
    bool save_snapshots{false};  ///< Enable wavefield snapshot saving
    int snapshot_interval{10};   ///< Time steps between snapshots
  } output;

  /**
   * @brief Validates all configuration parameters.
   *
   * Checks that all parameters are within valid ranges and logically
   * consistent. Throws std::runtime_error if validation fails.
   *
   * @throws std::runtime_error if any parameter is invalid
   */
  void Validate() const
  {
    ValidateGrid();
    ValidateStencil();
    ValidateSource();
    ValidateVelocity();
    ValidateTime();
    ValidateBoundary();
    ValidateOutput();
  }

  /**
   * @brief Binds configuration options to command-line parser.
   *
   * Configures a cxxopts::Options object to parse command-line arguments
   * and populate the provided FdtdOptions instance.
   *
   * @param opts The cxxopts::Options object to configure
   * @param options The FdtdOptions instance to populate with parsed values
   */
  static void BindCli(cxxopts::Options& opts, FdtdOptions& options)
  {
    // Grid options
    opts.add_options("Grid")("nx", "Number of grid points in X direction",
                             cxxopts::value<int>(options.grid.nx))(
        "ny", "Number of grid points in Y direction",
        cxxopts::value<int>(options.grid.ny))(
        "nz", "Number of grid points in Z direction",
        cxxopts::value<int>(options.grid.nz))(
        "dx", "Grid spacing in X direction (meters)",
        cxxopts::value<float>(options.grid.dx))(
        "dy", "Grid spacing in Y direction (meters)",
        cxxopts::value<float>(options.grid.dy))(
        "dz", "Grid spacing in Z direction (meters)",
        cxxopts::value<float>(options.grid.dz))(
        "mesh", "Mesh type (cartesian)",
        cxxopts::value<std::string>(options.grid.mesh));

    // Stencil options
    opts.add_options("Stencil")("lx", "Half-width of stencil in X direction",
                                cxxopts::value<int>(options.stencil.lx))(
        "ly", "Half-width of stencil in Y direction",
        cxxopts::value<int>(options.stencil.ly))(
        "lz", "Half-width of stencil in Z direction",
        cxxopts::value<int>(options.stencil.lz))(
        "implem", "Stencil implementation method (remez|taylor)",
        cxxopts::value<std::string>(options.stencil.implem));

    // Source options
    opts.add_options("Source")("xs", "Source X position (-1 for grid center)",
                               cxxopts::value<int>(options.source.xs))(
        "ys", "Source Y position (-1 for grid center)",
        cxxopts::value<int>(options.source.ys))(
        "zs", "Source Z position (-1 for grid center)",
        cxxopts::value<int>(options.source.zs))(
        "f0", "Source peak frequency (Hz)",
        cxxopts::value<float>(options.source.f0))(
        "sourceOrder", "Source time derivative order",
        cxxopts::value<int>(options.source.source_order));

    // Velocity options
    opts.add_options("Velocity")("vmin", "Minimum velocity in model (m/s)",
                                 cxxopts::value<float>(options.velocity.vmin))(
        "vmax", "Maximum velocity in model (m/s)",
        cxxopts::value<float>(options.velocity.vmax))(
        "fileModel", "Path to velocity model file",
        cxxopts::value<std::string>(options.velocity.file_model));

    // Time stepping options
    opts.add_options("Time")("timeStep",
                             "Time step size (0 for auto-compute from CFL)",
                             cxxopts::value<float>(options.time.time_step))(
        "timeMax", "Maximum simulation time (seconds)",
        cxxopts::value<float>(options.time.time_max))(
        "method", "Time stepping method (FDTD)",
        cxxopts::value<std::string>(options.time.method));

    // Boundary options
    opts.add_options("Boundary")(
        "usePML", "Enable PML absorbing boundaries",
        cxxopts::value<bool>(options.boundary.use_pml))(
        "pmlSize", "PML layer thickness (grid points)",
        cxxopts::value<int>(options.boundary.pml_size))(
        "useSponge", "Enable sponge absorbing boundaries",
        cxxopts::value<bool>(options.boundary.use_sponge))(
        "spongeSize", "Sponge layer thickness (grid points)",
        cxxopts::value<int>(options.boundary.sponge_size))(
        "spongeAlpha", "Sponge damping coefficient",
        cxxopts::value<float>(options.boundary.sponge_alpha));

    // Output options
    opts.add_options("Output")(
        "saveSnapShots", "Enable wavefield snapshot saving",
        cxxopts::value<bool>(options.output.save_snapshots))(
        "snapShotInterval", "Time steps between snapshots",
        cxxopts::value<int>(options.output.snapshot_interval));

    // Help option
    opts.add_options()("h,help", "Print usage information");
  }

 private:
  void ValidateGrid() const
  {
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0)
    {
      throw std::runtime_error("Grid dimensions must be positive");
    }
    if (grid.dx <= 0.f || grid.dy <= 0.f || grid.dz <= 0.f)
    {
      throw std::runtime_error("Grid spacing must be positive");
    }
  }

  void ValidateStencil() const
  {
    if (stencil.lx <= 0 || stencil.ly <= 0 || stencil.lz <= 0)
    {
      throw std::runtime_error("Stencil half-widths must be positive");
    }
    if (stencil.implem != "remez" && stencil.implem != "taylor")
    {
      throw std::runtime_error(
          "Stencil implementation must be 'remez' or 'taylor'");
    }
  }

  void ValidateSource() const
  {
    if (source.f0 <= 0.f)
    {
      throw std::runtime_error("Source frequency must be positive");
    }
    if (source.source_order < 1)
    {
      throw std::runtime_error("Source order must be >= 1");
    }
  }

  void ValidateVelocity() const
  {
    if (velocity.vmin <= 0.f || velocity.vmax <= 0.f)
    {
      throw std::runtime_error("Velocities must be positive");
    }
    if (velocity.vmin >= velocity.vmax)
    {
      throw std::runtime_error("Minimum velocity must be less than maximum");
    }
  }

  void ValidateTime() const
  {
    if (time.time_max <= 0.f)
    {
      throw std::runtime_error("Maximum simulation time must be positive");
    }
    if (time.time_step < 0.f)
    {
      throw std::runtime_error("Time step cannot be negative");
    }
  }

  void ValidateBoundary() const
  {
    if (boundary.pml_size < 0)
    {
      throw std::runtime_error("PML size cannot be negative");
    }
    if (boundary.sponge_size < 0)
    {
      throw std::runtime_error("Sponge size cannot be negative");
    }
    if (boundary.sponge_alpha < 0.f)
    {
      throw std::runtime_error("Sponge damping coefficient cannot be negative");
    }
    if (!boundary.use_pml && !boundary.use_sponge)
    {
      throw std::runtime_error(
          "At least one absorbing boundary condition must be enabled");
    }
  }

  void ValidateOutput() const
  {
    if (output.snapshot_interval <= 0)
    {
      throw std::runtime_error("Snapshot interval must be positive");
    }
  }
};

}  // namespace options
}  // namespace fdtd

#endif  // SRC_MAIN_FD_INCLUDE_FDTD_OPTIONS_HPP_
