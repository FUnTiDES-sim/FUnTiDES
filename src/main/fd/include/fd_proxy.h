//************************************************************************
// Finite Difference Time Domain (FDTD) Acoustic Simulation
// Version 0.0.1
//
// fdtd_proxy.hpp: Main interface for the FDTD proxy application
//
// This header defines the primary orchestration class for finite difference
// time domain acoustic wave propagation simulations. The FdtdProxy class
// manages initialization, coordinate grid setup, time stepping, and
// output generation for forward modeling applications in geophysics.
//
// Copyright (c) 2025
// License: [Specify license here]
//************************************************************************

#ifndef SRC_MAIN_FD_INCLUDE_FDTD_PROXY_HPP_
#define SRC_MAIN_FD_INCLUDE_FDTD_PROXY_HPP_

#include <chrono>
#include <memory>

#include "args_parse.h"
#include "fd_abckernels.h"
#include "fd_grids.h"
#include "fd_io.h"
#include "fd_kernels.h"
#include "fd_options.h"
#include "fd_solver.h"
#include "fd_source_receivers.h"
#include "fd_stencils.h"
#include "utils.h"

namespace fdtd
{
/**
 * @class FdtdProxy
 * @brief Main orchestration class for FDTD acoustic wave propagation
 *        simulation.
 *
 * This class coordinates all components of the FDTD simulation including:
 * - Grid initialization and spatial discretization
 * - Source and receiver configuration
 * - Time-stepping loop management
 * - I/O operations for results
 *
 * Typical usage:
 * @code
 *   fdtd::FdtdOptions options = ParseCommandLineArgs(argc, argv);
 *   FdtdProxy proxy(options);
 *   proxy.InitFdtd();
 *   proxy.Run();
 * @endcode
 */
class FdtdProxy
{
 public:
  /**
   * @brief Constructs an FDTD proxy with the given simulation options.
   * @param opt Configuration options for the simulation including grid size,
   *            time stepping parameters, and I/O settings.
   */
  explicit FdtdProxy(const fdtd::options::FdtdOptions& opt);

  /**
   * @brief Destructor for the FdtdProxy class.
   */
  ~FdtdProxy() = default;

  /**
   * @brief Initializes the FDTD simulation environment.
   *
   * Performs the following initialization steps:
   * - Allocates and initializes computational grids
   * - Configures finite difference stencils
   * - Sets up source and receiver geometries
   * - Prepares I/O subsystems
   *
   * @pre Constructor has been called with valid options.
   * @post All internal state is initialized and ready for Run().
   */
  void InitFdtd();

  /**
   * @brief Executes the main time-stepping loop of the simulation.
   *
   * Advances the wave equation solution through time using the configured
   * numerical scheme. Writes output at specified intervals according to
   * the I/O configuration.
   *
   * @pre InitFdtd() must be called before this method.
   * @post Simulation results are written to disk; performance metrics
   *       may be printed to stdout.
   */
  void Run();

 private:
  /**
   * @brief Initializes the computational grid geometry.
   */
  void InitializeGrid();

  /**
   * @brief Initializes finite difference stencil coefficients.
   */
  void InitializeStencils();

  /**
   * @brief Initializes velocity model parameters and time stepping.
   */
  void InitializeVelocityModel();

  /**
   * @brief Initializes model arrays (velocity, density, etc.).
   */
  void InitializeModelArrays();

  /**
   * @brief Allocates and initializes wavefield arrays.
   */
  void InitializeWavefieldArrays();

  /**
   * @brief Configures seismic source parameters and position.
   */
  void InitializeSource();

  /**
   * @brief Initializes absorbing boundary conditions.
   */
  void InitializeBoundaries();

  /**
   * @brief Initializes the seismic source wavelet.
   *
   * Computes the source time function (e.g., Ricker wavelet) and stores
   * it for injection during time stepping.
   *
   * @post Source wavelet is computed and ready for injection.
   */
  void InitSource();

  /**
   * @brief Prints performance metrics after simulation completion.
   * @param total_compute_time Accumulated computation time
   * @param total_output_time Accumulated I/O time
   */
  void PrintPerformanceMetrics(
      const std::chrono::nanoseconds& total_compute_time,
      const std::chrono::nanoseconds& total_output_time) const;

  // Simulation configuration
  fdtd::options::FdtdOptions opt_;

  // Grid indexing for time integration (current and next time level)
  int time_index_current_ = 0;
  int time_index_next_ = 1;

  // Finite difference stencil sizes
  int num_coefs_x_;
  int num_coefs_y_;
  int num_coefs_z_;

  // Time integration parameters
  int num_time_samples_;
  float time_step_;
  float time_max_;

  // Source configuration
  int source_order_;        ///< Spatial derivative order for source injection
  float source_frequency_;  ///< Dominant frequency (Hz) of source wavelet
  float velocity_min_;      ///< Minimum velocity in model (m/s)
  float velocity_max_;      ///< Maximum velocity in model (m/s)
  float wavelength_max_;    ///< Maximum wavelength for stability analysis

  // Source location (grid indices)
  int source_x_ = -1;
  int source_y_ = -1;
  int source_z_ = -1;

  // Core simulation components
  model::fdgrid::FdtdGrids grids_;
  fdtd::stencils::FdtdStencils stencils_;
  fdtd::kernel::FdtdKernels kernels_;
  fdtd::abckernel::FdtdAbcKernels abckernels_;
  fdtd::io::FdtdIo io_;
  SolverUtils utils_;
  FdtdSolver solver_;
  FdtdSourceReceivers source_receivers_;
};

}  // namespace fdtd
#endif  // SRC_MAIN_FD_INCLUDE_FDTD_PROXY_HPP_
