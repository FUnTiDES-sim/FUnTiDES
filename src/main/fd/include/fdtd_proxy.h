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

#include <memory>

#include "args_parse.h"
#include "fdtd_grids.h"
#include "fdtd_io.h"
#include "fdtd_kernels.h"
#include "fdtd_options.h"
#include "fdtd_solver.h"
#include "fdtd_source_receivers.h"
#include "fdtd_stencils.h"
#include "read_sepfile.h"
#include "utils.h"

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
 *   FdtdOptions options = ParseCommandLineArgs(argc, argv);
 *   FdtdProxy proxy(options);
 *   proxy.InitFdtd();
 *   proxy.Run();
 * @endcode
 */
class FdtdProxy {
 public:
  /**
   * @brief Constructs an FDTD proxy with the given simulation options.
   * @param opt Configuration options for the simulation including grid size,
   *            time stepping parameters, and I/O settings.
   */
  explicit FdtdProxy(const FdtdOptions& opt);

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
   * @brief Initializes the seismic source term.
   *
   * Configures the source wavelet (e.g., Ricker wavelet) and its spatial
   * location within the computational domain.
   *
   * @post Source parameters are set and ready for injection during
   *       time stepping.
   */
  void InitSource();

  // Simulation configuration
  FdtdOptions opt_;

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
  int source_order_;       ///< Spatial derivative order for source injection
  float source_frequency_; ///< Dominant frequency (Hz) of source wavelet
  float velocity_min_;     ///< Minimum velocity in model (m/s)
  float velocity_max_;     ///< Maximum velocity in model (m/s)
  float wavelength_max_;   ///< Maximum wavelength for stability analysis

  // Source location (grid indices)
  int source_x_ = -1;
  int source_y_ = -1;
  int source_z_ = -1;

  // Core simulation components
  model::fdgrid::FdtdGrids grids_;
  FdtdStencils stencils_;
  FdtdKernels kernels_;
  FdtdIo io_;
  SolverUtils utils_;
  FdtdSolver solver_;
  FdtdSourceReceivers source_receivers_;
};

#endif  // SRC_MAIN_FD_INCLUDE_FDTD_PROXY_HPP_
