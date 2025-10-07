/**
 * @file fdtd_grids.h
 * @brief FDTD grid and model management
 *
 * This file provides the core grid infrastructure for finite-difference
 * time-domain (FDTD) acoustic wave simulations. It manages:
 * - 3D computational grid geometry and spacing
 * - Physical property models (velocity, density)
 * - Boundary condition regions
 * - Model initialization from files or synthetic models
 *
 * @license
 * Copyright (C) 2025.
 * This software is governed by the CeCILL license under French law.
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_

#include <data_type.h>
#include <fdtd_options.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "fdtd_macros.h"

namespace model
{
namespace fdgrid
{

/**
 * @brief FDTD computational grid and model container
 *
 * Manages the complete 3D computational domain for FDTD simulations including:
 * - Grid dimensions (nx, ny, nz) and spatial sampling (dx, dy, dz)
 * - Physical property models (P-wave velocity vp, S-wave velocity vs, density
 * rho)
 * - Precomputed spatial derivative coefficients for performance
 * - Absorbing boundary layer (PML/sponge) geometry
 *
 * **Usage Pattern:**
 * @code
 * FdtdGrids grids;
 * grids.InitGrid(options);        // Step 1: Set up geometry
 * grids.InitModelArrays(options); // Step 2: Allocate and initialize models
 * @endcode
 *
 * **Coordinate System:**
 * - Storage order: column-major for Z (fastest), then Y, then X (slowest)
 * - Index calculation: index = nz*ny*i + nz*j + k
 *
 * @note Grid must be initialized via InitGrid() before InitModelArrays()
 * @note Thread-safe for read access after initialization
 */
class FdtdGrids
{
 public:
  FdtdGrids() = default;

  /**
   * @brief Initialize grid geometry from configuration options
   *
   * Sets up the computational domain dimensions and spatial sampling.
   * If a model file is specified in options, grid parameters are loaded
   * from that file, overriding the values in the options structure.
   *
   * Precomputes spatial derivative coefficients:
   * - hdx_2 = 1/(4*dx²) for centered 2nd-order finite differences
   * - The factor of 4 comes from the stencil: (u[i+1] - 2*u[i] + u[i-1])/(dx²)
   *   when using a centered difference approximation
   *
   * @param opt FDTD configuration options containing grid dimensions,
   *            spacing, and optional model file path
   * @throws std::runtime_error if model file loading fails
   *
   * @post Grid dimensions (nx_, ny_, nz_) and spacing (dx_, dy_, dz_) are set
   * @post Derivative coefficients (hdx_2_, hdy_2_, hdz_2_) are computed
   */
  void InitGrid(const fdtd_options& opt)
  {
    // Grid dimensions
    nx_ = opt.grid.nx;
    ny_ = opt.grid.ny;
    nz_ = opt.grid.nz;

    // Grid sampling
    dx_ = opt.grid.dx;
    dy_ = opt.grid.dy;
    dz_ = opt.grid.dz;

    // Override with model file parameters if provided
    LoadModelInfo(nx_, ny_, nz_, dx_, dy_, dz_, opt.velocity.fileModel);

    // Precompute inverse squared spacing factors for finite differences
    // Factor of 4 accounts for centered difference stencil geometry
    hdx_2_ = 1.0f / (4.0f * dx_ * dx_);
    hdy_2_ = 1.0f / (4.0f * dy_ * dy_);
    hdz_2_ = 1.0f / (4.0f * dz_ * dz_);
  }

  /**
   * @brief Load model dimensions and spacing from binary file
   *
   * Reads grid metadata from a binary model file. The file format expects:
   * - 3 integers: nx, ny, nz (grid dimensions)
   * - 3 floats: dx, dy, dz (spatial sampling in meters)
   *
   * This allows models to define their own grid geometry, which takes
   * precedence over configuration file settings.
   *
   * @param[in,out] nx Grid size in X direction (modified if file loaded)
   * @param[in,out] ny Grid size in Y direction (modified if file loaded)
   * @param[in,out] nz Grid size in Z direction (modified if file loaded)
   * @param[in,out] dx Spatial sampling in X (m) (modified if file loaded)
   * @param[in,out] dy Spatial sampling in Y (m) (modified if file loaded)
   * @param[in,out] dz Spatial sampling in Z (m) (modified if file loaded)
   * @param file_model Path to binary model file (empty string to skip)
   *
   * @throws std::runtime_error if file cannot be opened or is corrupted
   *
   * @note If file_model is empty, function returns immediately without changes
   * @note File format is platform-dependent (native endianness)
   */
  void LoadModelInfo(int& nx, int& ny, int& nz, float& dx, float& dy, float& dz,
                     const std::string& file_model)
  {
    if (file_model.empty())
    {
      return;
    }

    std::ifstream infile(file_model, std::ios::in | std::ios::binary);
    if (!infile)
    {
      throw std::runtime_error("Error opening model file: " + file_model);
    }

    // Read grid dimensions
    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));

    // Read spatial sampling
    infile.read(reinterpret_cast<char*>(&dx), sizeof(float));
    infile.read(reinterpret_cast<char*>(&dy), sizeof(float));
    infile.read(reinterpret_cast<char*>(&dz), sizeof(float));

    if (!infile)
    {
      throw std::runtime_error(
          "Error reading metadata from file: " + file_model +
          " (file may be corrupted or truncated)");
    }

    infile.close();
    printf("Model geometry loaded from file: %s\n", file_model.c_str());
    printf("  Grid: %d x %d x %d, Spacing: %.3f x %.3f x %.3f m\n", nx, ny, nz,
           dx, dy, dz);
  }

  /**
   * @brief Allocate and initialize model property arrays
   *
   * Allocates memory for the P-wave velocity model (vp_) and initializes it
   * with either:
   * - A synthetic two-layer model (default)
   * - Values loaded from a file (if usefilemodel is true, future
   * implementation)
   *
   * The velocity is stored as v² * dt² (squared velocity scaled by time step)
   * for efficient use in the wave equation solver, avoiding repeated
   * multiplications during time stepping.
   *
   * **Memory Layout:** Column-major 3D array, size = nx * ny * nz
   *
   * @param opt FDTD configuration options containing:
   *            - velocity.vmin: Minimum velocity (m/s) for layer 1
   *            - velocity.usefilemodel: Whether to load from file (not yet
   * implemented)
   *
   * @throws std::runtime_error if memory allocation fails
   *
   * @pre InitGrid() must be called first to set grid dimensions
   * @post vp_ is allocated and initialized with nx*ny*nz elements
   *
   * @note Currently only implements synthetic two-layer initialization
   * @note S-wave velocity (vs_) and density (rho_) are not yet implemented
   *
   * @todo Implement file-based velocity model loading
   * @todo Add support for vs_ and rho_ initialization
   */
  void InitModelArrays(const fdtd_options& opt)
  {
    int model_volume = nx_ * ny_ * nz_;

    try
    {
      vp_ = allocateVector<vectorReal>(model_volume, "vp");
    }
    catch (const std::exception& e)
    {
      throw std::runtime_error("Failed to allocate velocity model: " +
                               std::string(e.what()));
    }

    // Initialize velocity model with v² * dt² for wave equation efficiency
    // TODO: Make time step configurable instead of hardcoded
    constexpr float kTimeStep = 0.001f;  // 1 ms time step
    float init_vp_value =
        opt.velocity.vmin * opt.velocity.vmin * kTimeStep * kTimeStep;
    printf("Initializing velocity model: vmin=%f m/s, v²*dt²=%f\n",
           opt.velocity.vmin, init_vp_value);
    InitModel(init_vp_value, opt.velocity.usefilemodel);
  }

  /**
   * @brief Initialize velocity model with synthetic two-layer configuration
   *
   * Creates a simple two-layer model suitable for testing:
   * - Upper half (k < nz/2): base velocity
   * - Lower half (k >= nz/2): doubled velocity (simulates velocity increase
   * with depth)
   *
   * This synthetic model is useful for:
   * - Validating wave propagation code
   * - Testing reflection/transmission at interfaces
   * - Performance benchmarking
   *
   * @param init_vp_value Base velocity value (already scaled as v² * dt²)
   * @param from_file If true, load from file instead (not yet implemented)
   *
   * @throws std::logic_error if from_file=true (feature not implemented)
   *
   * @note Uses OpenMP parallelization with collapse(3) for efficient
   * initialization
   * @note The velocity ratio of 2.0 creates a significant impedance contrast
   *
   * @todo Implement file-based model loading for from_file=true
   */
  void InitModel(float init_vp_value, bool from_file = false)
  {
    if (from_file)
    {
      throw std::logic_error(
          "Loading velocity model from file is not yet implemented. "
          "Set velocity.usefilemodel=false in configuration.");
    }

    // Layer 1: Initialize entire volume with base velocity
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx_; i++)
    {
      for (int j = 0; j < ny_; j++)
      {
        for (int k = 0; k < nz_; k++)
        {
          vp_[Index3D(i, j, k)] = init_vp_value;
        }
      }
    }

    // Layer 2: Overwrite lower half with doubled velocity (simulates depth
    // increase)
    constexpr float kVelocityRatio = 2.0f;  // Interface impedance contrast
    const float layer2_vp = kVelocityRatio * init_vp_value;
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx_; i++)
    {
      for (int j = 0; j < ny_; j++)
      {
        for (int k = nz_ / 2; k < nz_; k++)
        {
          vp_[Index3D(i, j, k)] = layer2_vp;
        }
      }
    }
  }

  /// @name Grid Dimension Accessors
  /// @{

  /**
   * @brief Get number of grid points in X direction
   * @return Grid size in X (slowest varying index)
   */
  [[nodiscard]] int nx() const noexcept { return nx_; }

  /**
   * @brief Get number of grid points in Y direction
   * @return Grid size in Y (middle varying index)
   */
  [[nodiscard]] int ny() const noexcept { return ny_; }

  /**
   * @brief Get number of grid points in Z direction
   * @return Grid size in Z (fastest varying index)
   */
  [[nodiscard]] int nz() const noexcept { return nz_; }

  /**
   * @brief Get spatial sampling in X direction
   * @return Grid spacing in X (meters)
   */
  [[nodiscard]] float dx() const noexcept { return dx_; }

  /**
   * @brief Get spatial sampling in Y direction
   * @return Grid spacing in Y (meters)
   */
  [[nodiscard]] float dy() const noexcept { return dy_; }

  /**
   * @brief Get spatial sampling in Z direction
   * @return Grid spacing in Z (meters)
   */
  [[nodiscard]] float dz() const noexcept { return dz_; }

  /**
   * @brief Get precomputed X-direction finite difference coefficient
   * @return 1/(4*dx²) for centered difference stencil
   */
  [[nodiscard]] float hdx_2() const noexcept { return hdx_2_; }

  /**
   * @brief Get precomputed Y-direction finite difference coefficient
   * @return 1/(4*dy²) for centered difference stencil
   */
  [[nodiscard]] float hdy_2() const noexcept { return hdy_2_; }

  /**
   * @brief Get precomputed Z-direction finite difference coefficient
   * @return 1/(4*dz²) for centered difference stencil
   */
  [[nodiscard]] float hdz_2() const noexcept { return hdz_2_; }
  /// @}

  /// @name Boundary Region Accessors
  /// @{

  /**
   * @brief Get taper layer thickness in X direction
   * @return Number of grid points in X taper/sponge layer
   * @note Currently unused, reserved for future PML/sponge implementation
   */
  [[nodiscard]] int ntaperx() const noexcept { return ntaperx_; }

  /**
   * @brief Get taper layer thickness in Y direction
   * @return Number of grid points in Y taper/sponge layer
   * @note Currently unused, reserved for future PML/sponge implementation
   */
  [[nodiscard]] int ntapery() const noexcept { return ntapery_; }

  /**
   * @brief Get taper layer thickness in Z direction
   * @return Number of grid points in Z taper/sponge layer
   * @note Currently unused, reserved for future PML/sponge implementation
   */
  [[nodiscard]] int ntaperz() const noexcept { return ntaperz_; }

  /**
   * @brief Get damping layer thickness in X direction
   * @return Number of grid points in X damping layer
   * @note Currently unused, reserved for future boundary condition
   * implementation
   */
  [[nodiscard]] int ndampx() const noexcept { return ndampx_; }

  /**
   * @brief Get damping layer thickness in Y direction
   * @return Number of grid points in Y damping layer
   * @note Currently unused, reserved for future boundary condition
   * implementation
   */
  [[nodiscard]] int ndampy() const noexcept { return ndampy_; }

  /**
   * @brief Get damping layer thickness in Z direction
   * @return Number of grid points in Z damping layer
   * @note Currently unused, reserved for future boundary condition
   * implementation
   */
  [[nodiscard]] int ndampz() const noexcept { return ndampz_; }

  // Boundary region limits (reserved for future absorbing boundary
  // implementation)
  [[nodiscard]] int x1() const noexcept
  {
    return x1_;
  }  ///< X-direction boundary limit 1
  [[nodiscard]] int x2() const noexcept
  {
    return x2_;
  }  ///< X-direction boundary limit 2
  [[nodiscard]] int x3() const noexcept
  {
    return x3_;
  }  ///< X-direction boundary limit 3
  [[nodiscard]] int x4() const noexcept
  {
    return x4_;
  }  ///< X-direction boundary limit 4
  [[nodiscard]] int x5() const noexcept
  {
    return x5_;
  }  ///< X-direction boundary limit 5
  [[nodiscard]] int x6() const noexcept
  {
    return x6_;
  }  ///< X-direction boundary limit 6

  [[nodiscard]] int y1() const noexcept
  {
    return y1_;
  }  ///< Y-direction boundary limit 1
  [[nodiscard]] int y2() const noexcept
  {
    return y2_;
  }  ///< Y-direction boundary limit 2
  [[nodiscard]] int y3() const noexcept
  {
    return y3_;
  }  ///< Y-direction boundary limit 3
  [[nodiscard]] int y4() const noexcept
  {
    return y4_;
  }  ///< Y-direction boundary limit 4
  [[nodiscard]] int y5() const noexcept
  {
    return y5_;
  }  ///< Y-direction boundary limit 5
  [[nodiscard]] int y6() const noexcept
  {
    return y6_;
  }  ///< Y-direction boundary limit 6

  [[nodiscard]] int z1() const noexcept
  {
    return z1_;
  }  ///< Z-direction boundary limit 1
  [[nodiscard]] int z2() const noexcept
  {
    return z2_;
  }  ///< Z-direction boundary limit 2
  [[nodiscard]] int z3() const noexcept
  {
    return z3_;
  }  ///< Z-direction boundary limit 3
  [[nodiscard]] int z4() const noexcept
  {
    return z4_;
  }  ///< Z-direction boundary limit 4
  [[nodiscard]] int z5() const noexcept
  {
    return z5_;
  }  ///< Z-direction boundary limit 5
  [[nodiscard]] int z6() const noexcept
  {
    return z6_;
  }  ///< Z-direction boundary limit 6
  /// @}

  /// @name Physical Property Model Accessors
  /// @{

  /**
   * @brief Get read-only access to P-wave velocity model
   * @return Const reference to vp array (v² * dt² values)
   * @note Values are pre-scaled by dt² for computational efficiency
   */
  [[nodiscard]] const vectorReal& vp() const noexcept { return vp_; }

  /**
   * @brief Get read-only access to S-wave velocity model
   * @return Const reference to vs array
   * @note Currently unimplemented, returns empty vector
   */
  [[nodiscard]] const vectorReal& vs() const noexcept { return vs_; }

  /**
   * @brief Get read-only access to density model
   * @return Const reference to rho array
   * @note Currently unimplemented, returns empty vector
   */
  [[nodiscard]] const vectorReal& rho() const noexcept { return rho_; }

  /**
   * @brief Get mutable access to P-wave velocity model
   * @return Non-const reference to vp array for modification
   */
  [[nodiscard]] vectorReal& vp() noexcept { return vp_; }

  /**
   * @brief Get mutable access to S-wave velocity model
   * @return Non-const reference to vs array for modification
   * @note Currently unimplemented, returns empty vector
   */
  [[nodiscard]] vectorReal& vs() noexcept { return vs_; }

  /**
   * @brief Get mutable access to density model
   * @return Non-const reference to rho array for modification
   * @note Currently unimplemented, returns empty vector
   */
  [[nodiscard]] vectorReal& rho() noexcept { return rho_; }
  /// @}

 private:
  /**
   * @brief Convert 3D grid coordinates to linear array index
   *
   * Implements column-major (Fortran-style) indexing where Z is the
   * fastest varying dimension, followed by Y, then X (slowest).
   * This memory layout is cache-friendly for Z-direction sweeps.
   *
   * @param i X-direction grid index [0, nx-1]
   * @param j Y-direction grid index [0, ny-1]
   * @param k Z-direction grid index [0, nz-1]
   * @return Linear array index for accessing 1D storage
   *
   * @note Bounds checking is NOT performed for performance
   * @note Formula: index = nz*ny*i + nz*j + k
   */
  [[nodiscard]] size_t Index3D(int i, int j, int k) const noexcept
  {
    return static_cast<size_t>(nz_) * ny_ * i + nz_ * j + k;
  }

  // Grid dimensions
  int nx_{0};  ///< Number of grid points in X direction
  int ny_{0};  ///< Number of grid points in Y direction
  int nz_{0};  ///< Number of grid points in Z direction

  float dx_{0.0f};  ///< Spatial sampling in X direction (meters)
  float dy_{0.0f};  ///< Spatial sampling in Y direction (meters)
  float dz_{0.0f};  ///< Spatial sampling in Z direction (meters)

  // Sponge / PML parameters (reserved for future absorbing boundary
  // implementation)
  int ntaperx_{0};  ///< Taper layer thickness in X (grid points)
  int ntapery_{0};  ///< Taper layer thickness in Y (grid points)
  int ntaperz_{0};  ///< Taper layer thickness in Z (grid points)
  int ndampx_{0};   ///< Damping layer thickness in X (grid points)
  int ndampy_{0};   ///< Damping layer thickness in Y (grid points)
  int ndampz_{0};   ///< Damping layer thickness in Z (grid points)

  // Precomputed inverse spacing factors for finite difference operators
  float hdx_2_{0.0f};  ///< 1/(4*dx²) for centered 2nd derivative in X
  float hdy_2_{0.0f};  ///< 1/(4*dy²) for centered 2nd derivative in Y
  float hdz_2_{0.0f};  ///< 1/(4*dz²) for centered 2nd derivative in Z

  // Boundary region limits (reserved for future absorbing boundary
  // implementation)
  int x1_{0}, x2_{0}, x3_{0}, x4_{0}, x5_{0},
      x6_{0};  ///< X-direction boundary markers
  int y1_{0}, y2_{0}, y3_{0}, y4_{0}, y5_{0},
      y6_{0};  ///< Y-direction boundary markers
  int z1_{0}, z2_{0}, z3_{0}, z4_{0}, z5_{0},
      z6_{0};  ///< Z-direction boundary markers

  // Physical property models (column-major 3D arrays stored as 1D)
  vectorReal vp_;   ///< P-wave velocity (stored as v² * dt²)
  vectorReal vs_;   ///< S-wave velocity (not yet implemented)
  vectorReal rho_;  ///< Density (not yet implemented)
};

}  // namespace fdgrid
}  // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_
