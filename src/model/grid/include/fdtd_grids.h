/**
 * @file fdtd_grids.h
 * @brief FDTD grid and model management
 *
 * @license
 * Copyright (C) 2025.
 * This software is governed by the CeCILL license under French law.
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_

#include <data_type.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "fdtd_macros.h"
#include <fdtd_options.h>

namespace model {
namespace fdgrid {
/**
 * @brief FDTD computational grid and model container
 *
 * Manages grid geometry, boundary conditions, and physical properties
 * for finite-difference time-domain simulations.
 */
class FdtdGrids
{
 public:
  FdtdGrids() = default;

  /**
   * @brief Initialize grid geometry from options
   * @param opt FDTD configuration options
   * @throws std::runtime_error if file loading fails
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

    // Get model info if from file
    LoadModelInfo(nx_, ny_, nz_, dx_, dy_, dz_, opt.velocity.fileModel);

    // Precompute inverse squared spacing factors
    hdx_2_ = 1.0f / (4.0f * dx_ * dx_);
    hdy_2_ = 1.0f / (4.0f * dy_ * dy_);
    hdz_2_ = 1.0f / (4.0f * dz_ * dz_);
  }

  /**
   * @brief Load model dimensions from binary file
   * @param nx Reference to grid size X
   * @param ny Reference to grid size Y
   * @param nz Reference to grid size Z
   * @param dx Reference to spacing X
   * @param dy Reference to spacing Y
   * @param dz Reference to spacing Z
   * @param file_model Path to model file
   * @throws std::runtime_error if file cannot be opened or read
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
      throw std::runtime_error("Error opening file: " + file_model);
    }

    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&dx), sizeof(float));
    infile.read(reinterpret_cast<char*>(&dy), sizeof(float));
    infile.read(reinterpret_cast<char*>(&dz), sizeof(float));

    if (!infile)
    {
      throw std::runtime_error("Error reading data from file: " + file_model);
    }

    infile.close();
    printf("Model read from file %s\n", file_model.c_str());
  }

  /**
   * @brief Initialize model arrays and velocity field
   * @param opt FDTD configuration options
   * @throws std::runtime_error if allocation fails
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

    // Initialize velocity model
    // Using squared velocity scaled by time step for wave equation
    constexpr float kTimeStep = 0.001f;
    float init_vp_value =
        opt.velocity.vmin * opt.velocity.vmin * kTimeStep * kTimeStep;
    printf("init_vp_value=%f\n", init_vp_value);
    InitModel(init_vp_value, opt.velocity.usefilemodel);
  }

  /**
   * @brief Initialize two-layer velocity model
   * @param init_vp_value Base velocity value (v² * dt²)
   * @param from_file Whether to load from file (not yet implemented)
   */
  void InitModel(float init_vp_value, bool from_file = false)
  {
    if (!from_file)
    {
      // Layer 1: Initialize upper half with base velocity
#pragma omp parallel for collapse(3)
      for (int i = 0; i < nx_; i++)
      {
        for (int j = 0; j < ny_; j++)
        {
          for (int k = 0; k < nz_; k++)
          {
            vp_[nz_ * ny_ * i + nz_ * j + k] = init_vp_value;
          }
        }
      }

      // Layer 2: Initialize lower half with doubled velocity
      constexpr float kVelocityRatio = 2.0f;
#pragma omp parallel for collapse(3)
      for (int i = 0; i < nx_; i++)
      {
        for (int j = 0; j < ny_; j++)
        {
          for (int k = nz_ / 2; k < nz_; k++)
          {
            vp_[nz_ * ny_ * i + nz_ * j + k] = kVelocityRatio * init_vp_value;
          }
        }
      }
    }
    // TODO: Implement loading velocity model from file
  }

  /// @name Grid Dimension Accessors
  /// @{
  [[nodiscard]] int nx() const noexcept { return nx_; }
  [[nodiscard]] int ny() const noexcept { return ny_; }
  [[nodiscard]] int nz() const noexcept { return nz_; }
  [[nodiscard]] float dx() const noexcept { return dx_; }
  [[nodiscard]] float dy() const noexcept { return dy_; }
  [[nodiscard]] float dz() const noexcept { return dz_; }
  [[nodiscard]] float hdx_2() const noexcept { return hdx_2_; }
  [[nodiscard]] float hdy_2() const noexcept { return hdy_2_; }
  [[nodiscard]] float hdz_2() const noexcept { return hdz_2_; }
  /// @}

  /// @name Boundary Accessors
  /// @{
  [[nodiscard]] int ntaperx() const noexcept { return ntaperx_; }
  [[nodiscard]] int ntapery() const noexcept { return ntapery_; }
  [[nodiscard]] int ntaperz() const noexcept { return ntaperz_; }
  [[nodiscard]] int ndampx() const noexcept { return ndampx_; }
  [[nodiscard]] int ndampy() const noexcept { return ndampy_; }
  [[nodiscard]] int ndampz() const noexcept { return ndampz_; }

  [[nodiscard]] int x1() const noexcept { return x1_; }
  [[nodiscard]] int x2() const noexcept { return x2_; }
  [[nodiscard]] int x3() const noexcept { return x3_; }
  [[nodiscard]] int x4() const noexcept { return x4_; }
  [[nodiscard]] int x5() const noexcept { return x5_; }
  [[nodiscard]] int x6() const noexcept { return x6_; }

  [[nodiscard]] int y1() const noexcept { return y1_; }
  [[nodiscard]] int y2() const noexcept { return y2_; }
  [[nodiscard]] int y3() const noexcept { return y3_; }
  [[nodiscard]] int y4() const noexcept { return y4_; }
  [[nodiscard]] int y5() const noexcept { return y5_; }
  [[nodiscard]] int y6() const noexcept { return y6_; }

  [[nodiscard]] int z1() const noexcept { return z1_; }
  [[nodiscard]] int z2() const noexcept { return z2_; }
  [[nodiscard]] int z3() const noexcept { return z3_; }
  [[nodiscard]] int z4() const noexcept { return z4_; }
  [[nodiscard]] int z5() const noexcept { return z5_; }
  [[nodiscard]] int z6() const noexcept { return z6_; }
  /// @}

  /// @name Model Accessors
  /// @{
  [[nodiscard]] const vectorReal& vp() const noexcept { return vp_; }
  [[nodiscard]] const vectorReal& vs() const noexcept { return vs_; }
  [[nodiscard]] const vectorReal& rho() const noexcept { return rho_; }

  [[nodiscard]] vectorReal& vp() noexcept { return vp_; }
  [[nodiscard]] vectorReal& vs() noexcept { return vs_; }
  [[nodiscard]] vectorReal& rho() noexcept { return rho_; }
  /// @}

 private:
  // Grid dimensions
  int nx_{0};
  int ny_{0};
  int nz_{0};
  float dx_{0.0f};
  float dy_{0.0f};
  float dz_{0.0f};

  // Sponge / PML parameters
  int ntaperx_{0};
  int ntapery_{0};
  int ntaperz_{0};
  int ndampx_{0};
  int ndampy_{0};
  int ndampz_{0};

  // Precomputed inverse spacing factors
  float hdx_2_{0.0f};
  float hdy_2_{0.0f};
  float hdz_2_{0.0f};

  // Boundary limits
  int x1_{0}, x2_{0}, x3_{0}, x4_{0}, x5_{0}, x6_{0};
  int y1_{0}, y2_{0}, y3_{0}, y4_{0}, y5_{0}, y6_{0};
  int z1_{0}, z2_{0}, z3_{0}, z4_{0}, z5_{0}, z6_{0};

  // Physical properties
  vectorReal vp_;   ///< P-wave velocity
  vectorReal vs_;   ///< S-wave velocity
  vectorReal rho_;  ///< Density
};


} // namespace fdgrid
} // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_
