/**
 * @file fdtd_grid_geometry.h
 * @brief Grid geometry and indexing utilities
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_GRID_GEOMETRY_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_GRID_GEOMETRY_H_

#include <fd_options.h>

#include <cstddef>

namespace model
{
namespace fdgrid
{

/**
 * @brief 3D grid geometry and spatial discretization
 *
 * Manages the computational domain geometry including:
 * - Grid dimensions (nx, ny, nz)
 * - Spatial sampling (dx, dy, dz)
 * - 3D-to-1D index conversion
 * - Precomputed finite difference coefficients
 *
 * This is a lightweight, copyable class focused solely on grid structure.
 * It contains no large data arrays and can be safely passed by value.
 *
 * **Thread-safe:** All operations are const after construction
 */
class GridGeometry
{
 public:
  GridGeometry() = default;

  /**
   * @brief Construct grid geometry from dimensions and spacing
   * @param nx Number of grid points in X direction
   * @param ny Number of grid points in Y direction
   * @param nz Number of grid points in Z direction
   * @param dx Spatial sampling in X (meters)
   * @param dy Spatial sampling in Y (meters)
   * @param dz Spatial sampling in Z (meters)
   */
  GridGeometry(int nx, int ny, int nz, float dx, float dy, float dz)
      : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz)
  {
    ComputeDerivativeCoefficients();
  }

  /**
   * @brief Construct from FDTD options
   * @param opt Configuration options
   */
  explicit GridGeometry(const fdtd::options::FdtdOptions& opt)
      : GridGeometry(opt.grid.nx, opt.grid.ny, opt.grid.nz, opt.grid.dx,
                     opt.grid.dy, opt.grid.dz)
  {
  }

  /**
   * @brief Convert 3D grid coordinates to linear array index
   * @param i X-direction index [0, nx-1]
   * @param j Y-direction index [0, ny-1]
   * @param k Z-direction index [0, nz-1]
   * @return Linear index for column-major storage
   */
  [[nodiscard]] size_t Index3D(int i, int j, int k) const noexcept
  {
    return static_cast<size_t>(nz_) * ny_ * i + nz_ * j + k;
  }

  /**
   * @brief Get total number of grid points
   * @return nx * ny * nz
   */
  [[nodiscard]] size_t TotalPoints() const noexcept
  {
    return static_cast<size_t>(nx_) * ny_ * nz_;
  }

  /**
   * @brief Check if indices are within grid bounds
   * @param i X index
   * @param j Y index
   * @param k Z index
   * @return true if all indices are valid
   */
  [[nodiscard]] bool IsValidIndex(int i, int j, int k) const noexcept
  {
    return i >= 0 && i < nx_ && j >= 0 && j < ny_ && k >= 0 && k < nz_;
  }

  // Accessors
  [[nodiscard]] int nx() const noexcept { return nx_; }
  [[nodiscard]] int ny() const noexcept { return ny_; }
  [[nodiscard]] int nz() const noexcept { return nz_; }
  [[nodiscard]] float dx() const noexcept { return dx_; }
  [[nodiscard]] float dy() const noexcept { return dy_; }
  [[nodiscard]] float dz() const noexcept { return dz_; }
  [[nodiscard]] float hdx_2() const noexcept { return hdx_2_; }
  [[nodiscard]] float hdy_2() const noexcept { return hdy_2_; }
  [[nodiscard]] float hdz_2() const noexcept { return hdz_2_; }

 private:
  void ComputeDerivativeCoefficients()
  {
    hdx_2_ = 1.0f / (4.0f * dx_ * dx_);
    hdy_2_ = 1.0f / (4.0f * dy_ * dy_);
    hdz_2_ = 1.0f / (4.0f * dz_ * dz_);
  }

  int nx_{0}, ny_{0}, nz_{0};
  float dx_{0.0f}, dy_{0.0f}, dz_{0.0f};
  float hdx_2_{0.0f}, hdy_2_{0.0f}, hdz_2_{0.0f};
};

}  // namespace fdgrid
}  // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_GRID_GEOMETRY_H_
