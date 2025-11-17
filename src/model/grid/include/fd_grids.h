/**
 * @file fdtd_grids.h
 * @brief FDTD grid facade (backward compatibility wrapper)
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_

#include <fd_options.h>

#include "fd_boundary.h"
#include "fd_grid_geometry.h"
#include "fd_model_io.h"
#include "fd_velocity_model.h"

namespace model
{
namespace fdgrid
{

/**
 * @brief FDTD grid facade maintaining backward compatibility
 *
 * This class provides the same interface as the original monolithic
 * FdtdGrids class, but internally delegates to specialized components:
 * - GridGeometry: dimensions and indexing
 * - BoundaryLayers: PML/sponge configuration
 * - VelocityModel: physical property storage
 * - ModelIO: file I/O operations
 *
 * **Migration Path:**
 * Existing code continues to work unchanged. New code should use
 * the individual components directly for better encapsulation.
 */
class FdtdGrids
{
 public:
  FdtdGrids() = default;

  /**
   * @brief Initialize grid geometry (backward compatible)
   * @param opt FDTD configuration options
   */
  void InitGrid(const fdtd::options::FdtdOptions& opt)
  {
    // Load geometry from file if specified
    if (!opt.velocity.file_model.empty())
    {
      geom_ = ModelIO::ReadGeometry(opt.velocity.file_model);
    }
    else
    {
      geom_ = GridGeometry(opt);
    }

    // Initialize boundary layers
    boundary_.Initialize(opt);
  }

  /**
   * @brief Initialize model arrays (backward compatible)
   * @param opt FDTD configuration options
   */
  void InitModelArrays(const fdtd::options::FdtdOptions& opt)
  {
    constexpr float kTimeStep = 0.001f;  // TODO: Make configurable
    model_ = std::make_unique<VelocityModel>(geom_);
    model_->Initialize(opt, kTimeStep);
  }

  /**
   * @brief Legacy method - no longer used
   * @deprecated Use ModelIO::ReadGeometry() instead
   */
  void LoadModelInfo(int& nx, int& ny, int& nz, float& dx, float& dy, float& dz,
                     const std::string& file_model)
  {
    if (file_model.empty()) return;
    auto geom = ModelIO::ReadGeometry(file_model);
    nx = geom.nx();
    ny = geom.ny();
    nz = geom.nz();
    dx = geom.dx();
    dy = geom.dy();
    dz = geom.dz();
  }

  /**
   * @brief Legacy method - use VelocityModel::Initialize() instead
   * @deprecated
   */
  void InitModel(float init_vp_value, bool from_file = false)
  {
    if (from_file)
    {
      throw std::logic_error("File loading not supported in legacy method");
    }
    model_->InitializeSyntheticModel(init_vp_value);
  }

  // Geometry accessors (delegate to GridGeometry)
  [[nodiscard]] int nx() const noexcept { return geom_.nx(); }
  [[nodiscard]] int ny() const noexcept { return geom_.ny(); }
  [[nodiscard]] int nz() const noexcept { return geom_.nz(); }
  [[nodiscard]] float dx() const noexcept { return geom_.dx(); }
  [[nodiscard]] float dy() const noexcept { return geom_.dy(); }
  [[nodiscard]] float dz() const noexcept { return geom_.dz(); }
  [[nodiscard]] float hdx_2() const noexcept { return geom_.hdx_2(); }
  [[nodiscard]] float hdy_2() const noexcept { return geom_.hdy_2(); }
  [[nodiscard]] float hdz_2() const noexcept { return geom_.hdz_2(); }

  // Boundary accessors (delegate to BoundaryLayers)
  [[nodiscard]] int ntaperx() const noexcept { return boundary_.ntaperx(); }
  [[nodiscard]] int ntapery() const noexcept { return boundary_.ntapery(); }
  [[nodiscard]] int ntaperz() const noexcept { return boundary_.ntaperz(); }
  [[nodiscard]] int ndampx() const noexcept { return boundary_.ndampx(); }
  [[nodiscard]] int ndampy() const noexcept { return boundary_.ndampy(); }
  [[nodiscard]] int ndampz() const noexcept { return boundary_.ndampz(); }

  [[nodiscard]] int x1() const noexcept { return boundary_.x1(); }
  [[nodiscard]] int x2() const noexcept { return boundary_.x2(); }
  [[nodiscard]] int x3() const noexcept { return boundary_.x3(); }
  [[nodiscard]] int x4() const noexcept { return boundary_.x4(); }
  [[nodiscard]] int x5() const noexcept { return boundary_.x5(); }
  [[nodiscard]] int x6() const noexcept { return boundary_.x6(); }

  [[nodiscard]] int y1() const noexcept { return boundary_.y1(); }
  [[nodiscard]] int y2() const noexcept { return boundary_.y2(); }
  [[nodiscard]] int y3() const noexcept { return boundary_.y3(); }
  [[nodiscard]] int y4() const noexcept { return boundary_.y4(); }
  [[nodiscard]] int y5() const noexcept { return boundary_.y5(); }
  [[nodiscard]] int y6() const noexcept { return boundary_.y6(); }

  [[nodiscard]] int z1() const noexcept { return boundary_.z1(); }
  [[nodiscard]] int z2() const noexcept { return boundary_.z2(); }
  [[nodiscard]] int z3() const noexcept { return boundary_.z3(); }
  [[nodiscard]] int z4() const noexcept { return boundary_.z4(); }
  [[nodiscard]] int z5() const noexcept { return boundary_.z5(); }
  [[nodiscard]] int z6() const noexcept { return boundary_.z6(); }

  // Model accessors (delegate to VelocityModel)
  [[nodiscard]] const vectorReal& vp() const noexcept { return model_->vp(); }
  [[nodiscard]] const vectorReal& vs() const noexcept { return model_->vs(); }
  [[nodiscard]] const vectorReal& rho() const noexcept { return model_->rho(); }
  [[nodiscard]] vectorReal& vp() noexcept { return model_->vp(); }
  [[nodiscard]] vectorReal& vs() noexcept { return model_->vs(); }
  [[nodiscard]] vectorReal& rho() noexcept { return model_->rho(); }

  // Direct access to components for new code
  [[nodiscard]] const GridGeometry& geometry() const noexcept { return geom_; }
  [[nodiscard]] const BoundaryLayers& boundary() const noexcept
  {
    return boundary_;
  }
  [[nodiscard]] const VelocityModel& model() const noexcept { return *model_; }
  [[nodiscard]] VelocityModel& model() noexcept { return *model_; }

 private:
  GridGeometry geom_;
  BoundaryLayers boundary_;
  std::unique_ptr<VelocityModel> model_;
};

}  // namespace fdgrid
}  // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_GRIDS_H_
