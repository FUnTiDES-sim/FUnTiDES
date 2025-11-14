/**
 * @file fdtd_velocity_model.h
 * @brief Velocity model storage and initialization
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_VELOCITY_MODEL_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_VELOCITY_MODEL_H_

#include <data_type.h>
#include <fd_options.h>

#include "fd_grid_geometry.h"
// #include "fdtd_macros.h"

namespace model
{
namespace fdgrid
{

/**
 * @brief Velocity model storage and initialization
 *
 * Manages P-wave and S-wave velocity fields and density.
 * Supports both synthetic model generation and file-based loading.
 */
class VelocityModel
{
 public:
  explicit VelocityModel(const GridGeometry& geom) : geom_(geom) {}

  /**
   * @brief Initialize model from configuration options
   * @param opt FDTD options
   * @param time_step Time step for wave equation (dt)
   * @throws std::runtime_error if allocation fails
   */
  void Initialize(const fdtd::options::FdtdOptions& opt, float time_step)
  {
    size_t model_volume = geom_.TotalPoints();

    try
    {
      vp_ = allocateVector<vectorReal>(model_volume, "vp");
    }
    catch (const std::exception& e)
    {
      throw std::runtime_error("Failed to allocate velocity model: " +
                               std::string(e.what()));
    }

    // Initialize with scaled velocity for wave equation
    float init_vp_value =
        opt.velocity.vmin * opt.velocity.vmin * time_step * time_step;

    if (opt.velocity.use_file_model && !opt.velocity.file_model.empty())
    {
      LoadFromFile(opt.velocity.file_model);
    }
    else
    {
      InitializeSyntheticModel(init_vp_value);
    }
  }

  /**
   * @brief Create two-layer synthetic model for testing
   * @param init_vp_value Base velocity (v² * dt²)
   */
  void InitializeSyntheticModel(float init_vp_value)
  {
    const int nx = geom_.nx();
    const int ny = geom_.ny();
    const int nz = geom_.nz();

    // Layer 1: Entire volume
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int k = 0; k < nz; k++)
        {
          vp_[geom_.Index3D(i, j, k)] = init_vp_value;
        }
      }
    }

    // Layer 2: Lower half with doubled velocity
    constexpr float kVelocityRatio = 2.0f;
    const float layer2_vp = kVelocityRatio * init_vp_value;
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int k = nz / 2; k < nz; k++)
        {
          vp_[geom_.Index3D(i, j, k)] = layer2_vp;
        }
      }
    }

    printf("Initialized two-layer synthetic model\n");
  }

  /**
   * @brief Load velocity model from file
   * @param file_path Path to model file
   * @throws std::runtime_error if loading fails
   */
  void LoadFromFile(const std::string& file_path)
  {
    ModelIO::ReadVelocityModel(file_path, geom_, vp_);
  }

  // Accessors
  [[nodiscard]] const vectorReal& vp() const noexcept { return vp_; }
  [[nodiscard]] const vectorReal& vs() const noexcept { return vs_; }
  [[nodiscard]] const vectorReal& rho() const noexcept { return rho_; }
  [[nodiscard]] vectorReal& vp() noexcept { return vp_; }
  [[nodiscard]] vectorReal& vs() noexcept { return vs_; }
  [[nodiscard]] vectorReal& rho() noexcept { return rho_; }

 private:
  GridGeometry geom_;
  vectorReal vp_;
  vectorReal vs_;
  vectorReal rho_;
};

}  // namespace fdgrid
}  // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_VELOCITY_MODEL_H_
