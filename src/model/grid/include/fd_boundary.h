/**
 * @file fd_boundary.h
 * @brief Absorbing boundary layer management
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_BOUNDARY_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_BOUNDARY_H_

namespace model
{
namespace fdgrid
{

/**
 * @brief Absorbing boundary layer configuration
 *
 * Manages PML (Perfectly Matched Layer) and sponge boundary parameters.
 * Stores the geometry of absorbing regions around the computational domain
 * to prevent artificial reflections from domain edges.
 */
class BoundaryLayers
{
 public:
  BoundaryLayers() = default;

  /**
   * @brief Configure boundary layers from options
   * @param opt Configuration containing boundary parameters
   */
  void Initialize(const fdtd::options::FdtdOptions& opt)
  {
    // TODO: Extract from opt.boundary or similar
    // For now, set to reasonable defaults
    float lambdamax_ = opt.velocity.vmin /
                       (2.5 * opt.source.f0);  // Maximum wavelength in model
    ntaperx_ = ntapery_ = ntaperz_ = opt.boundary.pml_size;
    if (opt.boundary.use_pml)
    {
      ndampx_ = ntaperx_ * lambdamax_ / opt.grid.dx;
      ndampy_ = ntapery_ * lambdamax_ / opt.grid.dy;
      ndampz_ = ntaperz_ * lambdamax_ / opt.grid.dz;
    }
    if (opt.boundary.use_sponge)
    {
      ndampx_ = 0;
      ndampy_ = 0;
      ndampz_ = 0;
    }

    ComputeBoundaryRegions(opt.boundary.use_sponge, opt.boundary.use_pml,
                           opt.grid.nx, opt.grid.ny, opt.grid.nz);
  }

  /**
   * @brief Compute boundary region limits
   * @param nx Grid size X
   * @param ny Grid size Y
   * @param nz Grid size Z
   */
  void ComputeBoundaryRegions(bool sponge, bool pml, int nx, int ny, int nz)
  {
    // Example: divide domain into 6 regions per axis
    // x1-x2: left boundary, x3-x4: interior, x5-x6: right boundary
    x1_ = 0;
    x2_ = ndampx_;
    x3_ = ndampx_;
    x4_ = nx - ndampx_;
    x5_ = nx - ndampx_;
    x6_ = nx;

    y1_ = 0;
    y2_ = ndampy_;
    y3_ = ndampy_;
    y4_ = ny - ndampy_;
    y5_ = ny - ndampy_;
    y6_ = ny;

    z1_ = 0;
    z2_ = ndampz_;
    z3_ = ndampz_;
    z4_ = nz - ndampz_;
    z5_ = nz - ndampz_;
    z6_ = nz;
  }

  // Accessors
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

 private:
  int ntaperx_{0}, ntapery_{0}, ntaperz_{0};
  int ndampx_{0}, ndampy_{0}, ndampz_{0};
  int x1_{0}, x2_{0}, x3_{0}, x4_{0}, x5_{0}, x6_{0};
  int y1_{0}, y2_{0}, y3_{0}, y4_{0}, y5_{0}, y6_{0};
  int z1_{0}, z2_{0}, z3_{0}, z4_{0}, z5_{0}, z6_{0};
};

}  // namespace fdgrid
}  // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_BOUNDARY_H_
