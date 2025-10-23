#ifndef FDTD_STENCILS_H
#define FDTD_STENCILS_H

#include <data_type.h>

#include <cmath>

#include "fd_options.h"

using namespace std;

namespace fdtd
{
namespace stencils
{
struct FdtdStencils
{
  int ncoefsX{0}, ncoefsY{0}, ncoefsZ{0};
  int lx{0}, ly{0}, lz{0};  // half stencil lenght
  double coef0{0.0};

  vectorReal coefx;
  vectorReal coefy;
  vectorReal coefz;

  void init_coef(int L, float dx, vectorReal &coef)
  {
    float dx2 = dx * dx;
    switch (L)
    {
      case 1:
        coef[0] = -2.f / dx2;
        coef[1] = 1.f / dx2;
        break;
      case 2:
        coef[0] = -5.f / 2.f / dx2;
        coef[1] = 4.f / 3.f / dx2;
        coef[2] = -1.f / 12.f / dx2;
        break;
      case 3:
        coef[0] = -49.f / 18.f / dx2;
        coef[1] = 3.f / 2.f / dx2;
        coef[2] = -3. / 20. / dx2;
        coef[3] = 1. / 90. / dx2;
        break;
      case 4:
        coef[0] = -205.f / 72.f / dx2;
        coef[1] = 8.f / 5.f / dx2;
        coef[2] = -1.f / 5.f / dx2;
        coef[3] = 8.f / 315.f / dx2;
        coef[4] = -1.f / 560.f / dx2;
        break;
      case 5:
        coef[0] = -5269.f / 1800.f / dx2;
        coef[1] = 5.f / 3.f / dx2;
        coef[2] = -5.f / 21.f / dx2;
        coef[3] = 5.f / 126.f / dx2;
        coef[4] = -5.f / 1008.f / dx2;
        coef[5] = 1.f / 3150.f / dx2;
        break;
      case 6:
        coef[0] = -5369.f / 1800.f / dx2;
        coef[1] = 12.f / 7.f / dx2;
        coef[2] = -15.f / 56.f / dx2;
        coef[3] = 10.f / 189.f / dx2;
        coef[4] = -1.f / 112.f / dx2;
        coef[5] = 2.f / 1925.f / dx2;
        coef[6] = 1.f / 1663.f / dx2;
        break;
    }
  }

  void initStencilsCoefficients(fdtd::options::FdtdOptions &m_opt, float dx,
                                float dy, float dz)
  {
    // half stencil lenght
    lx = m_opt.stencil.lx;
    ly = m_opt.stencil.ly;
    lz = m_opt.stencil.lz;

    // number of coefficients
    ncoefsX = lx + 1;
    ncoefsY = ly + 1;
    ncoefsZ = lz + 1;

    // allocate coef
    coefx = allocateVector<vectorReal>(ncoefsX, "coefx");
    coefy = allocateVector<vectorReal>(ncoefsY, "coefy");
    coefz = allocateVector<vectorReal>(ncoefsZ, "coefz");
    // init coef
    init_coef(lx, dx, coefx);
    init_coef(ly, dy, coefy);
    init_coef(lz, dz, coefz);

    // compute coef0 and all coefs such that sum(coef)=0
    float tmpX = 0;
    for (int i = 1; i < ncoefsX; i++)
    {
      tmpX += coefx[i];
    }
    float tmpY = 0;
    for (int i = 1; i < ncoefsY; i++)
    {
      tmpY += coefy[i];
    }
    float tmpZ = 0;
    for (int i = 1; i < ncoefsZ; i++)
    {
      tmpZ += coefz[i];
    }
    coef0 = -2. * (tmpX + tmpY + tmpZ);
  }

  // compute stable time step
  float compute_dt_sch(float vmax)
  {
    float ftmp = 0.;
    float cfl = 0.8;
    ftmp += fabsf(coefx[0]) + fabsf(coefy[0]) + fabsf(coefz[0]);
    for (int i = 1; i < coefx.extent(0); i++)
    {
      ftmp += 2.f * fabsf(coefx[i]);
    }
    for (int i = 1; i < coefy.extent(0); i++)
    {
      ftmp += 2.f * fabsf(coefy[i]);
    }
    for (int i = 1; i < coefz.extent(0); i++)
    {
      ftmp += 2.f * fabsf(coefz[i]);
    }
    return 2 * cfl / (sqrtf(ftmp) * vmax);
  }
};

}  // namespace stencils
}  // namespace fdtd

#endif  // FDTD_STENCILS_H
