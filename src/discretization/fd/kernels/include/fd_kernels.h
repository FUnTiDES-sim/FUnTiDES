#ifndef FDTD_KERNEL_HPP
#define FDTD_KERNEL_HPP

#include "data_type.h"
#include "fd_abckernels.h"
#include "fd_macros.h"

namespace fdtd
{
namespace kernel
{
struct FdtdKernels
{
  vectorReal phi;
  vectorReal RHSTerm;
  arrayReal pnGlobal;

  // allocate arrays and vectors
  void initFieldsArrays(int nx, int ny, int nz, int lx, int ly, int lz)
  {
    int extModelVolume = (nx + 2 * lx) * (ny + 2 * ly) * (nz + 2 * lz);
    int modelVolume = nx * ny * nz;
    pnGlobal = allocateArray2D<arrayReal>(extModelVolume, 2, "pnGlobal");
    phi = allocateVector<vectorReal>(modelVolume, "phi");

    for (int i = -lx; i < nx + lx; i++)
    {
      for (int j = -ly; j < ny + ly; j++)
      {
        for (int k = -lz; k < nz + lz; k++)
        {
          pnGlobal(IDX3_l(i, j, k), 0) = 0.000001;
          pnGlobal(IDX3_l(i, j, k), 1) = 0.000001;
        }
      }
    }
  }
  // add RHS term
  int addRHS(const int itSample, int &cb, int const &nx, int const &ny,
             int const &nz, int const &lx, int const &ly, int const &lz,
             int const &xs, int const &ys, int const &zs, vectorReal const &vp,
             vectorReal const &RHSTerm, arrayReal const &pnGlobal) const
  {
    // CREATEVIEWRHS
    LOOP3DHEAD(xs, ys, zs, xs + 1, ys + 1, zs + 1)
    pnGlobal(IDX3_l(i, j, k), cb) += vp[IDX3(i, j, k)] * RHSTerm[itSample];
    LOOP3DEND
    return (0);
  }

  // innerpoints
  int inner3D(int &ca, int &cb, int const &nx, int const &ny, int const &nz,
              int const &lx, int const &ly, int const &lz, const int x3,
              const int x4, const int y3, const int y4, const int z3,
              const int z4, const double coef0, vectorReal const &coefx,
              vectorReal const &coefy, vectorReal const &coefz,
              vectorReal const &vp, arrayReal const &pnGlobal) const
  {
    LOOP3DHEAD(x3, y3, z3, x4, y4, z4)
    float lapx = 0;
    for (int l = 1; l < coefx.extent(0); l++)
    {
      lapx += coefx[l] * (pnGlobal(IDX3_l(i + l, j, k), cb) +
                          pnGlobal(IDX3_l(i - l, j, k), cb));
    }
    float lapy = 0;
    for (int l = 1; l < coefy.extent(0); l++)
    {
      lapy += coefy[l] * (pnGlobal(IDX3_l(i, j + l, k), cb) +
                          pnGlobal(IDX3_l(i, j - l, k), cb));
    }
    float lapz = 0;
    for (int l = 1; l < coefz.extent(0); l++)
    {
      lapz += coefz[l] * (pnGlobal(IDX3_l(i, j, k + l), cb) +
                          pnGlobal(IDX3_l(i, j, k - l), cb));
    }
    pnGlobal(IDX3_l(i, j, k), ca) =
        2. * pnGlobal(IDX3_l(i, j, k), cb) - pnGlobal(IDX3_l(i, j, k), ca) +
        vp[IDX3(i, j, k)] *
            (coef0 * pnGlobal(IDX3_l(i, j, k), cb) + lapx + lapy + lapz);
    LOOP3DEND
    return 0;
  }

  // pml update
  int pml3D(int &ca, int &cb, int const &nx, int const &ny, int const &nz,
            int const &lx, int const &ly, int const &lz, const int x3,
            const int x4, const int y3, const int y4, const int z3,
            const int z4, const double coef0, const float &hdx_2,
            const float &hdy_2, const float &hdz_2, vectorReal const &coefx,
            vectorReal const &coefy, vectorReal const &coefz,
            vectorReal const &vp, vectorReal const &eta, vectorReal &phi,
            arrayReal &pnGlobal) const
  {
    LOOP3DHEAD(x3, y3, z3, x4, y4, z4)
    float lapx = 0;
    for (int l = 1; l < coefx.extent(0); l++)
    {
      lapx += coefx[l] * (pnGlobal(IDX3_l(i + l, j, k), cb) +
                          pnGlobal(IDX3_l(i - l, j, k), cb));
    }
    float lapy = 0;
    for (int l = 1; l < coefy.extent(0); l++)
    {
      lapy += coefy[l] * (pnGlobal(IDX3_l(i, j + l, k), cb) +
                          pnGlobal(IDX3_l(i, j - l, k), cb));
    }
    float lapz = 0;
    for (int l = 1; l < coefz.extent(0); l++)
    {
      lapz += coefz[l] * (pnGlobal(IDX3_l(i, j, k + l), cb) +
                          pnGlobal(IDX3_l(i, j, k - l), cb));
    }

    float lap = coef0 * pnGlobal(IDX3_l(i, j, k), cb) + lapx + lapy + lapz;

    pnGlobal(IDX3_l(i, j, k), ca) =
        ((2. - eta[IDX3_eta1(i, j, k)] * eta[IDX3_eta1(i, j, k)] +
          2. * eta[IDX3_eta1(i, j, k)]) *
             pnGlobal(IDX3_l(i, j, k), cb) -
         pnGlobal(IDX3_l(i, j, k), ca) +
         vp[IDX3(i, j, k)] * (phi[IDX3(i, j, k)] + lap)) /
        (1. + 2. * eta[IDX3_eta1(i, j, k)]);

    phi[IDX3(i, j, k)] =
        (phi[IDX3(i, j, k)] -
         ((eta[IDX3_eta1(i + 1, j, k)] - eta[IDX3_eta1(i - 1, j, k)]) *
              (pnGlobal(IDX3_l(i + 1, j, k), cb) -
               pnGlobal(IDX3_l(i - 1, j, k), cb)) *
              hdx_2 +
          (eta[IDX3_eta1(i, j + 1, k)] - eta[IDX3_eta1(i, j - 1, k)]) *
              (pnGlobal(IDX3_l(i, j + 1, k), cb) -
               pnGlobal(IDX3_l(i, j - 1, k), cb)) *
              hdy_2 +
          (eta[IDX3_eta1(i, j, k + 1)] - eta[IDX3_eta1(i, j, k - 1)]) *
              (pnGlobal(IDX3_l(i, j, k + 1), cb) -
               pnGlobal(IDX3_l(i, j, k - 1), cb)) *
              hdz_2)) /
        (1. + eta[IDX3_eta1(i, j, k)]);
    LOOP3DEND
    return (0);
  }

  // apply sponge boundary to wavefield
  int applySponge(int &ca, int &cb, int const &nx, int const &ny, int const &nz,
                  int const &lx, int const &ly, int const &lz, const int x3,
                  const int x4, const int y3, const int y4, const int z3,
                  const int z4, vectorReal const &spongeArray,
                  arrayReal const &pnGlobal) const
  {
    // CREATEVIEWSPONGE
    LOOP3DHEAD(x3, y3, z3, x4, y4, z4)
    pnGlobal(IDX3_l(i, j, k), ca) =
        pnGlobal(IDX3_l(i, j, k), ca);  // * spongeArray(IDX3(i, j, k));
    pnGlobal(IDX3_l(i, j, k), cb) =
        pnGlobal(IDX3_l(i, j, k), cb);  // * spongeArray(IDX3(i, j, k));
    LOOP3DEND
    return 0;
  }
};

}  // namespace kernel
}  // namespace fdtd
#endif  // FDTD_KERNELS_H
