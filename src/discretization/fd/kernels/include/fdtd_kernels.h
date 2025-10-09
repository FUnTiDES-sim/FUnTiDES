#ifndef FDTD_KERNEL_HPP
#define FDTD_KERNEL_HPP

#include "data_type.h"
#include "fdtd_macros.h"

struct FdtdKernels
{
  vectorReal spongeArray;
  vectorReal RHSTerm;
  vectorReal phi;
  vectorReal eta;
  arrayReal pnGlobal;

  // allocate arrays and vectors
  void initFieldsArrays(int nx, int ny, int nz, int lx, int ly, int lz)
  {
    int modelVolume = nx * ny * nz;
    int extModelVolume = (nx + 2 * lx) * (ny + 2 * ly) * (nz + 2 * lz);

    spongeArray = allocateVector<vectorReal>(modelVolume, "spongeArray");
    pnGlobal = allocateArray2D<arrayReal>(extModelVolume, 2, "pnGlobal");
#pragma omp parallel for collapse(3)
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
    // CREATEVIEWINNER
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

  // sponge boundary
  // define sponge boundary
  void defineSpongeBoundary(int nx, int ny, int nz)
  {
    const int L = 20;
    const float alpha = -0.00015;

    // compute sponge boundary terms
    // intailize to 1
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int i = 0; i < nx; i++)
        {
          spongeArray(IDX3(i, j, k)) = 1;
        }
        for (int i = 0; i < nx; i++)
        {
          spongeArray(IDX3(i, j, k)) = 1;
        }
      }
    }

    // X boundary
    for (int k = L; k < nz - L; k++)
    {
      for (int j = L; j < ny - L; j++)
      {
        for (int i = 0; i < L; i++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-i)*dx,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - i), 2));
        }
        for (int i = nx - L; i < nx; i++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-(nx-i))*dx,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - (nx - i)), 2));
        }
      }
    }

    // Y boundary
    for (int k = L; k < nz - L; k++)
    {
      for (int i = L; i < nx - L; i++)
      {
        for (int j = 0; j < L; j++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-j)*dy,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - j), 2));
        }
        for (int j = ny - L; j < ny; j++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-(ny-j))*dy,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - (ny - j)), 2));
        }
      }
    }

    // Z boundary
    for (int j = L; j < ny - L; j++)
    {
      for (int i = L; i < nx - L; i++)
      {
        for (int k = 0; k < L; k++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-k))*dz,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - k), 2));
        }
        for (int k = nz - L; k < nz; k++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-(nz-k))*dz,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - (nz - k)), 2));
        }
      }
    }
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
        pnGlobal(IDX3_l(i, j, k), ca) * spongeArray(IDX3(i, j, k));
    pnGlobal(IDX3_l(i, j, k), cb) =
        pnGlobal(IDX3_l(i, j, k), cb) * spongeArray(IDX3(i, j, k));
    LOOP3DEND
    return 0;
  }
};
#endif  // FDTD_KERNELS_H
