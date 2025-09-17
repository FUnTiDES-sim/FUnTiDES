#ifndef FDTDUTILS_HPP
#define FDTDUTILS_HPP

#include <data_type.h>

#include "FDTDdata.hpp"
#include "FDTDkernels.hpp"

using namespace std;

using namespace std;

struct FDTDUtils
{
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
  float compute_dt_sch(const float vmax, vectorReal &coefx, vectorReal &coefy,
                       vectorReal &coefz)
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
    printf(" coefs computed\n");
    return 2 * cfl / (sqrtf(ftmp) * vmax);
  }

  void pml_profile_init(vector<float> &profile, int i_min, int i_max,
                        int n_first, int n_last, float scale)
  {
    int n = i_max - i_min + 1;
    int shift = i_min - 1;

    int first_beg = 1 + shift;
    int first_end = n_first + shift;
    int last_beg = n - n_last + 1 + shift;
    int last_end = n + shift;

#pragma omp parallel for
    for (int i = i_min; i <= i_max; ++i)
    {
      profile[i] = 0.f;
    }

    float tmp = scale / POW2(first_end - first_beg + 1);
#pragma omp parallel for
    for (int i = 1; i <= first_end - first_beg + 1; ++i)
    {
      profile[first_end - i + 1] = POW2(i) * tmp;
    }

#pragma omp parallel for
    for (int i = 1; i <= last_end - last_beg + 1; ++i)
    {
      profile[last_beg + i - 1] = POW2(i) * tmp;
    }
  }

  void pml_profile_extend(int nx, int ny, int nz, vectorReal &eta,
                          const vector<float> &etax, const vector<float> &etay,
                          const vector<float> &etaz, int xbeg, int xend,
                          int ybeg, int yend, int zbeg, int zend)
  {
    const int n_ghost = 1;
#pragma omp parallel for collapse(3)
    for (int ix = xbeg - n_ghost; ix <= xend + n_ghost; ++ix)
    {
      for (int iy = ybeg - n_ghost; iy <= yend + n_ghost; ++iy)
      {
        for (int iz = zbeg - n_ghost; iz <= zend + n_ghost; ++iz)
        {
          eta[(nz + 2) * (ny + 2) * ix + (nz + 2) * (iy) + iz] =
              etax[ix] + etay[iy] + etaz[iz];
        }
      }
    }
  }

  void pml_profile_extend_all(int nx, int ny, int nz, vectorReal &eta,
                              const vector<float> &etax,
                              const vector<float> &etay,
                              const vector<float> &etaz, int xmin, int xmax,
                              int ymin, int ymax, int x1, int x2, int x5,
                              int x6, int y1, int y2, int y3, int y4, int y5,
                              int y6, int z1, int z2, int z3, int z4, int z5,
                              int z6)
  {
    // Top.
    if (z1 != -1)
      pml_profile_extend(nx, ny, nz, eta, etax, etay, etaz, xmin, xmax, ymin,
                         ymax, z1, z2);
    // Bottom.
    if (z5 != -5)
      pml_profile_extend(nx, ny, nz, eta, etax, etay, etaz, xmin, xmax, ymin,
                         ymax, z5, z6);
    // Front.
    if ((y1 != -1) && (z3 != -3))
      pml_profile_extend(nx, ny, nz, eta, etax, etay, etaz, xmin, xmax, y1, y2,
                         z3, z4);
    // Back.
    if ((y6 != -6) && (z3 != -3))
      pml_profile_extend(nx, ny, nz, eta, etax, etay, etaz, xmin, xmax, y5, y6,
                         z3, z4);
    // Left.
    if ((x1 != -1) && (y3 != -3) && (z3 != -3))
      pml_profile_extend(nx, ny, nz, eta, etax, etay, etaz, x1, x2, y3, y4, z3,
                         z4);
    // Right.
    if ((x6 != -6) && (y3 != -3) && (z3 != -3))
      pml_profile_extend(nx, ny, nz, eta, etax, etay, etaz, x5, x6, y3, y4, z3,
                         z4);
  }

  void init_eta(int nx, int ny, int nz, int ndampx, int ndampy, int ndampz,
                int x1, int x2, int x3, int x4, int x5, int x6, int y1, int y2,
                int y3, int y4, int y5, int y6, int z1, int z2, int z3, int z4,
                int z5, int z6, float dx, float dy, float dz, float dt_sch,
                float vmax, vectorReal &eta)
  {
#pragma omp parallel for collapse(3)
    for (int i = -1; i < nx + 1; ++i)
    {
      for (int j = -1; j < ny + 1; ++j)
      {
        for (int k = -1; k < nz + 1; ++k)
        {
          eta[(nz + 2) * (ny + 2) * (i + 1) + (nz + 2) * (j + 1) + (k + 1)] =
              0.f;
        }
      }
    }

    vector<float> etax(nx + 2);
    vector<float> etay(ny + 2);
    vector<float> etaz(nz + 2);

    // etax
    float param = dt_sch * 3.f * vmax * logf(1000.f) / (2.f * ndampx * dx);
    // printf("param=%f\n",param);
    pml_profile_init(etax, 0, nx + 1, ndampx, ndampx, param);

    // etay
    param = dt_sch * 3.f * vmax * logf(1000.f) / (2.f * ndampy * dy);
    // printf("param=%f\n",param);
    pml_profile_init(etay, 0, ny + 1, ndampy, ndampy, param);

    // etaz
    param = dt_sch * 3.f * vmax * logf(1000.f) / (2.f * ndampz * dz);
    // printf("param=%f\n",param);
    pml_profile_init(etaz, 0, nz + 1, ndampz, ndampz, param);

    (void)pml_profile_extend_all(
        nx, ny, nz, eta, etax, etay, etaz, 1, nx, 1, ny, x1 + 1, x2, x5 + 1, x6,
        y1 + 1, y2, y3 + 1, y4, y5 + 1, y6, z1 + 1, z2, z3 + 1, z4, z5 + 1, z6);
  }

  void output(FDTDGRIDS &myGrids, arrayReal const &pnGlobal, int itSample,
              const int &i1, const bool saveSnapShots)
  {
    if (itSample % 50 == 0)
    {
      FDFENCE

      printf("TimeStep=%d\t; Pressure value at source [%d %d %d] = %f\n",
             itSample, myGrids.xs, myGrids.ys, myGrids.zs,
             pnGlobal(IDX3_l(myGrids.xs, myGrids.ys, myGrids.zs), i1));
      // #ifdef FD_SAVE_SNAPSHOTS
      if (saveSnapShots)
        write_io(myGrids, 0, myGrids.nx, myGrids.ny / 2, myGrids.ny / 2, 0,
                 myGrids.nz, pnGlobal, itSample, i1);
      // #endif
    }
  }

  void write_io(FDTDGRIDS &myGrids, int x0, int x1, int y0, int y1, int z0,
                int z1, arrayReal const &pnGlobal, int istep, const int &i1)
  {
    char filename_buf[32];
    snprintf(filename_buf, sizeof(filename_buf), "snapshot_it_%d.H@", istep);
    FILE *snapshot_file = fopen(filename_buf, "wb");
    // printf("write snapshot for: x=[%d %d], y=[%d %d], z=[%d
    // %d]\n",x0,x1,y0,y1,z0,z1);

    for (int k = z0; k < z1; ++k)
    {
      for (int j = y0; j < y1 + 1; ++j)
      {
        for (int i = x0; i < x1; ++i)
        {
          fwrite(&pnGlobal(IDX3_l(i, j, k), i1), sizeof(float), 1,
                 snapshot_file);
        }
      }
    }
    // Clean up
    fclose(snapshot_file);
  }
};

#endif  // FDTDUTILS_HPP
