#ifndef FDTD_PML_HPP
#define FDTD_PML_HPP

#include "data_type.h"
#include "fd_macros.h"

using namespace std;

namespace fdtd
{
namespace abckernel
{
struct FdtdAbcKernels
{
  vectorReal eta;
  vectorReal spongeArray;
  //------------------------------------------------------------------
  // sponge boundary
  //------------------------------------------------------------------
  // define sponge boundary
  void defineSpongeBoundary(int nx, int ny, int nz)
  {
    const int L = 20;
    const float alpha = -0.00015;
    // const float alpha = -0.00035;
    //  compute sponge boundary terms
    //  intailize to 1
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
          // printf("spongeArray(%d,%d,%d)=%f\n",i,j,k,spongeArray(IDX3(i, j,
          // k)));
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

  //------------------------------------------------------------------
  // PML
  //------------------------------------------------------------------
  // Initialize the PML profile
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
};

}  // namespace abckernel
}  // namespace fdtd
#endif  // FDTD_PML_HPP
