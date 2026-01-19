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
  void defineSpongeBoundary(int nx, int ny, int nz,int L=20,float alpha=-0.0015)
  {
    // Allocate spongeArray if not already allocated
    if (spongeArray.extent(0) == 0)
    {
      spongeArray = allocateVector<vectorReal>(nx * ny * nz, "spongeArray");
    }

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
    for (int k =0; k < nz ; k++)
    {
      for (int j = 0 ; j < ny; j++)
      {
        for (int i = 0; i < L; i++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-i)*dx,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - i), 2));
           //printf("spongeArray(%d,%d,%d)=%f\n",i,j,k,spongeArray(IDX3(i, j,k)));
        }
        for (int i = nx - L; i < nx; i++)
        {
          // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-(nx-i))*dx,2));
          spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - (nx - i)), 2));
        }
      }
    }

    // Y boundary
    for (int k = 0; k < nz ; k++)
    {
      for (int i = 0; i < nx; i++)
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
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        //for (int k = 0; k < L; k++)
        //{
        //  // spongeArray(IDX3(i,j,k))= exp(alpha*pow((L-k))*dz,2));
        //  spongeArray(IDX3(i, j, k)) = exp(alpha * pow((L - k), 2));
        //}
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
  void pml_profile_init(Kokkos::View<float*> profile, int i_min, int i_max,
                        int n_first, int n_last, float scale)
  {
    int n = i_max - i_min + 1;
    int shift = i_min - 1;

    int first_beg = 1 + shift;
    int first_end = n_first + shift;
    int last_beg = n - n_last + 1 + shift;
    int last_end = n + shift;

    // Zero initialization using Kokkos
    Kokkos::parallel_for("pml_profile_init_zero",
      Kokkos::RangePolicy<>(i_min, i_max + 1),
      KOKKOS_LAMBDA(int i) {
        profile(i) = 0.f;
      });

    float tmp = scale / POW2(first_end - first_beg + 1);

    // First boundary
    Kokkos::parallel_for("pml_profile_init_first",
      Kokkos::RangePolicy<>(1, first_end - first_beg + 2),
      KOKKOS_LAMBDA(int i) {
        profile(first_end - i + 1) = POW2(i) * tmp;
      });

    // Last boundary
    Kokkos::parallel_for("pml_profile_init_last",
      Kokkos::RangePolicy<>(1, last_end - last_beg + 2),
      KOKKOS_LAMBDA(int i) {
        profile(last_beg + i - 1) = POW2(i) * tmp;
      });

    Kokkos::fence();
  }

  void pml_profile_extend(int nx, int ny, int nz, vectorReal &eta,
                          const Kokkos::View<float*> &etax,
                          const Kokkos::View<float*> &etay,
                          const Kokkos::View<float*> &etaz,
                          int xbeg, int xend,
                          int ybeg, int yend, int zbeg, int zend)
  {
    const int n_ghost = 1;

    // Skip if any input bounds are negative (disabled regions)
    if (xbeg < 0 || xend < 0 || ybeg < 0 || yend < 0 || zbeg < 0 || zend < 0) {
      return;
    }

    // Ensure non-negative bounds for MDRangePolicy
    // Clamp to valid ranges [0, nx+2] for x, [0, ny+2] for y, [0, nz+2] for z
    int ix_start = (xbeg - n_ghost) > 0 ? (xbeg - n_ghost) : 0;
    int iy_start = (ybeg - n_ghost) > 0 ? (ybeg - n_ghost) : 0;
    int iz_start = (zbeg - n_ghost) > 0 ? (zbeg - n_ghost) : 0;

    int ix_end = (xend + n_ghost + 1) < (nx + 2) ? (xend + n_ghost + 1) : (nx + 2);
    int iy_end = (yend + n_ghost + 1) < (ny + 2) ? (yend + n_ghost + 1) : (ny + 2);
    int iz_end = (zend + n_ghost + 1) < (nz + 2) ? (zend + n_ghost + 1) : (nz + 2);

    // Skip if any dimension has invalid range
    if (ix_start >= ix_end || iy_start >= iy_end || iz_start >= iz_end) {
      return;
    }

    // Use Kokkos parallel_for with MDRangePolicy for 3D iteration
    Kokkos::parallel_for("pml_profile_extend",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        {ix_start, iy_start, iz_start},
        {ix_end, iy_end, iz_end}
      ),
      KOKKOS_LAMBDA(int ix, int iy, int iz) {
        eta((nz + 2) * (ny + 2) * ix + (nz + 2) * (iy) + iz) =
            etax(ix) + etay(iy) + etaz(iz);
      });

    Kokkos::fence();
  }

  void pml_profile_extend_all(int nx, int ny, int nz, vectorReal &eta,
                              const Kokkos::View<float*> &etax,
                              const Kokkos::View<float*> &etay,
                              const Kokkos::View<float*> &etaz,
                              int xmin, int xmax,
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
    // Zero initialization using Kokkos parallel_for
    // MDRangePolicy requires non-negative indices, so we iterate from 0
    Kokkos::parallel_for("init_eta_zero",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        {0, 0, 0},
        {nx + 2, ny + 2, nz + 2}
      ),
      KOKKOS_LAMBDA(int i, int j, int k) {
        eta((nz + 2) * (ny + 2) * i + (nz + 2) * j + k) = 0.f;
      });

    Kokkos::fence();

    // Allocate Kokkos::View instead of std::vector
    Kokkos::View<float*> etax("etax", nx + 2);
    Kokkos::View<float*> etay("etay", ny + 2);
    Kokkos::View<float*> etaz("etaz", nz + 2);

    // etax
    float param = dt_sch * 3.f * vmax * logf(1000.f) / (2.f * ndampx * dx);
    pml_profile_init(etax, 0, nx + 1, ndampx, ndampx, param);

    // etay
    param = dt_sch * 3.f * vmax * logf(1000.f) / (2.f * ndampy * dy);
    pml_profile_init(etay, 0, ny + 1, ndampy, ndampy, param);

    // etaz
    param = dt_sch * 3.f * vmax * logf(1000.f) / (2.f * ndampz * dz);
    pml_profile_init(etaz, 0, nz + 1, ndampz, ndampz, param);

    (void)pml_profile_extend_all(
        nx, ny, nz, eta, etax, etay, etaz, 1, nx, 1, ny, x1 + 1, x2, x5 + 1, x6,
        y1 + 1, y2, y3 + 1, y4, y5 + 1, y6, z1 + 1, z2, z3 + 1, z4, z5 + 1, z6);
  }
};

}  // namespace abckernel
}  // namespace fdtd
#endif  // FDTD_PML_HPP
