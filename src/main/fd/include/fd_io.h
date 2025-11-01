#ifndef FDTD_IO_H_
#define FDTD_IO_H_

#include <cstdio>
#include <vector>

#include "fd_grids.h"
#include "fd_kernels.h"
#include "fd_options.h"
#include "fd_source_receivers.h"
#include "fd_stencils.h"

using namespace std;
namespace fdtd
{
namespace io
{
struct FdtdIo
{
  // writes  pn values at the source location and save snapshot
  void outputPnValues(int itSample, int i1, model::fdgrid::FdtdGrids &m_grids,
                      fdtd::kernel::FdtdKernels &m_kernels,
                      fdtd::stencils::FdtdStencils &m_stencils,
                      fdtd::options::FdtdOptions &m_opt,
                      FdtdSourceReceivers &m_src)
  {
    int nx = m_grids.nx();
    int ny = m_grids.ny();
    int nz = m_grids.nz();
    int lx = m_stencils.lx;
    int ly = m_stencils.ly;
    int lz = m_stencils.lz;
    int xs = m_src.xsrc;
    int ys = m_src.ysrc;
    int zs = m_src.zsrc;

    bool saveSnapShots = m_opt.output.save_snapshots;
    if (itSample % 50 == 0)
    {
      FDFENCE
      printf("TimeStep=%d\t; Pressure value at source [%d %d %d] =%f %f %f\n",
             itSample, xs, ys, zs,
             m_kernels.pnGlobal(IDX3_l(xs - 1, ys, zs), i1),
             m_kernels.pnGlobal(IDX3_l(xs, ys, zs), i1),
             m_kernels.pnGlobal(IDX3_l(xs + 1, ys, zs), i1));
      if (saveSnapShots)
        WriteSnapshot(0, nx, ny / 2, ny / 2, 0, nz, itSample, i1, m_grids,
                      m_kernels, m_stencils, m_opt);
    }
  }

  // write snapshot to file
  void WriteSnapshot(const int &x0, const int &x1, const int &y0, const int &y1,
                     const int &z0, const int &z1, const int istep, int i1,
                     model::fdgrid::FdtdGrids &m_grids,
                     fdtd::kernel::FdtdKernels &m_kernels,
                     fdtd::stencils::FdtdStencils &m_stencils,
                     fdtd::options::FdtdOptions &m_opt)
  {
    int ny = m_grids.ny();
    int nz = m_grids.nz();
    int lx = m_stencils.lx;
    int ly = m_stencils.ly;
    int lz = m_stencils.lz;
    char filename_buf[32];
    snprintf(filename_buf, sizeof(filename_buf), "snapshot_it_%d.H@", istep);

    FILE *snapshot_file = fopen(filename_buf, "wb");
    if (!snapshot_file)
    {
      fprintf(stderr, "Error: Could not open file %s for writing\n",
              filename_buf);
      return;
    }

    // Buffer for batch writing
    std::vector<float> buffer;
    buffer.reserve((x1 - x0) * (y1 - y0 + 1) * (z1 - z0));

    // Collect data into buffer
    for (int k = z0; k < z1; ++k)
    {
      for (int j = y0; j < y1 + 1; ++j)
      {
        for (int i = x0; i < x1; ++i)
        {
          buffer.push_back(m_kernels.pnGlobal(IDX3_l(i, j, k), i1));
        }
      }
    }

    // Write entire buffer at once
    size_t elements_written =
        fwrite(buffer.data(), sizeof(float), buffer.size(), snapshot_file);

    if (elements_written != buffer.size())
    {
      fprintf(stderr,
              "Error: Failed to write complete snapshot, wrote %zu of %zu "
              "elements\n",
              elements_written, buffer.size());
    }

    fclose(snapshot_file);
  }
};

}  // namespace io
}  // namespace fdtd

#endif  // FDTD_IO.H_
