#ifndef FDTD_IO_H_
#define FDTD_IO_H_

#include <cstdio>
#include <vector>
#include "fdtd_grids.h"
#include "fdtd_kernels.h"
#include "fdtd_stencils.h"
#include "fdtd_options.h"
#include "fdtd_source_receivers.h"

using namespace std;

struct fdtd_io
{
    //writes  pn values at the source location and save snapshot
    void outputPnValues(int itSample, int i1,
                        fdtd_grids &m_grids,
                        fdtd_kernels &m_kernels,
                        fdtd_stencils &m_stencils,
                        fdtd_options &m_opt,
                        fdtd_source_receivers & m_src)  
    {
        int nx=m_grids.nx;
        int ny=m_grids.ny;
        int nz=m_grids.nz;
        int lx=m_stencils.lx;
        int ly=m_stencils.ly;
        int lz=m_stencils.lz;
        int xs=m_src.xsrc;
        int ys=m_src.ysrc;
        int zs=m_src.zsrc;
  
        bool saveSnapShots=m_opt.output.saveSnapShots;
        if(itSample%50==0)
        {    
            FDFENCE
            printf(
                "TimeStep=%d\t; Pressure value at source [%d %d %d] =%f %f %f\n",
                 itSample,
                 xs, ys, zs, 
                 m_kernels.pnGlobal(IDX3_l(xs-1,ys,zs),i1), 
                 m_kernels.pnGlobal(IDX3_l(xs,ys,zs),i1),
                 m_kernels.pnGlobal(IDX3_l(xs+1,ys,zs),i1 )
                ); 
            if(saveSnapShots) write_snapshot( 0, nx, ny/2, ny/2, 0, nz, itSample,i1,
                                              m_grids, m_kernels ,m_stencils,m_opt);
        } 
    }

    // write snapshot to file
    void write_snapshot(const int &x0, const int &x1, 
                       const int &y0, const int &y1, 
                       const int &z0, const int &z1, 
                       const int istep, int i1,
                       fdtd_grids &m_grids, fdtd_kernels &m_kernels,
                       fdtd_stencils &m_stencils,fdtd_options &m_opt)
    {
        int ny=m_grids.ny;
        int nz=m_grids.nz;
        int lx=m_stencils.lx;
        int ly=m_stencils.ly;
        int lz=m_stencils.lz;
        char filename_buf[32];
        snprintf(filename_buf, sizeof(filename_buf), "snapshot_it_%d.H@", istep);
    
        FILE* snapshot_file = fopen(filename_buf, "wb");
        if (!snapshot_file) {
            fprintf(stderr, "Error: Could not open file %s for writing\n", filename_buf);
            return;
        }

        // Buffer for batch writing
        std::vector<float> buffer;
        buffer.reserve((x1-x0) * (y1-y0+1) * (z1-z0));

        // Collect data into buffer
        for (int k = z0; k < z1; ++k) {
            for (int j = y0; j < y1+1; ++j) {
                for (int i = x0; i < x1; ++i) {
                    buffer.push_back(m_kernels.pnGlobal(IDX3_l(i,j,k), i1));
                }
            }
        }

        // Write entire buffer at once
        size_t elements_written = fwrite(buffer.data(), sizeof(float), buffer.size(), snapshot_file);
    
        if (elements_written != buffer.size()) {
            fprintf(stderr, "Error: Failed to write complete snapshot, wrote %zu of %zu elements\n",
                    elements_written, buffer.size());
        }

        fclose(snapshot_file);
    }
};
#endif //FDTD_IO.H_
