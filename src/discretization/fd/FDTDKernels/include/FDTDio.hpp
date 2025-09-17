#ifndef FDTDIO_HPP
#define FDTDIO_HPP

#include <cstdio>
#include <vector>
#include "FDTDGrids.hpp"
#include "FDTDkernels.hpp"
#include "FDTDStencils.hpp"
#include "FDproxyOptions.hpp"

using namespace std;

struct FDTDio
{
    //writes  pn values at the source location and save snapshot
    void outputPnValues(int itSample, int i1,
                        FDTDGRIDS &myGrids,
                        FDTDKernels &myKernels,
                        FDTDStencils &myStencils,
                        FDProxyOptions &m_opt)  
    {
        int nx=myGrids.nx;
        int ny=myGrids.ny;
        int nz=myGrids.nz;
        int lx=myStencils.lx;
        int ly=myStencils.ly;
        int lz=myStencils.lz;
        int xs=100;
        int ys=100;
        int zs=100;
  
        bool saveSnapShots=m_opt.output.saveSnapShots;
        if(itSample%50==0)
        {    
            FDFENCE
            printf(
                "TimeStep=%d\t; Pressure value at source [%d %d %d] =%f %f %f\n",
                 itSample,
                 xs, ys, zs, 
                 myKernels.pnGlobal(IDX3_l(xs-1,ys,zs),i1), 
                 myKernels.pnGlobal(IDX3_l(xs,ys,zs),i1),
                 myKernels.pnGlobal(IDX3_l(xs+1,ys,zs),i1 )
                ); 
            if(saveSnapShots) write_snapshot( 0, nx, ny/2, ny/2, 0, nz, itSample,i1,
                                              myGrids, myKernels ,myStencils,m_opt);
        } 
    }

    // write snapshot to file
    void write_snapshot(const int &x0, const int &x1, 
                       const int &y0, const int &y1, 
                       const int &z0, const int &z1, 
                       const int istep, int i1,
                       FDTDGRIDS &myGrids, FDTDKernels &myKernels,
                       FDTDStencils &myStencils,FDProxyOptions &m_opt)
    {
        int ny=myGrids.ny;
        int nz=myGrids.nz;
        int lx=myStencils.lx;
        int ly=myStencils.ly;
        int lz=myStencils.lz;
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
                    buffer.push_back(myKernels.pnGlobal(IDX3_l(i,j,k), i1));
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
#endif //FDTDio.hpp