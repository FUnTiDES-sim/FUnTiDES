#ifndef FDTDDATA_HPP_
#define FDTDDATA_HPP_

#include "FDTDmacros.hpp"

struct FDTDGRIDS
{
  int nx, ny, nz;  // number of grid on the x,y,z direction
  int xs, ys, zs;  // source location on the x,y,z direction
  int lx, ly, lz;
  float dx, dy, dz;

  int ntaperx, ntapery, ntaperz;
  int ndampx, ndampy, ndampz;
  float hdx_2, hdy_2, hdz_2;

  int x1, x2, x3, x4, x5, x6;
  int y1, y2, y3, y4, y5, y6;
  int z1, z2, z3, z4, z5, z6;
};

struct FDTDMODELS
{
  double coef0;

  vectorReal coefx;
  vectorReal coefy;
  vectorReal coefz;
  vectorReal RHSTerm;
  vectorReal spongeArray;

  vectorReal vp;
  vectorReal phi;
  vectorReal eta;
  vectorReal pnp1;
  vectorReal pn;
  arrayReal pnGlobal;
};

#endif  // FDTDDATA_HPP_
