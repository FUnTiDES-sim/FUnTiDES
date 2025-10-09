#ifndef FDTDINIT_HPP
#define FDTDINIT_HPP

#include <utils.h>

#include "FDTDutils.hpp"

struct FDTDInit
{
  FDTDUtils myFDTDUtils;
  SolverUtils myUtils;

  int sourceOrder = 1;
  int ncoefsX;
  int ncoefsY;
  int ncoefsZ;

  float f0 = 15.;
  float fmax = 2.5 * f0;
  float timeMax = 1.0;
  float timeStep;
  int nSamples;

  float vmin = 1500;
  float vmax = 4500;

  int i1 = 0;
  int i2 = 1;

  bool saveSnapShots = false;
  bool usePML = false;

  vectorReal RHSTerm;

  void init_geometry(int argc, char *argv[], FDTDGRIDS &myGrids)
  {
    myGrids.nx = 150;
    myGrids.ny = myGrids.nx;
    myGrids.nz = myGrids.nx;

    myGrids.xs = myGrids.nx / 2;
    myGrids.ys = myGrids.ny / 2;
    myGrids.zs = myGrids.nz / 2;
    myGrids.lx = 4;
    myGrids.ly = 4;
    myGrids.lz = 4;
    ncoefsX = 5;
    ncoefsY = 5;
    ncoefsZ = 5;

    myGrids.dx = 10;
    myGrids.dy = 10;
    myGrids.dz = 10;
    if (argc > 1)
    {
      int i = 1;
      while (i < argc)
      {
        // dimensions nx,ny,nz
        std::string arg = argv[i];
        if (arg == "-nx")
        {
          myGrids.nx = atoi(argv[i + 1]);
          myGrids.ny = myGrids.nx;
          myGrids.nz = myGrids.nx;
          i = i + 2;
        }
        if (arg == "-ny")
        {
          myGrids.ny = atoi(argv[i + 1]);
          i = i + 2;
        }
        if (arg == "-nz")
        {
          myGrids.nz = atoi(argv[i + 1]);
          i = i + 2;
        }

        // spatial sampling dx dy dz
        if (arg == "-dx")
        {
          myGrids.dx = atoi(argv[i + 1]);
          myGrids.dy = myGrids.dx;
          myGrids.dz = myGrids.dx;
          i = i + 2;
        }
        if (arg == "-dy")
        {
          myGrids.dy = atoi(argv[i + 1]);
          i = i + 2;
        }
        if (arg == "-dz")
        {
          myGrids.dz = atoi(argv[i + 1]);
          i = i + 2;
        }

        // half stencil length lx ly lz
        if (arg == "-lx")
        {
          myGrids.lx = atoi(argv[i + 1]);
          myGrids.ly = myGrids.lx;
          myGrids.lz = myGrids.lx;
          i = i + 2;
        }
        if (arg == "-ly")
        {
          myGrids.ly = atoi(argv[i + 1]);
          i = i + 2;
        }
        if (arg == "-lz")
        {
          myGrids.lz = atoi(argv[i + 1]);
          i = i + 2;
        }
        ncoefsX = myGrids.lx + 1;
        ncoefsY = myGrids.ly + 1;
        ncoefsZ = myGrids.lz + 1;

        // source location xs ys zs
        if (arg == "-xs")
        {
          myGrids.xs = atoi(argv[i + 1]);
          myGrids.ys = myGrids.xs;
          myGrids.zs = myGrids.xs;
          i = i + 2;
        }
        if (arg == "-ys")
        {
          myGrids.ys = atoi(argv[i + 1]);
          i = i + 2;
        }
        if (arg == "-zs")
        {
          myGrids.zs = atoi(argv[i + 1]);
          i = i + 2;
        }

        // save snapshots
        if (arg == "-savesnapshots")
        {
          saveSnapShots = true;
          i = i + 1;
        }

        //  use PML
        if (arg == "-usePML")
        {
          usePML = true;
          i = i + 1;
        }
      }
    }
    else
    {
      printf(
          "usage fd_exe -nx XXX -ny YYY -nz ZZZZ -lx LX -ly LY -lz LZ -dx DX "
          "-dy DY -dz DZ  ");
      printf("-xs XS -ys YS -zs ZS -savesnapshots (false) -usePML (false) \n ");
      printf("if no argument running with default values\n");
    }

    printf("Number of grids: nx=%d, ny=%d, nz=%d\n", myGrids.nx, myGrids.ny,
           myGrids.nz);
    printf("Source location: xs=%d, ys=%d, zs=%d\n", myGrids.xs, myGrids.ys,
           myGrids.zs);
    float lambdamax = vmin / fmax;

    // init pml limits
    myGrids.ntaperx = 4;
    myGrids.ntapery = 4;
    myGrids.ntaperz = 4;

    myGrids.hdx_2 = 1. / (4. * myGrids.dx * myGrids.dx);
    myGrids.hdy_2 = 1. / (4. * myGrids.dy * myGrids.dy);
    myGrids.hdz_2 = 1. / (4. * myGrids.dz * myGrids.dz);

    printf("usePML=%d\n", usePML);
    if (usePML)
    {
      myGrids.ndampx = myGrids.ntaperx * lambdamax / myGrids.dx;
      myGrids.ndampy = myGrids.ntapery * lambdamax / myGrids.dy;
      myGrids.ndampz = myGrids.ntaperz * lambdamax / myGrids.dz;
    }
    else
    {
      myGrids.ndampx = 0;
      myGrids.ndampy = 0;
      myGrids.ndampz = 0;
    }
    myGrids.x1 = 0;
    myGrids.x2 = myGrids.ndampx;
    myGrids.x3 = myGrids.ndampx;
    myGrids.x4 = myGrids.nx - myGrids.ndampx;
    myGrids.x5 = myGrids.nx - myGrids.ndampx;
    myGrids.x6 = myGrids.nx;

    myGrids.y1 = 0;
    myGrids.y2 = myGrids.ndampy;
    myGrids.y3 = myGrids.ndampy;
    myGrids.y4 = myGrids.ny - myGrids.ndampy;
    myGrids.y5 = myGrids.ny - myGrids.ndampy;
    myGrids.y6 = myGrids.ny;

    myGrids.z1 = 0;
    myGrids.z2 = myGrids.ndampz;
    myGrids.z3 = myGrids.ndampz;
    myGrids.z4 = myGrids.nz - myGrids.ndampz;
    myGrids.z5 = myGrids.nz - myGrids.ndampz;
    myGrids.z6 = myGrids.nz;
  }

  void init_coefficients(FDTDGRIDS &myGrids, FDTDMODELS &myModels)
  {
    myModels.coefx = allocateVector<vectorReal>(ncoefsX, "coefx");
    myModels.coefy = allocateVector<vectorReal>(ncoefsY, "coefy");
    myModels.coefz = allocateVector<vectorReal>(ncoefsZ, "coefz");

    myFDTDUtils.init_coef(myGrids.lx, myGrids.dx, myModels.coefx);
    myFDTDUtils.init_coef(myGrids.ly, myGrids.dy, myModels.coefy);
    myFDTDUtils.init_coef(myGrids.lz, myGrids.dz, myModels.coefz);

    float tmpX = 0;
    for (int i = 1; i < ncoefsX; i++)
    {
      tmpX += myModels.coefx[i];
    }
    float tmpY = 0;
    for (int i = 1; i < ncoefsY; i++)
    {
      tmpY += myModels.coefy[i];
    }
    float tmpZ = 0;
    for (int i = 1; i < ncoefsZ; i++)
    {
      tmpZ += myModels.coefz[i];
    }
    myModels.coef0 = -2. * (tmpX + tmpY + tmpZ);

    timeStep = myFDTDUtils.compute_dt_sch(vmax, myModels.coefx, myModels.coefy,
                                          myModels.coefz);
    nSamples = timeMax / timeStep;
    printf("init coefs done\n");
  }

  void init_source(FDTDMODELS &myModels)
  {
    // compute source term
    myModels.RHSTerm = allocateVector<vectorReal>(nSamples, "RHSTerm");

    std::vector<float> sourceTerm =
        myUtils.computeSourceTerm(nSamples, timeStep, f0, sourceOrder);
    for (int i = 0; i < nSamples; i++)
    {
      myModels.RHSTerm[i] = sourceTerm[i];
      // cout<<"sample "<<i<<"\t: sourceTerm = "<<sourceTerm[i]<< endl;
    }
  }

  void init_models(FDTDGRIDS &myGrids, FDTDMODELS &myModels)
  {
    int modelVolume = myGrids.nx * myGrids.ny * myGrids.nz;
    int extModelVolume = (myGrids.nx + 2 * myGrids.lx) *
                         (myGrids.ny + 2 * myGrids.ly) *
                         (myGrids.nz + 2 * myGrids.lz);
    int etaModelVolume = (myGrids.nx + 2) * (myGrids.ny + 2) * (myGrids.nz + 2);

    myModels.spongeArray =
        allocateVector<vectorReal>(modelVolume, "spongeArray");
    myModels.vp = allocateVector<vectorReal>(modelVolume, "vp");
    myModels.pnp1 = allocateVector<vectorReal>(extModelVolume, "pnp1");
    myModels.pn = allocateVector<vectorReal>(extModelVolume, "pn");
    myModels.pnGlobal =
        allocateArray2D<arrayReal>(extModelVolume, 2, "pnGlobal");

    if (usePML)
    {
      printf("ici\n");
      myModels.phi = allocateVector<vectorReal>(modelVolume, "phi");
      myModels.eta = allocateVector<vectorReal>(etaModelVolume, "eta");
      myFDTDUtils.init_eta(
          myGrids.nx, myGrids.ny, myGrids.nz, myGrids.ndampx, myGrids.ndampy,
          myGrids.ndampz, myGrids.x1, myGrids.x2, myGrids.x3, myGrids.x4,
          myGrids.x5, myGrids.x6, myGrids.y1, myGrids.y2, myGrids.y3,
          myGrids.y4, myGrids.y5, myGrids.y6, myGrids.z1, myGrids.z2,
          myGrids.z3, myGrids.z4, myGrids.z5, myGrids.z6, myGrids.dx,
          myGrids.dy, myGrids.dz, timeStep, vmax, myModels.eta);
      printf("ici\n");
    }
    float init_vp_value = vmin * vmin * timeStep * timeStep;

    // initialize vp and pressure field
    for (int i = 0; i < myGrids.nx; i++)
    {
      for (int j = 0; j < myGrids.ny; j++)
      {
        for (int k = 0; k < myGrids.nz; k++)
        {
          myModels.vp[IDX3(i, j, k)] = init_vp_value;
        }
      }
    }
    for (int i = 0; i < myGrids.nx; i++)
    {
      for (int j = 0; j < myGrids.ny; j++)
      {
        for (int k = myGrids.nz / 2; k < myGrids.nz; k++)
        {
          myModels.vp[IDX3(i, j, k)] = 2 * init_vp_value;
        }
      }
    }

    for (int i = -myGrids.lx; i < myGrids.nx + myGrids.lx; i++)
    {
      for (int j = -myGrids.ly; j < myGrids.ny + myGrids.ly; j++)
      {
        for (int k = -myGrids.lz; k < myGrids.nz + myGrids.lz; k++)
        {
          myModels.pnp1[IDX3_l(i, j, k)] = 0.000001;
          myModels.pn[IDX3_l(i, j, k)] = 0.000001;
          myModels.pnGlobal(IDX3_l(i, j, k), 0) = 0.000001;
          myModels.pnGlobal(IDX3_l(i, j, k), 1) = 0.000001;
        }
      }
    }
  }

  void defineSpongeBoundary(FDTDGRIDS &myGrids, vectorReal &spongeArray)
  {
    const int L = 20;
    const int nx = myGrids.nx;
    const int ny = myGrids.ny;
    const int nz = myGrids.nz;
    // const float dx=myGrids.dx;
    // const float dy=myGrids.dy;
    // const float dz=myGrids.dz;
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
};

#endif  // FDTDINIT_HPP
