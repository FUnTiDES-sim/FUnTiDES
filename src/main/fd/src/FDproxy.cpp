//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include <cxxopts.hpp>
#include "FDproxy.hpp"
#include "data_type.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>


FDproxy::FDproxy(const FDProxyOptions& opt):m_opt(opt) {}

// Initialize the simulation.
// @post run()  
void FDproxy::initFDElem() 
{
  printf("+======================================\n");
  printf("saveSnapshots=%d snapShotInterval=%d\n",m_opt.output.saveSnapShots,m_opt.output.snapShotInterval);
  printf("--------------------------------------\n");
  printf("\n");
  printf("geometry init\n");
  myGrids.initGrid(m_opt);
  printf("--------------------------------------\n");
  printf("dx=%f dy=%f dz=%f\n",myGrids.dx,myGrids.dy,myGrids.dz);
  printf("nx=%d ny=%d nz=%d\n",myGrids.nx,myGrids.ny,myGrids.nz);

  printf("stencil init\n");
  printf("--------------------------------------\n");
  myStencils.initStencilsCoefficients(m_opt,myGrids.dx,myGrids.dy,myGrids.dz);
  printf("stencil coefficients\n");
  printf("lx=%d ly=%d lz=%d\n",myStencils.lx,myStencils.ly,myStencils.lz);  
  printf("coef0=%f\n",myStencils.coef0);
  for( int i=0;i<myStencils.ncoefsX;i++ )
    printf("coefx[%d]=%f ",i,myStencils.coefx[i]);
  printf("\n");
  for( int i=0;i<myStencils.ncoefsY;i++ )
    printf("coefy[%d]=%f ",i,myStencils.coefy[i]);
  printf("\n");
  for( int i=0;i<myStencils.ncoefsZ;i++ )
    printf("coefz[%d]=%f ",i,myStencils.coefz[i]);
  printf("\n");

  // velocity model
  printf("\n");
  printf("velocity model init\n");
  printf("vmin=%f vmax=%f\n",m_opt.velocity.vmin,m_opt.velocity.vmax);
  printf("--------------------------------------\n");
  vmin=m_opt.velocity.vmin;
  vmax=m_opt.velocity.vmax;
  lambdamax=m_opt.velocity.vmax/(2.5*m_opt.source.f0);
  timeStep=m_opt.time.timeStep;
  timeMax=m_opt.time.timeMax;
  printf("user defined time step=%e\n",timeStep); 
  printf("user defined max time=%f\n",m_opt.time.timeMax);
  printf("--------------------------------------\n");
  if(timeStep==0)
  {
    timeStep=myStencils.compute_dt_sch(vmax);
    printf("compute time step from CFL condition\n");
  }
  else
    printf("user defined time step\n"); 
  nSamples=timeMax/timeStep;
  printf("timeStep=%e\n",timeStep);
  printf("nSamples=%d\n",nSamples);
  printf("--------------------------------------\n");

  // init model arrays
  printf("model init\n");
  myGrids.initModelArrays(m_opt);
  printf("model init done\n");
  printf("--------------------------------------\n");

  //allocate and initialize wavefields arrays
  myKernels.initFieldsArrays(myGrids.nx, myGrids.ny, myGrids.nz,
                   myStencils.lx, myStencils.ly, myStencils.lz);
  printf("arrays init done\n");
  printf("--------------------------------------\n"); 

  // set up source
  f0=m_opt.source.f0;
  sourceOrder=m_opt.source.sourceOrder;
  printf("central freq and source order\n");
  printf("f0=%f\n",f0);
  printf("sourceOrder=%d\n",sourceOrder);
  printf("--------------------------------------\n");
  xsrc=m_opt.source.xs;
  ysrc=m_opt.source.ys;
  zsrc=m_opt.source.zs;
  if(xsrc<0) xsrc=myGrids.nx/2;
  if(ysrc<0) ysrc=myGrids.ny/2;
  if(zsrc<0) zsrc=myGrids.nz/2;
  printf("source position\n");
  printf("xsrc=%d ysrc=%d zsrc=%d\n",xsrc,ysrc,zsrc);
  printf("--------------------------------------\n");
  initSource();
  printf("source init done\n");
  printf("--------------------------------------\n");
  
  // define sponge boundary
  myKernels.defineSpongeBoundary(myGrids.nx, myGrids.ny, myGrids.nz);
  printf("sponge boundary init done\n"); 
  printf("--------------------------------------\n");
}


// inti source term, e.g., Ricker wavelet
void FDproxy::initSource()
{
  // compute source term
  myKernels.RHSTerm = allocateVector< vectorReal >( nSamples, "RHSTerm" );

  std::vector< float > sourceTerm=myUtils.computeSourceTerm( nSamples, timeStep, f0, sourceOrder );
  for( int i=0; i<nSamples; i++ )
  {
    myKernels.RHSTerm[i]=sourceTerm[i];
    //cout<<"sample "<<i<<"\t: sourceTerm = "<<sourceTerm[i]<< endl;
  }
}

// compute one time step
void FDproxy::computeOneStep(int itime)
{

  int x3=myStencils.lx;
  int x4=myGrids.nx - myStencils.lx;
  int y3=myStencils.ly;
  int y4=myGrids.ny - myStencils.ly;
  int z3=myStencils.lz;
  int z4=myGrids.nz - myStencils.lz;
  // add source term
  myKernels.addRHS( itime, i2,
                   myGrids.nx, myGrids.ny, myGrids.nz,
                   myStencils.lx, myStencils.ly, myStencils.lz,
                   xsrc, ysrc, zsrc,
                   myGrids.vp,
                   myKernels.RHSTerm,
                   myKernels.pnGlobal );
  //printf("addRHS done\n");
  FDFENCE
  // inner points
  myKernels.inner3D( i1, i2,
                    myGrids.nx, myGrids.ny, myGrids.nz,
                    myStencils.lx, myStencils.ly, myStencils.lz,
                    x3, x4, y3, y4, z3, z4,
                    myStencils.coef0,
                    myStencils.coefx,
                    myStencils.coefy,
                    myStencils.coefz,
                    myGrids.vp,
                    myKernels.pnGlobal );
  //printf("inner3D done\n");
  FDFENCE
  // apply sponge boundary to wavefield
  myKernels.applySponge( i1, i2,
                        myGrids.nx, myGrids.ny, myGrids.nz,
                        myStencils.lx, myStencils.ly, myStencils.lz,
                        x3, x4, y3, y4, z3, z4,
                        myKernels.spongeArray,
                        myKernels.pnGlobal );
  //printf("applySponge done\n"); 
  FDFENCE
}

// run all time steps
void FDproxy::run() 
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime, totalOutputTime;
  for (int indexTimeSample = 0; indexTimeSample < nSamples; indexTimeSample++)
  {
    startComputeTime = system_clock::now();
    computeOneStep(indexTimeSample);
    totalComputeTime += system_clock::now() - startComputeTime;
    startOutputTime = system_clock::now();
    if (indexTimeSample % m_opt.output.snapShotInterval == 0)
    {
      myIO.outputPnValues(indexTimeSample, i1,myGrids,myKernels,myStencils,m_opt);
    }
    //swap(i1, i2);
    auto tmp = i1;
    i1 = i2;
    i2 = tmp;
    totalOutputTime += system_clock::now() - startOutputTime;
   
    fflush(stdout);
  }
  float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                            .time_since_epoch()
                            .count();
  float outputtime_ms=time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();
  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."<< endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."<< endl;
  cout << "------------------------------------------------ " << endl;
}

//************************************************************************
// End of file
//************************************************************************