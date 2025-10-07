//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include <cxxopts.hpp>
#include "fdtd_proxy.h"
#include "data_type.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>


fdtd_proxy::fdtd_proxy(const fdtd_options & opt)
    : m_opt(opt)
    , m_grids()
    , m_stencils()
    , m_kernels()
    , m_io()
    , m_utils()
    , m_solver(m_grids,m_kernels,m_stencils,m_source_receivers) {}

// Initialize the simulation.
// @post run()  
void fdtd_proxy::init_fdtd() 
{
  printf("+======================================\n");
  printf("saveSnapshots=%d snapShotInterval=%d\n",m_opt.output.saveSnapShots,
                                                  m_opt.output.snapShotInterval);
  printf("--------------------------------------\n");
  printf("\n");
  printf("geometry init\n");
  m_grids.InitGrid(m_opt);
  printf("--------------------------------------\n");
  printf("dx=%f dy=%f dz=%f\n",m_grids.dx(),m_grids.dy(),m_grids.dz());
  printf("nx=%d ny=%d nz=%d\n",m_grids.nx(),m_grids.ny(),m_grids.nz());

  printf("stencil init\n");
  printf("--------------------------------------\n");
  m_stencils.initStencilsCoefficients(m_opt,m_grids.dx(),m_grids.dy(),m_grids.dz());
  printf("stencil coefficients\n");
  printf("lx=%d ly=%d lz=%d\n",m_stencils.lx,m_stencils.ly,m_stencils.lz);  
  printf("coef0=%f\n",m_stencils.coef0);
  for( int i=0;i<m_stencils.ncoefsX;i++ )
    printf("coefx[%d]=%f ",i,m_stencils.coefx[i]);
  printf("\n");
  for( int i=0;i<m_stencils.ncoefsY;i++ )
    printf("coefy[%d]=%f ",i,m_stencils.coefy[i]);
  printf("\n");
  for( int i=0;i<m_stencils.ncoefsZ;i++ )
    printf("coefz[%d]=%f ",i,m_stencils.coefz[i]);
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
    timeStep=m_stencils.compute_dt_sch(vmax);
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
  m_grids.InitModelArrays(m_opt);
  printf("model init done\n");
  printf("--------------------------------------\n");

  //allocate and initialize wavefields arrays
  m_kernels.initFieldsArrays(m_grids.nx(), m_grids.ny(), m_grids.nz(),
                   m_stencils.lx, m_stencils.ly, m_stencils.lz);
  printf("arrays init done\n");
  printf("--------------------------------------\n"); 

  // set up source
  f0=m_opt.source.f0;
  sourceOrder=m_opt.source.sourceOrder;
  printf("central freq and source order\n");
  printf("f0=%f\n",f0);
  printf("sourceOrder=%d\n",sourceOrder);
  printf("--------------------------------------\n");
  m_source_receivers.xsrc=m_opt.source.xs;
  m_source_receivers.ysrc=m_opt.source.ys;
  m_source_receivers.zsrc=m_opt.source.zs;
  if(m_source_receivers.xsrc<0) m_source_receivers.xsrc=m_grids.nx()/2;
  if(m_source_receivers.ysrc<0) m_source_receivers.ysrc=m_grids.ny()/2;
  if(m_source_receivers.zsrc<0) m_source_receivers.zsrc=m_grids.nz()/2;
  printf("source position\n");
  printf("xsrc=%d ysrc=%d zsrc=%d\n",m_source_receivers.xsrc,m_source_receivers.ysrc,m_source_receivers.zsrc);
  printf("--------------------------------------\n");
  init_source();
  printf("source init done\n");
  printf("--------------------------------------\n");
  
  // define sponge boundary
  m_kernels.defineSpongeBoundary(m_grids.nx(), m_grids.ny(), m_grids.nz());
  printf("sponge boundary init done\n"); 
  printf("--------------------------------------\n");

  printf("solver initialization done\n");
  printf("--------------------------------------\n");
  }

// inti source term, e.g., Ricker wavelet
void fdtd_proxy::init_source()
{
  // compute source term
  m_kernels.RHSTerm = allocateVector< vectorReal >( nSamples, "RHSTerm" );

  std::vector< float > sourceTerm=m_utils.computeSourceTerm( nSamples, timeStep, f0, sourceOrder );
  for( int i=0; i<nSamples; i++ )
  {
    m_kernels.RHSTerm[i]=sourceTerm[i];
    //cout<<"sample "<<i<<"\t: sourceTerm = "<<sourceTerm[i]<< endl;
  }
}

// run all time steps
void fdtd_proxy::run() 
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime, totalOutputTime;
  for (int indexTimeSample = 0; indexTimeSample < nSamples; indexTimeSample++)
  {
    startComputeTime = system_clock::now();
    m_solver.compute_one_step(indexTimeSample,i1,i2);
    totalComputeTime += system_clock::now() - startComputeTime;
    startOutputTime = system_clock::now();
    if (indexTimeSample % m_opt.output.snapShotInterval == 0)
    {
      m_io.outputPnValues(indexTimeSample, i1,m_grids,m_kernels,m_stencils,m_opt,m_source_receivers);
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
