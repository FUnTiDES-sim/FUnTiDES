//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <chrono>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <omp.h>
#include <vector>

#define MemSpace Kokkos::SharedSpace
using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy< ExecSpace >;


// define data containers
typedef Kokkos::LayoutLeft layout;
typedef Kokkos::View< int *, layout, MemSpace >   vectorInt;
typedef Kokkos::View< float *, layout, MemSpace >   vectorReal;
typedef Kokkos::View< double *, layout, MemSpace >   vectorDouble;
typedef Kokkos::View< int * *, layout, MemSpace >   arrayInt;
typedef Kokkos::View< float * *, layout, MemSpace >   arrayReal;
typedef Kokkos::View< double * * *, layout, MemSpace >   array3DReal;
typedef Kokkos::View< double * *, layout, MemSpace >   arrayDouble;
typedef Kokkos::View< double * * *, layout, MemSpace >   array3DDouble;

template< class T >
T allocateVector( int n1 )
{
  std::cout<<"allocate vector of size "<<n1<<std::endl;
  T vect( "v", n1 );
  return vect;
}
template< class T >
T allocateArray2D( int n1, int n2 )
{
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
  T array( "a", n1, n2 );
  return array;
}
template< class T >
T allocateArray3D( int n1, int n2, int n3 )
{
  std::cout<<"allocate array of size "<<n1<<" ,"<<n2<<" ,"<<n3<<std::endl;
  T array( "a", n1, n2, n3 );
  return array;
}

//define source term

float evaluateRicker( float const & time_n, float const & f0, int order )
{
  float const o_tpeak = 1.0/f0;
  float pulse = 0.0;
  if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
  {
    return pulse;
  }

  constexpr float pi = M_PI;
  float const lam = (f0*pi)*(f0*pi);

  switch( order )
  {
    case 2:
    {
      pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)
              *exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
    }
    break;
    case 1:
    {
      pulse = -2.0*lam*(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
    }
    break;
    case 0:
    {
      pulse = -(time_n-o_tpeak)*exp( -2*lam*(time_n-o_tpeak)*(time_n-o_tpeak));
    }
    break;
    default:
      std::cout<<"This option is not supported yet, rickerOrder must be 0, 1 or 2"<<std::endl;
      break;
  }

  return pulse;
}

std::vector< float > computeSourceTerm( const int nSamples, const float timeStep, const float f0, const int order )
{
  std::vector< float > sourceTerm( nSamples );
  for( int i=0; i<nSamples; i++ )
  {
    float time_n=i*timeStep;
    sourceTerm[i]=evaluateRicker( time_n, f0, order );
  }
  return sourceTerm;
}

int main( int argc, char *argv[] )
{
  Kokkos::initialize( argc, argv );
  {
    const int n1=208;
    const int n2=208;
    const int n3=208;
    const float dx=10;

    int sourceOrder=1;
    int xs=n1/2;
    int ys=n2/2;
    int zs=n3/2;
    float f0=10.;
    float timeMax=2.0;

    const int ncoefs=5;
    vectorReal coef;
    coef=allocateVector< vectorReal >( ncoefs );
    float dx2 = dx*dx;
    coef[0] =-205./72.;
    coef[1] =8./5;
    coef[2] =-1./5.;
    coef[3] =8./315.;
    coef[4] =-1./560.;
    float vmax=1500;
    float cfl=0.80;

    double coef0=-6.*(coef[1]+coef[2]+coef[3]+coef[4]);
    float ftmp = 0;
    ftmp += fabsf( coef[0] ) + fabsf( coef[0] ) + fabsf( coef[0] );
    for( int i = 1; i < ncoefs; i++ )
    {
      ftmp += 2.f*fabsf( coef[i] );
      ftmp += 2.f*fabsf( coef[i] );
      ftmp += 2.f*fabsf( coef[i] );
    }
    float timeStep=2.*cfl*dx/(sqrtf( ftmp )*vmax);
    float timeStep2=timeStep*timeStep;
    const int nSamples=timeMax/timeStep;
    printf( "timeStep=%f\n", timeStep );

    vectorReal RHSTerm;
    RHSTerm=allocateVector< vectorReal >( nSamples );
    // compute source term
    std::vector< float >sourceTerm=computeSourceTerm( nSamples, timeStep, f0, sourceOrder );
    for( int i=0; i<nSamples; i++ )
    {
      RHSTerm[i]=sourceTerm[i];
    }

    // allocate vector and arrays
    array3DReal vp;
    vp=allocateArray3D< array3DReal >( n1, n1, n3 );
    array3DReal pnp1;
    pnp1=allocateArray3D< array3DReal >( n1, n2, n3 );
    array3DReal pn;
    pn=allocateArray3D< array3DReal >( n1, n2, n3 );
    array3DReal pnm1;
    pnm1=allocateArray3D< array3DReal >( n1, n2, n3 );

    printf( "memory used for vectra and arrays %d bytes\n", (4*n1*n2*n3+nSamples+ncoefs)*4 );

    Kokkos::Timer timer1;
    // initialize vp and pressure field
    Kokkos::parallel_for( Kokkos::MDRangePolicy< Kokkos::Rank< 3 > >( {0, 0, 0}, {n1, n2, n3} ), KOKKOS_LAMBDA ( int i, int j, int k )
    {
      vp( i, j, k )=1500.*1500./dx2;
      pnp1( i, j, k )=0.;
      pn( i, j, k )=0.;
      pnm1( i, j, k )=0.;
    } );

    for( int itSample=0; itSample<nSamples; itSample++ )
    {
      Kokkos::parallel_for( Kokkos::MDRangePolicy< Kokkos::Rank< 3 > >( {xs, xs, zs}, {xs+1, ys+1, zs+1} ), KOKKOS_LAMBDA ( int i, int j, int k )
      {
        pn( i, j, k )+=vp( i, j, k )*timeStep*timeStep*RHSTerm[itSample];
      } );
      Kokkos::parallel_for( Kokkos::MDRangePolicy< Kokkos::Rank< 3 > >( {4, 4, 4}, {n1-4, n2-4, n3-4} ), KOKKOS_LAMBDA ( int i, int j, int k )
      {
        float lapx=coef[1]*(pn( i+1, j, k )+pn( i-1, j, k ))
                    +coef[2]*(pn( i+2, j, k )+pn( i-2, j, k ))
                    +coef[3]*(pn( i+3, j, k )+pn( i-3, j, k ))
                    +coef[4]*(pn( i+4, j, k )+pn( i-4, j, k ));
        float lapy=coef[1]*(pn( i, j+1, k )+pn( i, j-1, k ))
                    +coef[2]*(pn( i, j+2, k )+pn( i, j-2, k ))
                    +coef[3]*(pn( i, j+3, k )+pn( i, j-3, k ))
                    +coef[4]*(pn( i, j+4, k )+pn( i, j-4, k ));
        float lapz=coef[1]*(pn( i, j, k+1 )+pn( i, j, k-1 ))
                    +coef[2]*(pn( i, j, k+2 )+pn( i, j, k-2 ))
                    +coef[3]*(pn( i, j, k+3 )+pn( i, j, k-3 ))
                    +coef[4]*(pn( i, j, k+4 )+pn( i, j, k-4 ));
        pnp1( i, j, k )=2.*pn( i, j, k )-pnm1( i, j, k )+timeStep2*vp( i, j, k )*(coef0*pn( i, j, k )+lapx+lapy+lapz);
        //if(i==xs && j==ys && k==zs)printf("%f %f %f\n",coef0*pn(i,j,k),lapx+lapy+lapz,pn(i,j,k));
        //if(i==xs && j==ys && k==zs)printf("exact %f\n",coef0+6*(coef[1]+coef[2]+coef[3]+coef[4]));
      } );
      Kokkos::fence();
      if( itSample%50==0 )
      {
        printf( "result 1 %f\n", pnp1( xs, ys, zs ));
      }
      Kokkos::parallel_for( Kokkos::MDRangePolicy< Kokkos::Rank< 3 > >( {0, 0, 0}, {n1, n2, n3} ), KOKKOS_LAMBDA ( int i, int j, int k )
      {
        pnm1( i, j, k )=pn( i, j, k );
        pn( i, j, k )=pnp1( i, j, k );
      } );
    }
    Kokkos::fence();
    double time1=timer1.seconds();
    std::cout << "Elapsed Time parallel_for MDRANGE kokkos  loop : "<<time1 <<" seconds.\n\n";

  }
  Kokkos::finalize();
  return (0);
}
