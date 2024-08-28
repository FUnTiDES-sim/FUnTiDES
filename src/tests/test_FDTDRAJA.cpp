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
#include <vector>
#include "Array.hpp"
#include "RAJA/RAJA.hpp"
#include "Macros.hpp"
#include "ChaiBuffer.hpp"

// Create an 1D array of integers.
using vectorInt=LvArray::Array< int,
                                1,
                                camp::idx_seq< 0 >,
                                std::ptrdiff_t,
                                LvArray::ChaiBuffer >;
using vectorReal=LvArray::Array< float,
                                 1,
                                 camp::idx_seq< 0 >,
                                 std::ptrdiff_t,
                                 LvArray::ChaiBuffer >;
using vectorDouble=LvArray::Array< double,
                                   1,
                                   camp::idx_seq< 0 >,
                                   std::ptrdiff_t,
                                   LvArray::ChaiBuffer >;
using arrayInt=LvArray::Array< int,
                               2,
                               camp::idx_seq< 0, 1 >,
                               std::ptrdiff_t,
                               LvArray::ChaiBuffer >;
using arrayReal=LvArray::Array< float,
                                2,
                                camp::idx_seq< 0, 1 >,
                                std::ptrdiff_t,
                                LvArray::ChaiBuffer >;
using arrayDouble=LvArray::Array< double,
                                  2,
                                  camp::idx_seq< 0, 1 >,
                                  std::ptrdiff_t,
                                  LvArray::ChaiBuffer >;
using array3DInt=LvArray::Array< int,
                                 3,
                                 camp::idx_seq< 0, 1, 2 >,
                                 std::ptrdiff_t,
                                 LvArray::ChaiBuffer >;
using array3DReal=LvArray::Array< float,
                                  3,
                                  camp::idx_seq< 0, 1, 2 >,
                                  std::ptrdiff_t,
                                  LvArray::ChaiBuffer >;
using array3DDouble=LvArray::Array< double,
                                    3,
                                    camp::idx_seq< 0, 1, 2 >,
                                    std::ptrdiff_t,
                                    LvArray::ChaiBuffer >;

using vectorIntView=LvArray::ArrayView< int,
                                        1,
                                        0,
                                        std::ptrdiff_t,
                                        LvArray::ChaiBuffer >;
using vectorRealView=LvArray::ArrayView< float,
                                         1,
                                         0,
                                         std::ptrdiff_t,
                                         LvArray::ChaiBuffer >;
using arrayIntView=LvArray::ArrayView< int,
                                       2,
                                       0,
                                       std::ptrdiff_t,
                                       LvArray::ChaiBuffer >;
using arrayRealView=LvArray::ArrayView< float,
                                        2,
                                        0,
                                        std::ptrdiff_t,
                                        LvArray::ChaiBuffer >;
using arrayDoubleView=LvArray::ArrayView< double,
                                          2,
                                          0,
                                          std::ptrdiff_t,
                                          LvArray::ChaiBuffer >;
using array3DIntView=LvArray::ArrayView< int,
                                         3,
                                         2,
                                         std::ptrdiff_t,
                                         LvArray::ChaiBuffer >;
using array3DRealView=LvArray::ArrayView< float,
                                          3,
                                          2,
                                          std::ptrdiff_t,
                                          LvArray::ChaiBuffer >;
using array3DDoubleView=LvArray::ArrayView< double,
                                            3,
                                            2,
                                            std::ptrdiff_t,
                                            LvArray::ChaiBuffer >;
template< class T >
T allocateVector( int n1 )
{
  std::cout<<"allocate vector of size "<<n1<<std::endl;
  T vect( n1 );
  return vect;
}
template< class T >
T allocateArray2D( int n1, int n2 )
{
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
  T array( n1, n2 );
  return array;
}
template< class T >
T allocateArray3D( int n1, int n2, int n3 )
{
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<", "<<n3<<std::endl;
  T array( n1, n2, n3 );
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
  const int n1=208;
  const int n2=208;
  const int n3=208;
  const float dx=10;

  const int sourceOrder=1;
  const int xs=n1/2;
  const int ys=n2/2;
  const int zs=n3/2;
  const float f0=10.;
  const float timeMax=2.0;

  const int ncoefs=5;
  vectorReal h_coef;
  h_coef=allocateVector< vectorReal >( ncoefs );
  float dx2 = dx*dx;
  h_coef[0] =-205./72.;
  h_coef[1] =8./5;
  h_coef[2] =-1./5.;
  h_coef[3] =8./315.;
  h_coef[4] =-1./560.;
  float vmax=1500;
  float cfl=0.80;

  double coef0=-6.*(h_coef[1]+h_coef[2]+h_coef[3]+h_coef[4]);
  float ftmp = 0;
  ftmp += fabsf( h_coef[0] ) + fabsf( h_coef[0] ) + fabsf( h_coef[0] );
  for( int i = 1; i < ncoefs; i++ )
  {
    ftmp += 2.f*fabsf( h_coef[i] );
    ftmp += 2.f*fabsf( h_coef[i] );
    ftmp += 2.f*fabsf( h_coef[i] );
  }
  float timeStep=2.*cfl*dx/(sqrtf( ftmp )*vmax);
  float timeStep2=timeStep*timeStep;
  const int nSamples=timeMax/timeStep;
  printf( "timeStep=%f\n", timeStep );

  vectorReal h_RHSTerm;
  h_RHSTerm=allocateVector< vectorReal >( nSamples );
  // compute source term
  std::vector< float >sourceTerm=computeSourceTerm( nSamples, timeStep, f0, sourceOrder );
  for( int i=0; i<nSamples; i++ )
  {
    h_RHSTerm[i]=sourceTerm[i];
  }

  // allocate vector and arrays
  array3DReal h_vp;
  h_vp=allocateArray3D< array3DReal >( n1, n1, n3 );
  array3DReal h_pnp1;
  h_pnp1=allocateArray3D< array3DReal >( n1, n2, n3 );
  array3DReal h_pn;
  h_pn=allocateArray3D< array3DReal >( n1, n2, n3 );
  array3DReal h_pnm1;
  h_pnm1=allocateArray3D< array3DReal >( n1, n2, n3 );

  printf( "memory used for vectra and arrays %d bytes\n", (4*n1*n2*n3+nSamples+ncoefs)*4 );

#pragma omp parallel for collapse(3)
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      for( int k=0; k<n3; k++ )
      {
        h_vp( i, j, k )=1500.*1500./dx2;
        h_pnp1( i, j, k )=0.;
        h_pn( i, j, k )=0.;
        h_pnm1( i, j, k )=0.;
      }
    }
  }

  vectorRealView const coef=h_coef.toView();
  vectorRealView const RHSTerm=h_RHSTerm.toView();
  array3DRealView const vp=h_vp.toView();
  array3DRealView const pnp1=h_pnp1.toView();
  array3DRealView const pn=h_pn.toView();
  array3DRealView const pnm1=h_pnm1.toView();

  //RAJA_INDEX_VALUE_T(KIDX, int, "KIDX");
  //RAJA_INDEX_VALUE_T(JIDX, int, "JIDX");
  //RAJA_INDEX_VALUE_T(IIDX, int, "IIDX");
  constexpr int imins = xs;
  constexpr int imaxs = xs+1;
  constexpr int jmins = ys;
  constexpr int jmaxs = ys+1;
  constexpr int kmins = zs;
  constexpr int kmaxs = zs+1;
  RAJA::TypedRangeSegment< int > KRanges( kmins, kmaxs );
  RAJA::TypedRangeSegment< int > JRanges( jmins, jmaxs );
  RAJA::TypedRangeSegment< int > IRanges( imins, imaxs );
  constexpr int imini = 4;
  constexpr int imaxi = n1-4;
  constexpr int jmini = 4;
  constexpr int jmaxi = n2-4;
  constexpr int kmini = 4;
  constexpr int kmaxi = n3-4;
  RAJA::TypedRangeSegment< int > KRangei( kmini, kmaxi );
  RAJA::TypedRangeSegment< int > JRangei( jmini, jmaxi );
  RAJA::TypedRangeSegment< int > IRangei( imini, imaxi );
  constexpr int imin = 0;
  constexpr int imax = n1;
  constexpr int jmin = 0;
  constexpr int jmax = n2;
  constexpr int kmin = 0;
  constexpr int kmax = n3;
  RAJA::TypedRangeSegment< int > KRange( kmin, kmax );
  RAJA::TypedRangeSegment< int > JRange( jmin, jmax );
  RAJA::TypedRangeSegment< int > IRange( imin, imax );

  using EXEC_POL5 =
    RAJA::KernelPolicy<
      RAJA::statement::CudaKernel<
        RAJA::statement::For< 2, RAJA::cuda_thread_x_loop,      // k
                              RAJA::statement::For< 1, RAJA::cuda_thread_y_loop, // j
                                                    RAJA::statement::For< 0, RAJA::cuda_thread_z_loop, // i
                                                                          RAJA::statement::Lambda< 0 >
                                                                          >
                                                    >
                              >
        >
      >;

  for( int itSample=0; itSample<nSamples; itSample++ )
  {
    RAJA::kernel< EXEC_POL5 >( RAJA::make_tuple( IRanges, JRanges, KRanges ), [=] __device__ (int i, int j, int k)
    {
      pn( i, j, k )+=vp( i, j, k )*timeStep*timeStep*RHSTerm[itSample];
    } );

    RAJA::kernel< EXEC_POL5 >( RAJA::make_tuple( IRangei, JRangei, KRangei ), [=] __device__ (int i, int j, int k)
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
    if( itSample%50==0 )
    {
      RAJA::forall< RAJA::loop_exec >( RAJA::RangeSegment( 0, n1 ), [pnp1] ( int i)
      {} );
      printf( "result 1 %f\n", pnp1( xs, ys, zs ));
    }

    RAJA::kernel< EXEC_POL5 >( RAJA::make_tuple( IRange, JRange, KRange ), [=] __device__ (int i, int j, int k)
    {
      pnm1( i, j, k )=pn( i, j, k );
      pn( i, j, k )=pnp1( i, j, k );
    } );
  }

  return (0);
}
