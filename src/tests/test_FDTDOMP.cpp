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
#include <omp.h>
#include <vector>

template< class T > class Vector
{
public:
  Vector( int numRows ): data( numRows ) {}
  Vector(): data( 0 ) {}

  T & operator()( int index )
  {
    return data[index];
  }
  T & operator[]( int index )
  {
    return data[index];
  }
  int size(){return this->size();}

private:
  std::vector< T > data;
};

template< class T > class Array2D
{
public:
  Array2D( int numRows, int numCols ): data( numRows, std::vector< T >( numCols, 0 )) {}
  Array2D(): data( 0, std::vector< T >( 0 )) {}

  std::vector< T > & operator[]( int index ){return data[index];}
  T& operator()( size_t row, size_t col ) {return data[row][col];}

private:
  std::vector< std::vector< T > > data;
};
template< class T > class Array3D
{
public:
  Array3D( int X, int Y, int Z ): data( X, std::vector< std::vector< T > >( Y, std::vector< T >( Z ))) {}
  Array3D(): data( 0, std::vector< std::vector< T > >( 0 )) {}

  std::vector< T > & operator[]( int index ){return data[index];}
  T& operator()( size_t X, size_t Y, size_t Z ) {return data[X][Y][Z];}

private:
  std::vector< std::vector< std::vector< T > > > data;
};

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
using vectorInt=std::vector< int >;
using vectorReal=std::vector< float >;
using vectorDouble=std::vector< double >;
using arrayInt=Array2D< int >;
using arrayReal=Array2D< float >;
using arrayDouble=Array2D< double >;
using array3DInt=Array3D< int >;
using array3DReal=Array3D< float >;
using array3DDouble=Array3D< double >;

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

  float coef0=-6.*(coef[1]+coef[2]+coef[3]+coef[4]);
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

  // initialize vp and pressure field
    #pragma omp parallel for collapse(3)
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      for( int k=0; k<n3; k++ )
      {
        vp( i, j, k )=1500.*1500./dx2;
        pnp1( i, j, k ) = 0;
        pn( i, j, k ) = 0;
        pnm1( i, j, k ) = 0;
      }
    }
  }

  for( int itSample=0; itSample<nSamples; itSample++ )
  {
    pn( xs, ys, zs )+=vp( xs, ys, zs )*timeStep*timeStep*RHSTerm[itSample];

      #pragma omp parallel for collapse(3)
    for( int i=4; i<n1-4; i++ )
    {
      for( int j=4; j<n2-4; j++ )
      {
        for( int k=4; k<n3-4; k++ )
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
        }
      }
    }
    if( itSample%50==0 )
    {
      printf( "result 1 %f\n", pnp1( xs, ys, zs ));
    }
      #pragma omp parallel for collapse(3)
    for( int i=0; i<n1; i++ )
    {
      for( int j=0; j<n2; j++ )
      {
        for( int k=0; k<n3; k++ )
        {
          pnm1( i, j, k )=pn( i, j, k );
          pn( i, j, k )=pnp1( i, j, k );
        }
      }
    }

  }
  return (0);
}
