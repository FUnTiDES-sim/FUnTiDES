//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include "RAJA/RAJA.hpp"
#include "Macros.hpp"
#include "Array.hpp"
#include "MallocBuffer.hpp"
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
using vectorDoubleView=LvArray::ArrayView< double,
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



int main( int argc, char *argv[] )
{

  const int n1=500;
  const int n2=1000;
  const int n3=1000;
  const int ncoefs=4;
  float dt=0.001;

  vectorReal h_coef( ncoefs );
//  h_coef=allocateVector<vectorReal>(ncoefs);
  array3DReal h_pnp1( n1, n2, n3 );
//  h_pnp1=allocateArray3D<array3DReal>(n1,n2,n3);
  array3DReal h_pn( n1, n2, n3 );
//  h_pn=allocateArray3D<array3DReal>(n1,n2,n3);
  array3DReal h_pnm1( n1, n2, n3 );
//  h_pnm1=allocateArray3D<array3DReal>(n1,n2,n3);

 #pragma omp parallel for
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      for( int k=0; k<n3; k++ )
      {
        h_pnp1( i, j, k )=1;
        h_pn( i, j, k )=1;
        h_pnm1( i, j, k )=1;
      }
    }
  }

  h_coef( 0 )=4*dt*dt;
  h_coef( 1 )=1*dt*dt;
  h_coef( 2 )=1*dt*dt;
  h_coef( 3 )=1*dt*dt;

  std::chrono::time_point< std::chrono::system_clock > startTimeSequentialLoop = std::chrono::system_clock::now();
  for( int i=3; i<n1-3; i++ )
  {
    for( int j=3; j<n2-3; j++ )
    {
      for( int k=3; k<n3-3; k++ )
      {
        h_pnp1( i, j, k )=(2. +8*h_coef( 0 ))*h_pn( i, j, k )+h_coef( 1 )*(h_pn( i+1, j, k )+h_pn( i-1, j, k )+h_pn( i, j+1, k )+h_pn( i, j-1, k )+h_pn( i, j, k+1 )+h_pn( i, j, k-1 ))
                           +h_coef( 2 )*(h_pn( i+2, j, k )+h_pn( i-2, j, k )+h_pn( i, j+2, k )+h_pn( i, j-2, k )+h_pn( i, j, k+2 )+h_pn( i, j, k-2 ))
                           +h_coef( 3 )*(h_pn( i+3, j, k )+h_pn( i-3, j, k )+h_pn( i, j+3, k )+h_pn( i, j-3, k )+h_pn( i, j, k+3 )+h_pn( i, j, k-3 ))
                           -h_pnm1( i, j, k );
      }
    }
  }
  printf( "result 1 %f\n", h_pnp1( n1/2, n2/2, n3/2 ));
  std::cout << "Elapsed Time Sequential loop : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeSequentialLoop ).count() / 1000.0 <<" seconds.\n"<<std::endl;

 #pragma omp parallel for
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      for( int k=0; k<n3; k++ )
      {
        h_pnp1( i, j, k )=1;
        h_pn( i, j, k )=1;
        h_pnm1( i, j, k )=1;
      }
    }
  }

  std::chrono::time_point< std::chrono::system_clock > startTimeOMPLoop = std::chrono::system_clock::now();
 #pragma omp parallel for
  for( int i=3; i<n1-3; i++ )
  {
    for( int j=3; j<n2-3; j++ )
    {
      for( int k=3; k<n3-3; k++ )
      {
        h_pnp1( i, j, k )=(2. +8*h_coef( 0 ))*h_pn( i, j, k )+h_coef( 1 )*(h_pn( i+1, j, k )+h_pn( i-1, j, k )+h_pn( i, j+1, k )+h_pn( i, j-1, k )+h_pn( i, j, k+1 )+h_pn( i, j, k-1 ))
                           +h_coef( 2 )*(h_pn( i+2, j, k )+h_pn( i-2, j, k )+h_pn( i, j+2, k )+h_pn( i, j-2, k )+h_pn( i, j, k+2 )+h_pn( i, j, k-2 ))
                           +h_coef( 3 )*(h_pn( i+3, j, k )+h_pn( i-3, j, k )+h_pn( i, j+3, k )+h_pn( i, j-3, k )+h_pn( i, j, k+3 )+h_pn( i, j, k-3 ))
                           -h_pnm1( i, j, k );
      }
    }
  }
  printf( "result 2 %f\n", h_pnp1( n1/2, n2/2, n3/2 ));
  std::cout << "Elapsed Time OMP loop : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeOMPLoop ).count() / 1000.0 <<" seconds.\n"<<std::endl;

 #pragma omp parallel for
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      for( int k=0; k<n3; k++ )
      {
        h_pnp1( i, j, k )=1;
        h_pn( i, j, k )=1;
        h_pnm1( i, j, k )=1;
      }
    }
  }
  std::chrono::time_point< std::chrono::system_clock > startTimeOMPCollapseLoop = std::chrono::system_clock::now();
 #pragma omp parallel for collapse(3)
  for( int i=3; i<n1-3; i++ )
  {
    for( int j=3; j<n2-3; j++ )
    {
      for( int k=3; k<n3-3; k++ )
      {
        h_pnp1( i, j, k )=(2. +8*h_coef( 0 ))*h_pn( i, j, k )+h_coef( 1 )*(h_pn( i+1, j, k )+h_pn( i-1, j, k )+h_pn( i, j+1, k )+h_pn( i, j-1, k )+h_pn( i, j, k+1 )+h_pn( i, j, k-1 ))
                           +h_coef( 2 )*(h_pn( i+2, j, k )+h_pn( i-2, j, k )+h_pn( i, j+2, k )+h_pn( i, j-2, k )+h_pn( i, j, k+2 )+h_pn( i, j, k-2 ))
                           +h_coef( 3 )*(h_pn( i+3, j, k )+h_pn( i-3, j, k )+h_pn( i, j+3, k )+h_pn( i, j-3, k )+h_pn( i, j, k+3 )+h_pn( i, j, k-3 ))
                           -h_pnm1( i, j, k );
      }
    }
  }
  printf( "result 3 %f\n", h_pnp1( n1/2, n2/2, n3/2 ));
  std::cout << "Elapsed Time OMP collapse loop : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeOMPCollapseLoop ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  h_pnp1( n1/2, n2/2, n3/2 )=1;
  h_pn( n1/2, n2/2, n3/2 )=1;
  h_pnm1( n1/2, n2/2, n3/2 )=1;
  vectorRealView const coef= h_coef.toView();
  array3DRealView const pnp1= h_pnp1.toView();
  array3DRealView const pn=h_pn.toView();
  array3DRealView const pnm1=h_pnm1.toView();


  using execPolicy=RAJA::cuda_exec< 32 >;
  RAJA::forall< execPolicy >( RAJA::RangeSegment( 0, n1 ), [=] LVARRAY_HOST_DEVICE ( int i )
  {
    for( int j=0; j<n2; j++ )
    {
      for( int k=0; k<n3; k++ )
      {
        pnp1( i, j, k )=1;
        pn( i, j, k )=1;
        pnm1( i, j, k )=1;
      }
    }
  } );

  std::chrono::time_point< std::chrono::system_clock > startTimeRAJALoop = std::chrono::system_clock::now();
  RAJA::forall< execPolicy >( RAJA::RangeSegment( 3, n1-3 ), [=] LVARRAY_HOST_DEVICE ( int i )
  {
    for( int j=3; j<n2-3; j++ )
    {
      for( int k=3; k<n3-3; k++ )
      {
        pnp1( i, j, k )=(2. +8*coef( 0 ))*pn( i, j, k )+coef( 1 )*(pn( i+1, j, k )+pn( i-1, j, k )+pn( i, j+1, k )+pn( i, j-1, k )+pn( i, j, k+1 )+pn( i, j, k-1 ))
                         +coef( 2 )*(pn( i+2, j, k )+pn( i-2, j, k )+pn( i, j+2, k )+pn( i, j-2, k )+pn( i, j, k+2 )+pn( i, j, k-2 ))
                         +coef( 3 )*(pn( i+3, j, k )+pn( i-3, j, k )+pn( i, j+3, k )+pn( i, j-3, k )+pn( i, j, k+3 )+pn( i, j, k-3 ))
                         -pnm1( i, j, k );
      }
    }
  } );
  printf( "result 4 %f\n", h_pnp1( n1/2, n2/2, n3/2 ));
  std::cout << "Elapsed Time RAJA loop : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeRAJALoop ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  //RAJA_INDEX_VALUE_T(KIDX, int, "KIDX");
  //RAJA_INDEX_VALUE_T(JIDX, int, "JIDX");
  //RAJA_INDEX_VALUE_T(IIDX, int, "IIDX");
  constexpr int imin = 3;
  constexpr int imax = n1-3;
  constexpr int jmin = 3;
  constexpr int jmax = n2-3;
  constexpr int kmin = 3;
  constexpr int kmax = n3-3;
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
  RAJA::kernel< EXEC_POL5 >( RAJA::make_tuple( IRange, JRange, KRange ), [=] __device__ (int i, int j, int k) {
    //printf( " (%d, %d, %d) \n", (int)(*i), (int)(*j), (int)(*k));
    pnp1( i, j, k )=1;
    pn( i, j, k )=1;
    pnm1( i, j, k )=1;
  } );

  std::chrono::time_point< std::chrono::system_clock > startTimeRAJANestedLoop = std::chrono::system_clock::now();
  RAJA::kernel< EXEC_POL5 >( RAJA::make_tuple( IRange, JRange, KRange ), [=] __device__ (int i, int j, int k) {
    pnp1( i, j, k )=(2. +8*coef( 0 ))*pn( i, j, k )+coef( 1 )*(pn( i+1, j, k )+pn( i-1, j, k )+pn( i, j+1, k )+pn( i, j-1, k )+pn( i, j, k+1 )+pn( i, j, k-1 ))
                     +coef( 2 )*(pn( i+2, j, k )+pn( i-2, j, k )+pn( i, j+2, k )+pn( i, j-2, k )+pn( i, j, k+2 )+pn( i, j, k-2 ))
                     +coef( 3 )*(pn( i+3, j, k )+pn( i-3, j, k )+pn( i, j+3, k )+pn( i, j-3, k )+pn( i, j, k+3 )+pn( i, j, k-3 ))
                     -pnm1( i, j, k );
  } );
  RAJA::forall< RAJA::loop_exec >( RAJA::RangeSegment( 0, n1 ), [pnp1] ( int i)
  {} );
  printf( "result 5 %f\n", pnp1( n1/2, n2/2, n3/2 ));
  std::cout << "Elapsed Time RAJA Nested loop : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeRAJANestedLoop ).count() / 1000.0 <<" seconds.\n"<<std::endl;
}
