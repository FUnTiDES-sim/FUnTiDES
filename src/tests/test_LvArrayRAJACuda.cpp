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

// test lvArray
using arrayRealCB=LvArray::Array< float,
                                  2,
                                  camp::idx_seq< 1, 0 >,
                                  std::ptrdiff_t,
                                  LvArray::ChaiBuffer >;
using arrayRealViewCB=LvArray::ArrayView< float,
                                          2,
                                          0,
                                          std::ptrdiff_t,
                                          LvArray::ChaiBuffer >;



// test ChaiBuffer ChaiBuffer
// -------------------------------
LVARRAY_HOST_DEVICE inline int  computeCBCB( int const n1, int const n2, arrayRealViewCB const & arrayViewIn, arrayRealViewCB const & arrayViewOut )
{
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      arrayViewOut[i][j]=sqrt( 2*i*j )*arrayViewIn[i][j];
    }
  }
  return(0);
}

///*
LVARRAY_HOST_DEVICE inline arrayRealCB   compute1CBCB( int const n1, int const n2, arrayRealViewCB const arrayViewIn )
{
  arrayRealCB array( n1, n2 );
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      array[i][j]=sqrt( 2*i*j )*arrayViewIn[i][j];
    }
  }
  return(array);
}
//*/

int main( int argc, char *argv[] )
{

  const int nIter=10000000;
  const int n1=1000;
  const int n2=1000;

  arrayRealCB arrayOut( n1, n2 );
  arrayRealCB arrayIn( n1, n2 );
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      arrayIn[i][j]=sqrt( 2*i*j );
    }
  }
  std::chrono::time_point< std::chrono::system_clock > startTimeCudaLoop = std::chrono::system_clock::now();
  arrayRealViewCB const arrayViewIn=arrayIn.toView();
  arrayRealViewCB const arrayViewOut=arrayOut.toView();
  std::cout<<"start loop1 \n";
  //using execPolicy=RAJA::omp_parallel_for_exec;
  using execPolicy=RAJA::cuda_exec< 32 >;
  RAJA::forall< execPolicy >( RAJA::RangeSegment( 0, nIter ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
    for( int i=0; i<n1; i++ )
    {
      for( int j=0; j<n2; j++ )
      {
        arrayViewOut[i][j]=sqrt( 2*i*j )*arrayViewIn[i][j];
      }
    }
  } );
  std::cout << "Elapsed Time RAJA Cuda loop : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeCudaLoop ).count() / 1000.0 <<" seconds.\n"<<std::endl;


  // test ChaiBuffer ChaiBuffer RAJA
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeRAJACBCB = std::chrono::system_clock::now();
  //using execPolicy=RAJA::omp_parallel_for_exec;
  using execPolicy1=RAJA::cuda_exec< 32 >;
  RAJA::forall< execPolicy1 >( RAJA::RangeSegment( 0, nIter ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
    int i=computeCBCB( n1, n2, arrayViewIn, arrayViewOut );
    //arrayRealViewCB arrayViewOut=compute1CBCB(n1,n2,arrayView).toView() ;
    //if (e==nIter-1)
    //  std::cout <<"Result RAJA CBCB "<<arrayViewOut[n1/2][n2/2]<<std::endl;
  } );
  std::cout << "Elapsed Time RAJACuda function call  loop CBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeRAJACBCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;
}
