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
#include "Array.hpp"
#include "MallocBuffer.hpp"
#include "ChaiBuffer.hpp"

template< class T > class Array2D
{
public:
  Array2D( int numRows, int numCols ): data( numRows, std::vector< T >( numCols )) {}
  Array2D(): data( 0, std::vector< T >( 0 )) {}

  std::vector< T > & operator[]( int index ){return data[index];}
  T& operator()( size_t row, size_t col ) {return data[row][col];}

private:
  std::vector< std::vector< T > > data;
};

using arrayReal=Array2D< float >;
// test lvArray
using arrayRealMB=LvArray::Array< float,
                                  2,
                                  camp::idx_seq< 1, 0 >,
                                  std::ptrdiff_t,
                                  LvArray::MallocBuffer >;

using arrayRealCB=LvArray::Array< float,
                                  2,
                                  camp::idx_seq< 1, 0 >,
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

// test vector vector
//
// -------------------------------
void compute( int const n1, int const n2, arrayReal & arrayIn, arrayReal & array )
{
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      array( i, j )=sqrt( 2*i*j )*arrayIn[i][j];
    }
  }
}

// test MallocBuffer MallocBuffer
//
// -------------------------------
void computeMBMB( int const n1, int const n2, arrayRealMB const & arrayIn, arrayRealMB const & array )
{
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      array[i][j]=sqrt( 2*i*j )*arrayIn[i][j];
    }
  }
}

// test ChaiBuffer ChaiBuffer
// -------------------------------
void computeCBCB( int const n1, int const n2, arrayRealCB const & arrayIn, arrayRealCB const & array )
{
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      array[i][j]=sqrt( 2*i*j )*arrayIn[i][j];
    }
  }

}



int main( int argc, char *argv[] )
{


    << << <<< HEAD
    const int nIter=1000;
  =======
    const int nIter=10000;
  >> >> >>> feature/henri/allocateArray
  const int n1=1000;
  const int n2=1000;

  arrayReal arrayIn;
  arrayIn=allocateArray2D< arrayReal >( n1, n2 );
  arrayReal arrayOut;
  arrayOut=allocateArray2D< arrayReal >( n1, n2 );
  arrayRealCB arrayCBIn;
  arrayCBIn=allocateArray2D< arrayRealCB >( n1, n2 );
  arrayRealCB arrayCBOut;
  arrayCBOut=allocateArray2D< arrayRealCB >( n1, n2 );
  arrayRealMB arrayMBIn;
  arrayMBIn=allocateArray2D< arrayRealMB >( n1, n2 );
  arrayRealMB arrayMBOut;
  arrayMBOut=allocateArray2D< arrayRealMB >( n1, n2 );

  #pragma omp parallel for
  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      arrayIn[i][j]=sqrt( 2*i*j );
      arrayOut[i][j]=0;
    }
  }

  for( int i=0; i<n1; i++ )
  {
    for( int j=0; j<n2; j++ )
    {
      arrayCBIn[i][j]=sqrt( 2*i*j );
      arrayMBIn[i][j]=sqrt( 2*i*j );
      arrayCBOut[i][j]=0;
      arrayMBOut[i][j]=0;
    }
  }

  // test vector vector OMP
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTime = std::chrono::system_clock::now();
  #pragma omp parallel for
  for( int e=0; e<nIter; e++ )
  {
    compute( n1, n2, arrayIn, arrayOut );
    if( e==nIter-1 )
      std::cout <<"Result OMP vector "<<arrayOut[n1/2][n2/2]<<std::endl;
  }
  std::cout << "Elapsed Time OMP loop vector : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTime ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  // test vector vector OMP flat loop
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeflat  = std::chrono::system_clock::now();
  #pragma omp parallel for
  for( int e=0; e<nIter; e++ )
  {
    for( int i=0; i<n1; i++ )
    {
      for( int j=0; j<n2; j++ )
      {
        arrayOut( i, j )=sqrt( 2*i*j )*arrayIn[i][j];
      }
    }
    if( e==nIter-1 )
      std::cout <<"Result OMP vector "<<arrayOut[n1/2][n2/2]<<std::endl;
  }
  std::cout << "Elapsed Time OMP loop vector : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeflat ).count() / 1000.0 <<" seconds.\n"<<std::endl;


  // test MallocBuffer MallocBuffer OMP
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeMB = std::chrono::system_clock::now();
  #pragma omp parallel for
  for( int e=0; e<nIter; e++ )
  {
    computeMBMB( n1, n2, arrayMBIn, arrayMBOut );
    if( e==nIter-1 )
      std::cout <<"Result OMP MBMB "<<arrayMBOut[n1/2][n2/2]<<std::endl;
  }
  std::cout << "Elapsed Time OMP loop MBMB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeMB ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  // test ChaiBuffer ChaiBuffer OMP
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeCB = std::chrono::system_clock::now();
  #pragma omp parallel for
  for( int e=0; e<nIter; e++ )
  {
    computeCBCB( n1, n2, arrayCBIn, arrayCBOut );
    if( e==nIter-1 )
      std::cout <<"Result OMP CBCB "<<arrayCBOut[n1/2][n2/2]<<std::endl;
  }
  std::cout << "Elapsed Time OMP loop CBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;


  // test MallocBuffer MallocBuffer RAJA
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeRAJAMBMB = std::chrono::system_clock::now();
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, nIter ), [=] ( int e )
  {
    computeMBMB( n1, n2, arrayMBIn, arrayMBOut );
    if( e==nIter-1 )
      std::cout <<"Result RAJA MBMB "<<arrayMBOut[n1/2][n2/2]<<std::endl;
  } );
  std::cout << "Elapsed Time RAJA loop MBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeRAJAMBMB ).count() / 1000.0 <<" seconds.\n"<<std::endl;


  // test ChaiBuffer ChaiBuffer RAJA
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeRAJACBCB = std::chrono::system_clock::now();
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, nIter ), [=] ( int e )
  {
    computeCBCB( n1, n2, arrayCBIn, arrayCBOut );
    if( e==nIter-1 )
      std::cout <<"Result RAJA CBCB "<<arrayCBOut[n1/2][n2/2]<<std::endl;
  } );
  std::cout << "Elapsed Time RAJA loop CBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
    ( std::chrono::system_clock::now() - startTimeRAJACBCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;
}
