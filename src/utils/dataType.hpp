#ifndef DATATYPE_HPP_
#define DATATYPE_HPP_

#include <fstream>
#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>

#include "commonMacros.hpp"

#ifdef USE_VECTOR

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

  T & operator[]( int index ) const
  {
    return const_cast< T & >(data[index]);
  }
  T & operator=( const T & data )
  { return *this; };

  int size(){return this->size();};

private:
  std::vector< T > data;
};

template< class T > class Array2D
{
public:
  Array2D( int numRows, int numCols ): data( numRows, std::vector< T >( numCols, 0 )) {}
  Array2D(): data( 0, std::vector< T >( 0 )) {}

  std::vector< T > & operator[]( int index ){return data[index];}
  T& operator()( int row, int col ) {return data[row][col];}
  T& operator()( int row, int col ) const {return const_cast< T & >(data[row][col]);}
  T& operator=( const T & data ) { return *this; };

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

using vectorReal=Vector< float >;
using vectorInt=Vector< int >;
using vectorDouble=Vector< double >;

using arrayInt=Array2D< int >;
using arrayReal=Array2D< float >;
using arrayDouble=Array2D< double >;

using array3DInt=Array3D< int >;
using array3DReal=Array3D< float >;
using array3DDouble=Array3D< double >;

#endif //USE_VECTOR

#ifdef USE_LVARRAY

#include "RAJA/RAJA.hpp"
#include "Array.hpp"
#include "ChaiBuffer.hpp"
using hostExecPolicy=RAJA::omp_parallel_for_exec;
using hostAtomicPolicy=RAJA::omp_atomic;
#ifdef ENABLE_CUDA
using deviceAtomicPolicy=RAJA::cuda_atomic;
using deviceExecPolicy=RAJA::cuda_exec< 32 >;
#endif
#ifdef ENABLE_HIP
using deviceAtomicPolicy=RAJA::hip_atomic;
using deviceExecPolicy=RAJA::hip_exec< 32 >;
#endif
// define vectors and  arrays.
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
// defines Views
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
                                       1,
                                       std::ptrdiff_t,
                                       LvArray::ChaiBuffer >;
using arrayRealView=LvArray::ArrayView< float,
                                        2,
                                        1,
                                        std::ptrdiff_t,
                                        LvArray::ChaiBuffer >;
using arrayDoubleView=LvArray::ArrayView< double,
                                          2,
                                          1,
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

#endif //USE_LVARRAY

#ifdef USE_KOKKOS

#ifdef ENABLE_HIP
#define __HIP_PLATFORM_AMD__ 1
#endif

#include <Kokkos_Core.hpp>
#define MemSpace Kokkos::SharedSpace

#ifdef ENABLE_CUDA
using Layout=Kokkos::LayoutLeft;
#else
using Layout=Kokkos::LayoutRight;
#endif

typedef Kokkos::View< int *, Layout, MemSpace > vectorInt;
typedef Kokkos::View< float *, Layout, MemSpace > vectorReal;
typedef Kokkos::View< double *, Layout, MemSpace > vectorDouble;

typedef Kokkos::View< int * *, Layout, MemSpace > arrayInt;
typedef Kokkos::View< float * *, Layout, MemSpace > arrayReal;
typedef Kokkos::View< double * *, Layout, MemSpace > arrayDouble;

typedef Kokkos::View< int * * *, Layout, MemSpace > array3DInt;
typedef Kokkos::View< float * * *, Layout, MemSpace > array3DReal;
typedef Kokkos::View< double * * *, Layout, MemSpace > array3DDouble;

/*
   typedef Kokkos::View<int*,     Layout, DeviceMemorySpace> vectorIntViewView;
   typedef Kokkos::View<float*,   Layout,  DeviceMemorySpace> vectorRealViewView;
   typedef Kokkos::View<double*,  Layout, DeviceMemorySpace> vectorDoubleView;

   typedef Kokkos::View<int**,    Layout, DeviceMemorySpace> arrayIntView;
   typedef Kokkos::View<float**,  Layout, DeviceMemorySpace> arrayRealView;
   typedef Kokkos::View<double**, Layout, DeviceMemorySpace> arrayDoubleView;

   typedef Kokkos::View<int***,    Layout, DeviceMemorySpace> array3DIntView;
   typedef Kokkos::View<float***,  Layout, DeviceMemorySpace> array3DRealView;
   typedef Kokkos::View<double***, Layout, DeviceMemorySpace> array3DDoubleView;
 */

#endif //USE_KOKKOS

template< class T >
T allocateVector( int n1 )
{
  #ifdef PRINT_ALLOC_INFO
  std::cout<<"allocate vector of size "<<n1<<std::endl;
  #endif
  T vect( KOKKOSNAME n1 );
  return vect;
}
template< class T >
T allocateVector( int n1, const char * name )
{
  #ifdef PRINT_ALLOC_INFO
  std::cout<<"allocate vector: "<<name<<" of size: "<<n1<<std::endl;
  #endif
  T vect( KOKKOSNAME n1 );
  return vect;
}
template< class T >
T allocateArray2D( int n1, int n2 )
{
  #ifdef PRINT_ALLOC_INFO
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
  #endif
  T array( KOKKOSNAME n1, n2 );
  return array;
}
template< class T >
T allocateArray2D( int n1, int n2, const char * name )
{
  #ifdef PRINT_ALLOC_INFO
  std::cout<<"allocate array : "<<name<<" of size: ("<<n1<<", "<<n2<<")"<<std::endl;
  #endif
  T array( KOKKOSNAME n1, n2 );
  return array;
}
template< class T >
T allocateArray3D( int n1, int n2, int n3 )
{
  #ifdef PRINT_ALLOC_INFO
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<", "<<n3<<std::endl;
  #endif
  T array( KOKKOSNAME n1, n2, n3 );
  return array;
}


#endif //DATATYPE_HPP_
