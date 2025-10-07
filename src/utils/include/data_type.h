#ifndef DATATYPE_HPP_
#define DATATYPE_HPP_

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "common_macros.h"

using real_t = float;

using namespace std;

#ifdef USE_DOUBLE
using real_t = double;
#else
using real_t = float;
#endif

#ifdef USE_VECTOR

template <class T>
class Vector
{
 public:
  Vector(int numRows) : data_(numRows) {}
  Vector() : data_(0) {}

  // Iterator range constructor for copying from other containers
  template <typename Iterator>
  Vector(Iterator first, Iterator last) : data_(first, last)
  {
  }

  // Element access
  T &operator()(int index) { return data_[index]; }
  T &operator[](int index) { return data_[index]; }
  const T &operator[](int index) const { return data_[index]; }

  // Assignment operator (fixed - was incorrect)
  Vector &operator=(const Vector &other)
  {
    if (this != &other)
    {
      data_ = other.data_;
    }
    return *this;
  }

  // Size/extent
  size_t extent(int dim) const { return data_.size(); }

  // Data access
  T *data() { return data_.data(); }
  const T *data() const { return data_.data(); }

 private:
  std::vector<T> data_;
};

template <class T>
class Array2D
{
 public:
  Array2D(int numRows, int numCols) : data(numRows, std::vector<T>(numCols, 0))
  {
  }
  Array2D() : data(0, std::vector<T>(0)) {}

  std::vector<T> &operator[](int index) { return data[index]; }
  T &operator()(int row, int col) { return data[row][col]; }
  T &operator()(int row, int col) const
  {
    return const_cast<T &>(data[row][col]);
  }
  T &operator=(const T &data) { return *this; };

  size_t extent(int dim) const
  {
    if (dim == 0) return data.size();
    if (dim == 1 && !data.empty()) return data[0].size();
    return 0;
  }

  std::vector<T> getColumn(int colIndex) const
  {
    if (data.empty() || colIndex >= data[0].size())
    {
      return {};  // Empty vector
    }

    std::vector<T> column(data.size());
    for (size_t i = 0; i < data.size(); ++i)
    {
      column[i] = data[i][colIndex];
    }
    return column;
  }

 private:
  std::vector<std::vector<T>> data;
};

template <class T>
class Array3D
{
 public:
  Array3D(int X, int Y, int Z)
      : data(X, std::vector<std::vector<T>>(Y, std::vector<T>(Z)))
  {
  }
  Array3D() : data(0, std::vector<std::vector<T>>(0)) {}

  std::vector<T> &operator[](int index) { return data[index]; }
  T &operator()(size_t X, size_t Y, size_t Z) const { return data[X][Y][Z]; }

  size_t extent(int dim) const
  {
    if (dim == 0) return data.size();
    if (dim == 1 && !data.empty()) return data[0].size();
    if (dim == 2 && !data.empty() && !data[0].empty()) return data[0][0].size();
    return 0;
  }

 private:
  std::vector<std::vector<std::vector<T>>> data;
};

using vectorReal = Vector<float>;
using vectorInt = Vector<int>;
using vectorDouble = Vector<double>;

using arrayInt = Array2D<int>;
using arrayReal = Array2D<float>;
using arrayDouble = Array2D<double>;

using array3DInt = Array3D<int>;
using array3DReal = Array3D<float>;
using array3DDouble = Array3D<double>;

#endif  // USE_VECTOR

#ifdef USE_KOKKOS

#include "Kokkos_Core_fwd.hpp"

#ifdef ENABLE_HIP
#define __HIP_PLATFORM_AMD__ 1
#endif

#include <Kokkos_Core.hpp>
// #define MemSpace Kokkos::SharedSpace
//  #define MemSpace Kokkos::CudaUVMSpace

#ifdef ENABLE_CUDA
#define MemSpace Kokkos::SharedSpace
using Layout = Kokkos::LayoutLeft;
#else
#define MemSpace Kokkos::HostSpace
using Layout = Kokkos::LayoutRight;
#endif

typedef Kokkos::View<int *, Layout, MemSpace> vectorInt;
typedef Kokkos::View<float *, Layout, MemSpace> vectorReal;
typedef Kokkos::View<double *, Layout, MemSpace> vectorDouble;

typedef Kokkos::View<int **, Layout, MemSpace> arrayInt;
typedef Kokkos::View<float **, Layout, MemSpace> arrayReal;
typedef Kokkos::View<double **, Layout, MemSpace> arrayDouble;

typedef Kokkos::View<int ***, Layout, MemSpace> array3DInt;
typedef Kokkos::View<float ***, Layout, MemSpace> array3DReal;
typedef Kokkos::View<double ***, Layout, MemSpace> array3DDouble;

#endif  // USE_KOKKOS

template <class T>
T allocateVector(int n1)
{
#ifdef PRINT_ALLOC_INFO
  std::cout << "allocate vector of size " << n1 << std::endl;
#endif
  T vect(KOKKOSNAME n1);
  return vect;
}
template <class T>
T allocateVector(int n1, const char *name)
{
#ifdef PRINT_ALLOC_INFO
  std::cout << "allocate vector: " << name << " of size: " << n1 << std::endl;
#endif
  T vect(KOKKOSNAME n1);
  return vect;
}
template <class T>
T allocateArray2D(int n1, int n2)
{
#ifdef PRINT_ALLOC_INFO
  std::cout << "allocate array of size " << n1 << ", " << n2 << std::endl;
#endif
  T array(KOKKOSNAME n1, n2);
  return array;
}
template <class T>
T allocateArray2D(int n1, int n2, const char *name)
{
#ifdef PRINT_ALLOC_INFO
  std::cout << "allocate array : " << name << " of size: (" << n1 << ", " << n2
            << ")" << std::endl;
#endif
  T array(KOKKOSNAME n1, n2);
  return array;
}
template <class T>
T allocateArray3D(int n1, int n2, int n3)
{
#ifdef PRINT_ALLOC_INFO
  std::cout << "allocate array of size " << n1 << ", " << n2 << ", " << n3
            << std::endl;
#endif
  T array(KOKKOSNAME n1, n2, n3);
  return array;
}

// to be placed somewhere else
template <typename T, typename... Args>
void printJMatrix(const int &element, T &J, string matrixname, Args... args)
{
  if (element < 2)
  {
    (cout << ... << args) << '\n';
    printf("%s at element %d\n", matrixname.c_str(), element);
#ifdef USE_SHIVA
    for (int l = 0; l < 3; l++)
      printf("%f, %f, %f\n", J(l, 0), J(l, 1), J(l, 2));
#else
    for (int l = 0; l < 3; l++)
      printf("%f, %f, %f\n", J[l][0], J[l][1], J[l][2]);
#endif
  }
}

template <typename T>
void printBMatrix(const int &element, T &B)
{
  // print B matrix at the first element
  if (element < 2)
  {
    printf("\nB matrix at element %d\n", element);
    printf("%f, %f, %f\n", B[0], B[5], B[4]);
    printf("%f, %f, %f\n", B[5], B[1], B[3]);
    printf("%f, %f, %f\n\n", B[4], B[3], B[2]);
  }
}

// float jacobianTime=0;
// float detJTime=0;
// float massMatrixTime=0;
// float BTime=0;
// float gradPhiBGradPhiTime=0;
// float stiffnessTime=0;

#define timewatch(timepoint)                                \
  chrono::time_point<std::chrono::system_clock> timepoint = \
      chrono::system_clock::now();
#define accumtime(accumulatedtime, starttime) \
  accumulatedtime += (chrono::system_clock::now() - starttime).count();

#endif  // DATATYPE_HPP_
