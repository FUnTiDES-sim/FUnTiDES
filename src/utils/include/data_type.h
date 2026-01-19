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
#include "data_type_vector.h"
#endif  // USE_VECTOR

#ifdef USE_KOKKOS
#include "data_type_kokkos.h"
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
