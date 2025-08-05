//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "SEMproxy.hpp"
#include "dataType.hpp"
#ifdef USE_EZV
#include "ezvLauncher.hpp"
#endif // USE_EZV
#include <iomanip>
#include <iostream>
#include <sstream>

SEMproxy::SEMproxy(int argc, char *argv[]) {
  int ex = (cmdOptionExists(argv, argv + argc, "-ex"))
               ? std::stoi(getCmdOption(argv, argv + argc, "-ex"))
               : 100;
  int ey = (cmdOptionExists(argv, argv + argc, "-ey"))
               ? std::stoi(getCmdOption(argv, argv + argc, "-ey"))
               : ex;
  ;
  int ez = (cmdOptionExists(argv, argv + argc, "-ez"))
               ? std::stoi(getCmdOption(argv, argv + argc, "-ez"))
               : ex;
  ;

  float lx = (cmdOptionExists(argv, argv + argc, "-lx"))
                 ? std::stof(getCmdOption(argv, argv + argc, "-lx"))
                 : 2000;
  float ly = (cmdOptionExists(argv, argv + argc, "-ly"))
                 ? std::stof(getCmdOption(argv, argv + argc, "-ly"))
                 : lx;
  float lz = (cmdOptionExists(argv, argv + argc, "-lz"))
                 ? std::stof(getCmdOption(argv, argv + argc, "-lz"))
                 : lx;

  Mesh simpleMesh(ex, ey, ez, lx, ly, lz, order);
  myMesh = simpleMesh;
}

void SEMproxy::run() {
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  for (int indexTimeSample = 0; indexTimeSample < myNumSamples;
       indexTimeSample++) {
    startComputeTime = system_clock::now();
    mySolver.computeOneStep(indexTimeSample, myTimeStep, i1, i2, myRHSTerm, pnGlobal,
                            rhsElement, rhsWeights);
    totalComputeTime += system_clock::now() - startComputeTime;

    startOutputTime = system_clock::now();

    FENCE

    // if (indexTimeSample == 0)
    // {
    //   for (int i = 0; i < (10*3+1) *(10*3+1) * (10*3+1)  ; i++)
    //   {
    //     cout << pnGlobal(i1, i) << " ";
    //     if (i % 10 == 0) cout << endl;
    //   }
    // }
    if (indexTimeSample == 1400)
    {
      VECTOR_REAL_VIEW slice = extractXYSlice(Kokkos::subview(pnGlobal, i1, Kokkos::ALL), 10*3+1, (10*3+1)/2);
      saveSlice(slice, 10*3+1, "slice.dat");
    }

    swap(i1, i2);
    totalOutputTime += system_clock::now() - startOutputTime;
  }

  float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                            .time_since_epoch()
                            .count();
  float outputtime_ms =
      time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();

  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
       << endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
       << endl;
  cout << "------------------------------------------------ " << endl;
}

// Initialize arrays
void SEMproxy::init_arrays() {
  cout << "Allocate host memory for source and pressure values ..." << endl;
  myRHSTerm =
      allocateArray2D<arrayReal>(myNumberOfRHS, myNumSamples, "RHSTerm");
  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(myNumberOfRHS, myMesh.getNumberOfPointsPerElement(), "RHSWeight");
  pnGlobal =
      allocateArray2D<arrayReal>(myMesh.getNumberOfNodes(), 2, "pnGlobal");
}

// Initialize sources
void SEMproxy::init_source() {
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // set number of rhs and location
  myRHSLocation(0, 0) = myMesh.domainSize(0) / 2;
  myRHSLocation(0, 1) = myMesh.domainSize(1) / 2;
  myRHSLocation(0, 2) = myMesh.domainSize(2) / 2;
  cout << "\nSource 1 location: " << myRHSLocation(0, 0) << ", "
       << myRHSLocation(0, 1) << ", " << myRHSLocation(0, 2) << endl;
  cout << "Corresponding to element id "
       << myMesh.elementFromCoordinate(myRHSLocation(0,0), myRHSLocation(0,1),myRHSLocation(0,2)) << endl;
  for (int i = 0; i < 1; i++) {
    // extract element number for current rhs
    rhsElement[i] = myMesh.elementFromCoordinate(myRHSLocation(i,0), myRHSLocation(i,1),myRHSLocation(i,2));
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(myNumSamples, myTimeStep, f0, sourceOrder);
  for (int j = 0; j < myNumSamples; j++) {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }
  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;
  // Setting the weight for source ponderation on node.
  for (int i = 0; i < myNumberOfRHS; i++)
  {
    for (int j = 0; j < myMesh.getNumberOfPointsPerElement(); j++)
    {
      rhsWeights(i, j) = 1.0/ myMesh.getNumberOfPointsPerElement();
    }
  }
}

std::string formatSnapshotFilename(int id, int width = 5) {
  std::ostringstream oss;
  // oss << "snapshot" << std::setw(width) << std::setfill('0') << id;
  oss << "snapshot" << id;
  return oss.str();
}

#ifndef USE_KOKKOS
VECTOR_REAL_VIEW SEMproxy::extractXYSlice(const VECTOR_REAL_VIEW& array, int size, int z)
{
  int expected_size = size * size * size;

  // Calculate slice parameters
  int slice_size = size * size;
  int start_index = z * slice_size;

  // Create output view
  VECTOR_REAL_VIEW xy_slice("xy_slice", slice_size);

  // Extract the slice using parallel_for
  // Kokkos::parallel_for("extract_xy_slice", slice_size, KOKKOS_LAMBDA(int i) {
  LOOPHEAD(slice_size, i)
      xy_slice(i) = array_1d(start_index + i);
  LOOPEND
  // });

  return xy_slice;
}

#else // USE_KOKKOS
VECTOR_REAL_VIEW SEMproxy::extractXYSlice(const VECTOR_REAL_VIEW& array, int size, int z) {
    // Validate inputs
    if (z < 0 || z >= size) {
        Kokkos::abort("Z index out of bounds");
    }

    int slice_size = size * size;
    int start_index = z * slice_size;
    int end_index = start_index + slice_size;

    // Create subview (zero-copy operation)
    return Kokkos::subview(array, Kokkos::make_pair(start_index, end_index));
}
#endif // USE_KOKKOS

/**
 * Save slice in gnuplot matrix format (default - best for gnuplot)
 * Format: space-separated matrix with blank lines between rows for 3D plotting
 */
void SEMproxy::saveSlice(const VECTOR_REAL_VIEW& host_slice,
               int size, const std::string& filepath) {
    std::ofstream file(filepath);

    file << std::fixed << std::setprecision(6);

    // Standard matrix format - works with 'plot "file" matrix with image'
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            file << host_slice(y * size + x);
            if (x < size - 1) file << " ";
        }
        file << "\n";
    }

    file.close();
}
