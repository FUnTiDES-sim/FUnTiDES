//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "SEMsolver.hpp"
#ifdef USE_EZV
#include "ezvLauncher.hpp"
#include <cstdlib>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // USE_KOKKOS
#endif // USE_EZV

void SEMsolver::computeFEInit(SEMinfo &myInfo, Mesh mesh) {
  order = myInfo.myOrderNumber;
  allocateFEarrays(myInfo);
  initFEarrays(myInfo, mesh);
}

// compute one step of the time dynamic wave equation solver
void SEMsolver::computeOneStep(const int &timeSample, const int &order,
                               const int &nPointsPerElement, const int &i1,
                               const int &i2, SEMinfo &myInfo,
                               const arrayReal &rhsTerm,
                               arrayReal const &pnGlobal,
                               const vectorInt &rhsElement) {
#ifdef USE_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  LOOPHEAD(myInfo.numberOfNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND

#ifdef USE_CALIPER
  CALI_MARK_BEGIN("updateNodeRHS");
#endif
  // update pnGLobal with right hade side
  LOOPHEAD(myInfo.myNumberOfRHS, i)
  int nodeRHS = globalNodesList(rhsElement[i], 0);
  pnGlobal(nodeRHS, i2) += myInfo.myTimeStep * myInfo.myTimeStep *
                           model[rhsElement[i]] * model[rhsElement[i]] *
                           rhsTerm(i, timeSample);
  LOOPEND
#ifdef USE_CALIPER
  CALI_MARK_END("updateNodeRHS");
  CALI_MARK_BEGIN("mainloop");
#endif
  // start main parallel section
  MAINLOOPHEAD(myInfo.numberOfElements, elementNumber)

  if (elementNumber < myInfo.numberOfElements) {
    float massMatrixLocal[ROW];
    float pnLocal[ROW];
    float Y[ROW];

    // get pnGlobal to pnLocal
    for (int i = 0; i < nPointsPerElement; i++) {
      int localToGlobal = globalNodesList(elementNumber, i);
      pnLocal[i] = pnGlobal(localToGlobal, i2);
    }

#ifdef USE_SEMCLASSIC
    myQkIntegrals.computeMassMatrixAndStiffnessVector(
        elementNumber, order, nPointsPerElement, globalNodesCoordsX,
        globalNodesCoordsY, globalNodesCoordsZ, weights,
        derivativeBasisFunction1D, massMatrixLocal, pnLocal, Y);
#endif

#ifdef USE_SEMOPTIM
    constexpr int ORDER = SEMinfo::myOrderNumber;
    myQkIntegrals.computeMassMatrixAndStiffnessVector<ORDER>(
        elementNumber, nPointsPerElement, globalNodesCoordsX,
        globalNodesCoordsY, globalNodesCoordsZ, massMatrixLocal, pnLocal, Y);
#endif

#ifdef USE_SHIVA
    constexpr int ORDER = SEMinfo::myOrderNumber;
    myQkIntegrals.computeMassMatrixAndStiffnessVector<ORDER>(
        elementNumber, nPointsPerElement, globalNodesCoordsX,
        globalNodesCoordsY, globalNodesCoordsZ, massMatrixLocal, pnLocal, Y);
#endif

    // compute global mass Matrix and global stiffness vector
    for (int i = 0; i < nPointsPerElement; i++) {
      int gIndex = globalNodesList(elementNumber, i);
      massMatrixLocal[i] /= (model[elementNumber] * model[elementNumber]);
      ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
      ATOMICADD(yGlobal[gIndex], Y[i]);
    }
  }
  MAINLOOPEND

#ifdef USE_CALIPER
  CALI_MARK_END("mainloop");
  CALI_MARK_BEGIN("pressureloop");
#endif // USE_CALIPER

#ifdef USE_EZV
  //  Kokkos::View<float *, Kokkos::CudaSpace> ezv_device_data("EZV Device",
  //  pnGlobal.size());
  auto ezv_host_data = Kokkos::create_mirror(pnGlobal);
  Kokkos::fence();
#endif

  // update pressure
  LOOPHEAD(myInfo.numberOfInteriorNodes, i)
  int I = listOfInteriorNodes[i];
  pnGlobal(I, i1) =
      2 * pnGlobal(I, i2) - pnGlobal(I, i1) -
      myInfo.myTimeStep * myInfo.myTimeStep * yGlobal[I] / massMatrixGlobal[I];
  LOOPEND

#ifdef USE_EZV
  Kokkos::fence();
  Kokkos::deep_copy(ezv_host_data, pnGlobal);
  auto nb_x = pnGlobal.extent(0);
  auto nb_y = pnGlobal.extent(1);
  Kokkos::fence();
  float *ezv_data = (float *)malloc((nb_x) * sizeof(float));
  // copy kokkos data into heap
  for (auto i = 0; i < nb_x; i++) {
    auto idx = i;
    ezv_data[idx] = ezv_host_data(i, i1);
  }

  ezv_thr_push_data_colors(get_ezv_ctx()[0], ezv_data);
#endif // USE_EZV

#ifdef USE_CALIPER
  CALI_MARK_END("pressureloop");
#endif // USE_CALIPER

  FENCE
}

void SEMsolver::outputPnValues(Mesh mesh, const int &indexTimeStep, int &i1,
                               int &myElementSource,
                               const arrayReal &pnGlobal) {
  // writes debugging ascii file.
  if (indexTimeStep % 50 == 0) {
    cout << "TimeStep=" << indexTimeStep
         << ";  pnGlobal @ elementSource location " << myElementSource
         << " after computeOneStep = "
         << pnGlobal(globalNodesList(myElementSource, 0), i1) << endl;
#ifdef SEM_SAVE_SNAPSHOTS
    mesh.saveSnapShot(indexTimeStep, i1, pnGlobal);
#endif
  }
}

void SEMsolver::initFEarrays(SEMinfo &myInfo, Mesh mesh) {
  // interior elements
  mesh.globalNodesList(myInfo.numberOfElements, globalNodesList);
  mesh.getListOfInteriorNodes(myInfo.numberOfInteriorNodes,
                              listOfInteriorNodes);
  // mesh coordinates
  mesh.nodesCoordinates(globalNodesCoordsX, globalNodesCoordsZ,
                        globalNodesCoordsY);
#endif
  // get model
  mesh.getModel(myInfo.numberOfElements, model);
  // get quadrature points
#ifdef USE_SEMCLASSIC
  myQkBasis.gaussLobattoQuadraturePoints(order, quadraturePoints);
  // get gauss-lobatto weights
  myQkBasis.gaussLobattoQuadratureWeights(order, weights);
  // get basis function and corresponding derivatives
  myQkBasis.getDerivativeBasisFunction1D(order, quadraturePoints,
                                         derivativeBasisFunction1D);
#endif
}

void SEMsolver::allocateFEarrays(SEMinfo &myInfo) {
  int nbQuadraturePoints = (order + 1) * (order + 1) * (order + 1);
  // interior elements
  cout << "Allocate host memory for arrays in the solver ..." << endl;
  globalNodesList = allocateArray2D<arrayInt>(myInfo.numberOfElements,
                                              myInfo.numberOfPointsPerElement,
                                              "globalNodesList");
  listOfInteriorNodes = allocateVector<vectorInt>(myInfo.numberOfInteriorNodes,
                                                  "listOfInteriorNodes");

  // global coordinates
  globalNodesCoordsX=allocateArray2D< arrayReal >( myInfo.numberOfElements,
    nbQuadraturePoints, "globalNodesCoordsX");
  globalNodesCoordsY=allocateArray2D< arrayReal >( myInfo.numberOfElements,
    nbQuadraturePoints, "globalNodesCoordsY");
  globalNodesCoordsZ=allocateArray2D< arrayReal >( myInfo.numberOfElements,
    nbQuadraturePoints, "globalNodesCoordsZ");

  model = allocateVector<vectorReal>(myInfo.numberOfElements, "model");

  quadraturePoints =
      allocateVector<vectorDouble>(order + 1, "quadraturePoints");

  weights = allocateVector<vectorDouble>(order + 1, "weights");

  derivativeBasisFunction1D = allocateArray2D<arrayDouble>(
      order + 1, order + 1, "derivativeBasisFunction1D");

  // shared arrays
  massMatrixGlobal =
      allocateVector<vectorReal>(myInfo.numberOfNodes, "massMatrixGlobal");
  yGlobal = allocateVector<vectorReal>(myInfo.numberOfNodes, "yGlobal");
}
