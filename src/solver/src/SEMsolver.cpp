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
  // myInfo = myInfo;
  myMesh = mesh;
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

  LOOPHEAD(myInfo.numberOfNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND

  // update pnGLobal with right hade side
  LOOPHEAD(myInfo.myNumberOfRHS, i)
  int nodeRHS = globalNodesList(rhsElement[i], 0);
  pnGlobal(nodeRHS, i2) += myInfo.myTimeStep * myInfo.myTimeStep *
                           model[rhsElement[i]] * model[rhsElement[i]] *
                           rhsTerm(i, timeSample);
  LOOPEND
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
#endif // USE_SEMCLASSIC

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
#endif // USE_SHIVA

#ifdef USE_SEMGEOS
    myQkIntegrals.computeMassMatrixAndStiffnessVector(
        elementNumber, nPointsPerElement, globalNodesCoordsX,
        globalNodesCoordsY, globalNodesCoordsZ, massMatrixLocal, pnLocal, Y);
#endif // USE_SEMGEOS

    // compute global mass Matrix and global stiffness vector
    for (int i = 0; i < nPointsPerElement; i++) {
      int gIndex = globalNodesList(elementNumber, i);
      massMatrixLocal[i] /= (model[elementNumber] * model[elementNumber]);
      ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
      ATOMICADD(yGlobal[gIndex], Y[i]);
    }
  }
  MAINLOOPEND

  // update pressure
  LOOPHEAD(myInfo.numberOfInteriorNodes, i)
  int I = listOfInteriorNodes[i];
  pnGlobal(I, i1) =
      2 * pnGlobal(I, i2) - pnGlobal(I, i1) -
      myInfo.myTimeStep * myInfo.myTimeStep * yGlobal[I] / massMatrixGlobal[I];
  LOOPEND

  // update pressure for sponge
  spongeUpdate(pnGlobal, i1, i2);

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
#endif // SEM_SAVE_SNAPSHOTS
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

  // get model
  mesh.getModel(myInfo.numberOfElements, model);

  // get minimal wavespeed
  double min;
  auto model_ = this->model; // Avoid implicit capture
#ifdef USE_KOKKOS
  Kokkos::parallel_reduce(
      "vMinFind", myInfo.numberOfElements,
      KOKKOS_LAMBDA(const int &e, double &lmin) {
        double val = model_[e];
        if (val < lmin)
          lmin = val;
      },
      Kokkos::Min<double>(min));
  vMin = min;
#else
  vMin = 1500;
#endif // USE_KOKKOS

  // get quadrature points
#ifdef USE_SEMCLASSIC
  myQkBasis.gaussLobattoQuadraturePoints(order, quadraturePoints);
  // get gauss-lobatto weights
  myQkBasis.gaussLobattoQuadratureWeights(order, weights);
  // get basis function and corresponding derivatives
  myQkBasis.getDerivativeBasisFunction1D(order, quadraturePoints,
                                         derivativeBasisFunction1D);
#endif // USE_SEMCLASSIC

  // Sponge boundaries
  initSpongeValues(mesh, myInfo, vMin, 0.001);
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
  listOfDampingNodes = allocateVector<vectorInt>(myInfo.numberOfDampingNodes,
                                                 "listOfDampingNodes");
  cout << "### Allocating " << myInfo.numberOfSpongeNodes << " sponge nodes"
       << endl;

  // global coordinates
  globalNodesCoordsX = allocateArray2D<arrayReal>(
      myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsX");
  globalNodesCoordsY = allocateArray2D<arrayReal>(
      myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsY");
  globalNodesCoordsZ = allocateArray2D<arrayReal>(
      myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsZ");

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

  // sponge allocation
  spongeTaperCoeff =
      allocateVector<vectorReal>(myInfo.numberOfNodes, "spongeTaperCoeff");
}

void SEMsolver::initSpongeValues(Mesh &mesh, SEMinfo &myInfo, const float vMin,
                                 const float r) {
  // Init all taper to 1 (default value)
  for (int i = 0; i < myInfo.numberOfNodes; i++) {
    spongeTaperCoeff(i) = 1;
  }

  int n = 0;
  int alpha = -0.0001;
  int spongeSize = mesh.getSpongeSize();
  int nx = mesh.getNx();
  int ny = mesh.getNy();
  int nz = mesh.getNz();

  // Update X boundaries
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      // lower x
      for (int i = 0; i <= spongeSize; i++) {
        n = mesh.ijktoI(i, j, k);
        int value = spongeSize - i;
        spongeTaperCoeff(n) =
            std::exp(alpha * static_cast<double>(value * value));
      }
      // upper x
      for (int i = nx - spongeSize - 1; i < nx; i++) {
        n = mesh.ijktoI(i, j, k);
        int value = spongeSize - (nx - i);
        spongeTaperCoeff(n) =
            std::exp(alpha * static_cast<double>(value * value));
      }
    }
  }

  // Update Y boundaries
  for (int k = 0; k < nz; k++) {
    for (int i = 0; i < nx; i++) {
      // lower y
      for (int j = 0; j <= spongeSize; j++) {
        n = mesh.ijktoI(i, j, k);
        int value = spongeSize - j;
        spongeTaperCoeff(n) =
            std::exp(alpha * static_cast<double>(value * value));
      }
      // upper y
      for (int j = ny - spongeSize - 1; j < ny; j++) {
        n = mesh.ijktoI(i, j, k);
        int value = spongeSize - (ny - j);
        spongeTaperCoeff(n) =
            std::exp(alpha * static_cast<double>(value * value));
      }
    }
  }

  // Update Z boundaries
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      // lower z
      for (int k = 0; k <= spongeSize; k++) {
        n = mesh.ijktoI(i, j, k);
        int value = spongeSize - k;
        spongeTaperCoeff(n) =
            std::exp(alpha * static_cast<double>(value * value));
      }
      // upper z
      for (int k = nz - spongeSize - 1; k < nz; k++) {
        n = mesh.ijktoI(i, j, k);
        int value = spongeSize - (nz - k);
        spongeTaperCoeff(n) =
            std::exp(alpha * static_cast<double>(value * value));
      }
    }
  }
}

void SEMsolver::spongeUpdate(const arrayReal &pnGlobal, const int i1,
                             const int i2) {
  for (int i = 0; i < myInfo.numberOfNodes; i++) {
    pnGlobal(i, i1) *= spongeTaperCoeff(i);
    pnGlobal(i, i2) *= spongeTaperCoeff(i);
  }
}
