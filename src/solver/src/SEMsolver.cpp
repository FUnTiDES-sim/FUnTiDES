//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "SEMsolver.hpp"
#include "dataType.hpp"
#ifdef USE_EZV
#include "ezvLauncher.hpp"
#include <cstdlib>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // USE_KOKKOS
#endif // USE_EZV

void SEMsolver::computeFEInit(SEMinfo &myInfo_in, Mesh mesh) {
  this->myInfo = &myInfo_in;
  myMesh = mesh;
  order = myInfo_in.myOrderNumber;
  allocateFEarrays(myInfo_in);
  initFEarrays(myInfo_in, mesh);
}

void SEMsolver::computeOneStep(const int &timeSample, const int &order,
                               const int &nPointsPerElement, const int &i1,
                               const int &i2, SEMinfo &myInfo,
                               const arrayReal &rhsTerm, arrayReal &pnGlobal,
                               const vectorInt &rhsElement) {
  resetGlobalVectors(myInfo.numberOfNodes);
  applyRHSTerm(timeSample, i2, rhsTerm, rhsElement, myInfo, pnGlobal);
  FENCE
  computeElementContributions(order, nPointsPerElement, myInfo, i2, pnGlobal);
  FENCE
  updatePressureField(i1, i2, myInfo, pnGlobal);
  FENCE
}

void SEMsolver::resetGlobalVectors(int numNodes) {
  LOOPHEAD(numNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND
}

void SEMsolver::applyRHSTerm(int timeSample, int i2, const arrayReal &rhsTerm,
                             const vectorInt &rhsElement, SEMinfo &myInfo,
                             arrayReal &pnGlobal) 
{
  float const dt2 = myInfo.myTimeStep * myInfo.myTimeStep;
  LOOPHEAD(myInfo.myNumberOfRHS, i)
    int nodeRHS = globalNodesList(rhsElement[i], 0);
    float scale = dt2 * model[rhsElement[i]] * model[rhsElement[i]];
    pnGlobal(nodeRHS, i2) += scale * rhsTerm(i, timeSample);
  LOOPEND
}

void SEMsolver::computeElementContributions(int order, int nPointsPerElement,
                                            SEMinfo &myInfo, int i2,
                                            const arrayReal &pnGlobal) {
  MAINLOOPHEAD(myInfo.numberOfElements, elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= myInfo.numberOfElements)
    return;

  float massMatrixLocal[SEMinfo::nPointsPerElement];
  float pnLocal[SEMinfo::nPointsPerElement];
  float Y[SEMinfo::nPointsPerElement] = {0};

  for (int i = 0; i < nPointsPerElement; ++i) 
  {
    int const globalIdx = globalNodesList(elementNumber, i);
    pnLocal[i] = pnGlobal(globalIdx, i2);
  }

#if defined(USE_SEMCLASSIC)
  myQkIntegrals.computeMassMatrixAndStiffnessVector(
                elementNumber, order, nPointsPerElement, globalNodesCoordsX,
                globalNodesCoordsY, globalNodesCoordsZ, weights,
                derivativeBasisFunction1D, massMatrixLocal, pnLocal, Y);
#else
  myQkIntegrals.computeMassMatrixAndStiffnessVector( elementNumber, 
                                                     nPointsPerElement, 
                                                     globalNodesCoordsX, 
                                                     globalNodesCoordsY,
                                                     globalNodesCoordsZ, 
                                                     massMatrixLocal, 
                                                     pnLocal, 
                                                     Y);
#endif

  auto const inv_model2 = 1.0f / (model[elementNumber] * model[elementNumber]);
  for (int i = 0; i < SEMinfo::nPointsPerElement; ++i) 
  {
    int const gIndex = globalNodesList(elementNumber, i);
    massMatrixLocal[i] *= inv_model2;
    ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
    ATOMICADD(yGlobal[gIndex], Y[i]);
  }

  MAINLOOPEND
}

void SEMsolver::updatePressureField(int i1, int i2, SEMinfo &myInfo,
                                    arrayReal &pnGlobal) 
{

  float const dt2 = myInfo.myTimeStep * myInfo.myTimeStep;
  LOOPHEAD(myInfo.numberOfNodes, I)
    pnGlobal(I, i1) = 2 * pnGlobal(I, i2) - pnGlobal(I, i1) - dt2 * yGlobal[I] / massMatrixGlobal[I];
    pnGlobal(I, i1) *= spongeTaperCoeff(I);
    pnGlobal(I, i2) *= spongeTaperCoeff(I);
  LOOPEND
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
  initSpongeValues(mesh, myInfo);
  Kokkos::fence();
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
  
  #ifdef USE_SEMCLASSIC
  quadraturePoints =
      allocateVector<vectorDouble>(order + 1, "quadraturePoints");

  weights = allocateVector<vectorDouble>(order + 1, "weights");

  derivativeBasisFunction1D = allocateArray2D<arrayDouble>(
      order + 1, order + 1, "derivativeBasisFunction1D");
  #endif // USE_SEMCLASSIC

  // shared arrays
  massMatrixGlobal =
      allocateVector<vectorReal>(myInfo.numberOfNodes, "massMatrixGlobal");
  yGlobal = allocateVector<vectorReal>(myInfo.numberOfNodes, "yGlobal");

  // sponge allocation
  spongeTaperCoeff =
      allocateVector<vectorReal>(myInfo.numberOfNodes, "spongeTaperCoeff");
}

void SEMsolver::initSpongeValues(Mesh &mesh, SEMinfo &myInfo) {
  // Init all taper to 1 (default value)
  LOOPHEAD(myInfo.numberOfNodes, i)
  spongeTaperCoeff(i) = 1;
  LOOPEND

  // int n = 0;
  double alpha = -0.0001;
  int spongeSize = mesh.getSpongeSize();
  int nx = mesh.getNx();
  int ny = mesh.getNy();
  int nz = mesh.getNz();

  // Update X boundaries
  LOOPHEAD(nz, k)
  for (int j = 0; j < ny; j++) {
    // lower x
    for (int i = 0; i <= spongeSize; i++) {
      int n = mesh.ijktoI(i, j, k);
      double value = spongeSize - i;
      spongeTaperCoeff(n) = std::exp(alpha * static_cast<float>(value * value));
    }
    // upper x
    for (int i = nx - spongeSize - 1; i < nx; i++) {
      int n = mesh.ijktoI(i, j, k);
      double value = spongeSize - (nx - i);
      spongeTaperCoeff(n) = std::exp(alpha * static_cast<float>(value * value));
    }
  }
  LOOPEND

  // Update Y boundaries
  // for (int k = 0; k < nz; k++) {
  LOOPHEAD(nz, k)
  int n;
  for (int i = 0; i < nx; i++) {
    // lower y
    for (int j = 0; j <= spongeSize; j++) {
      n = mesh.ijktoI(i, j, k);
      double value = spongeSize - j;
      spongeTaperCoeff(n) =
          std::exp(alpha * static_cast<double>(value * value));
    }
    // upper y
    for (int j = ny - spongeSize - 1; j < ny; j++) {
      n = mesh.ijktoI(i, j, k);
      double value = spongeSize - (ny - j);
      spongeTaperCoeff(n) =
          std::exp(alpha * static_cast<double>(value * value));
    }
  }
  LOOPEND

  // Update Z boundaries
  LOOPHEAD(ny, j)
  int n;
  for (int i = 0; i < nx; i++) {
    // lower z
    for (int k = 0; k <= spongeSize; k++) {
      n = mesh.ijktoI(i, j, k);
      double value = spongeSize - k;
      spongeTaperCoeff(n) =
          std::exp(alpha * static_cast<double>(value * value));
    }
    // upper z
    for (int k = nz - spongeSize - 1; k < nz; k++) {
      n = mesh.ijktoI(i, j, k);
      double value = spongeSize - (nz - k);
      spongeTaperCoeff(n) =
          std::exp(alpha * static_cast<double>(value * value));
    }
  }
  LOOPEND

  FENCE
}

void SEMsolver::spongeUpdate(const arrayReal &pnGlobal, const int i1,
                             const int i2) {
  // for (int i = 0; i < myInfo->numberOfNodes; i++) {
  LOOPHEAD(myInfo->numberOfNodes, i)
  pnGlobal(i, i1) *= spongeTaperCoeff(i);
  pnGlobal(i, i2) *= spongeTaperCoeff(i);
  LOOPEND
}
// }
