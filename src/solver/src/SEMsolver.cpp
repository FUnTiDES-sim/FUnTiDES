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
  spongeUpdate(myMesh, myInfo, pnGlobal, i1, i2);

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
  // sponge element id
  mesh.getListOfSpongeNodes(listOfSpongeNodes);

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

  getSpongeValues(mesh, myInfo, vMin, 0.00000001);
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
  listOfSpongeNodes = allocateVector<vectorInt>(myInfo.numberOfSpongeNodes,
                                                "listOfSpongeNodes");
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
  spongeTaperCoeff = allocateArray2D<arrayReal>(
      myInfo.numberOfElements, nbQuadraturePoints, "spongeTaperCoeff");
}

void SEMsolver::getSpongeValues(Mesh &mesh, SEMinfo &myInfo, const float vMin,
                                const float r) {

  MinMax3D globalMinMax;

  // Kokkos not allowing 'this' capture in lambda
  auto coordsX = globalNodesCoordsX;
  auto coordsY = globalNodesCoordsY;
  auto coordsZ = globalNodesCoordsZ;

  // Collects min and max value for each coordinate x y z
  // at node level
#ifdef USE_KOKKOS
  Kokkos::parallel_reduce(
      "Sponge Value MinMax", myInfo.numberOfElements,
      KOKKOS_LAMBDA(int e, MinMax3D &local) {
        // Loop over nodes in one element
        for (int i = 0; i < myInfo.numberOfPointsPerElement; i++) {
          const auto globalX = coordsX(e, i);
          const auto globalY = coordsY(e, i);
          const auto globalZ = coordsZ(e, i);

          if (globalX < local.min_x)
            local.min_x = globalX;
          if (globalX > local.max_x)
            local.max_x = globalX;

          if (globalY < local.min_y)
            local.min_y = globalY;
          if (globalY > local.max_y)
            local.max_y = globalY;

          if (globalZ < local.min_z)
            local.min_z = globalZ;
          if (globalZ > local.max_z)
            local.max_z = globalZ;
        }
      },
      globalMinMax);
#else
  throw std::logic_error(
      "Function MinMax not yet implemented outside Kokkos env.");
#endif // USE_KOKKOS

  // Compute sponge distance for every elements' nodes
  LOOPHEAD(myInfo.numberOfElements, e)
  for (int i = 0; i < myInfo.numberOfPointsPerElement; i++) {
    float taper = 0;
    float dist = 0;
    const auto nodeX = coordsX(e, i);
    const auto nodeY = coordsY(e, i);
    const auto nodeZ = coordsZ(e, i);

    float distXmin = nodeX - globalMinMax.min_x;
    float distYmin = nodeY - globalMinMax.min_y;

    float distXmax = nodeX - globalMinMax.max_x;
    float distYmax = nodeY - globalMinMax.max_y;
    float distZmax = nodeZ - globalMinMax.max_z;

    dist = std::min({distXmin, distXmax, distYmin, distYmax, distZmax});
    taper = std::exp(
        (((3 * vMin) / (2 * mesh.getSpongeSize())) * std::log(r) *
         std::pow((mesh.getSpongeSize() - dist) / mesh.getSpongeSize(), 2)) *
        myInfo.myTimeStep);

    spongeTaperCoeff(e, i) = taper;
  }
  LOOPEND
}

void SEMsolver::spongeUpdate(Mesh mesh, SEMinfo &myInfo,
                             const arrayReal &pnGlobal, const int i1,
                             const int i2) {
  LOOPHEAD(myInfo.numberOfSpongeNodes, i)
  int I = listOfSpongeNodes[i];
  int e, ii;
  int ret = myMesh.ItoEi(I, &e, &ii);
  assert(ret == 0);
  for (int pn : {i1}) {
    pnGlobal(I, pn) = 2 * pnGlobal(I, i2) - pnGlobal(I, i1) -
                      myInfo.myTimeStep * myInfo.myTimeStep * yGlobal[I] /
                          massMatrixGlobal[I] * spongeTaperCoeff(e, ii);
  }

  LOOPEND
}
