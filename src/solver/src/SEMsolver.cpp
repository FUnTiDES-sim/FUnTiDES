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
#include <cstdlib>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // USE_KOKKOS

void SEMsolver::computeFEInit(Mesh mesh_in) {
  myMesh = mesh_in;
  allocateFEarrays();
  initFEarrays();
}

void SEMsolver::computeOneStep(const int &timeSample, const int &i1,
                               const int &i2, const ARRAY_REAL_VIEW &rhsTerm,
                               const ARRAY_REAL_VIEW &pnGlobal,
                               const VECTOR_INT_VIEW &rhsElement) {
  FENCE
  resetGlobalVectors(myMesh.getNumberOfNodes());
  applyRHSTerm(timeSample, i2, rhsTerm, rhsElement, pnGlobal);
  FENCE
  computeElementContributions(i2, pnGlobal);
  FENCE
  updatePressureField(i1, i2, pnGlobal);
  FENCE
}

void SEMsolver::resetGlobalVectors(int numNodes) {
  LOOPHEAD(numNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND
}

void SEMsolver::applyRHSTerm(int timeSample, int i2,
                             const ARRAY_REAL_VIEW &rhsTerm,
                             const VECTOR_INT_VIEW &rhsElement,
                             const ARRAY_REAL_VIEW &pnGlobal) {
  float const dt2 = myTimeStep * myTimeStep;
  LOOPHEAD(rhsElement.size(), i)
  int nodeRHS = myMesh.globalNodeIndex(rhsElement[i], 0, 0, 0);
  float scale =
      dt2 * myMesh.getModel(rhsElement[i]) * myMesh.getModel(rhsElement[i]);
  pnGlobal(nodeRHS, i2) += scale * rhsTerm(i, timeSample);
  LOOPEND
}

void SEMsolver::computeElementContributions(int i2,
                                            const ARRAY_REAL_VIEW &pnGlobal) {
  MAINLOOPHEAD(myMesh.getNumberOfElements(), elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= myMesh.getNumberOfElements())
    return;

  float massMatrixLocal[myMesh.getNumberOfPointsPerElement()] = {0};
  float pnLocal[myMesh.getNumberOfPointsPerElement()] = {0};
  float Y[myMesh.getNumberOfPointsPerElement()] = {0};

  int dim = myMesh.getOrder() + 1;
  for (int i = 0; i < myMesh.getNumberOfPointsPerElement(); ++i) {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const globalIdx = myMesh.globalNodeIndex(elementNumber, x, y, z);
    pnLocal[i] = pnGlobal(globalIdx, i2);
  }

#ifdef USE_SEMCLASSIC
  myQkIntegrals.computeMassMatrixAndStiffnessVector(
      elementNumber, myMesh.getOrder(), myMesh.getNumberOfPointsPerElement(),
      globalNodesCoordsX, globalNodesCoordsY, globalNodesCoordsZ, weights,
      derivativeBasisFunction1D, massMatrixLocal, pnLocal, Y);
#else
  myQkIntegrals.computeMassMatrixAndStiffnessVector(
      elementNumber, myMesh.getNumberOfPointsPerElement(), globalNodesCoordsX,
      globalNodesCoordsY, globalNodesCoordsZ, massMatrixLocal, pnLocal, Y);
#endif

  auto const inv_model2 =
      1.0f / (myMesh.getModel(elementNumber) * myMesh.getModel(elementNumber));
  for (int i = 0; i < myMesh.getNumberOfPointsPerElement(); ++i) {
    int x = i % dim;
    int z = (i / dim) % dim;
    int y = i / (dim * dim);
    int const gIndex = myMesh.globalNodeIndex(elementNumber, x, y, z);
    massMatrixLocal[i] *= inv_model2;
    ATOMICADD(massMatrixGlobal[gIndex], massMatrixLocal[i]);
    ATOMICADD(yGlobal[gIndex], Y[i]);
  }

  MAINLOOPEND
}

void SEMsolver::updatePressureField(int i1, int i2,
                                    const ARRAY_REAL_VIEW &pnGlobal) {

  float const dt2 = myTimeStep * myTimeStep;
  LOOPHEAD(myMesh.getNumberOfNodes(), I)
  pnGlobal(I, i1) = 2 * pnGlobal(I, i2) - pnGlobal(I, i1) -
                    dt2 * yGlobal[I] / massMatrixGlobal[I];
  pnGlobal(I, i1) *= spongeTaperCoeff(I);
  pnGlobal(I, i2) *= spongeTaperCoeff(I);
  LOOPEND
}

void SEMsolver::outputPnValues(Mesh mesh, const int &indexTimeStep, int &i1,
                               int &myElementSource,
                               const ARRAY_REAL_VIEW &pnGlobal) {
  // writes debugging ascii file.
  if (indexTimeStep % 50 == 0) {
    cout << "TimeStep=" << indexTimeStep
         << ";  pnGlobal @ elementSource location " << myElementSource
         << " after computeOneStep = "
         << pnGlobal(myMesh.globalNodeIndex(myElementSource, 0, 0, 0), i1)
         << endl;
#ifdef SEM_SAVE_SNAPSHOTS
    mesh.saveSnapShot(indexTimeStep, i1, pnGlobal);
#endif // SEM_SAVE_SNAPSHOTS
  }
}

void SEMsolver::initFEarrays() {

  // get quadrature points
#ifdef USE_SEMCLASSIC
  myQkBasis.gaussLobattoQuadraturePoints(myMesh.getOrder(), quadraturePoints);
  // get gauss-lobatto weights
  myQkBasis.gaussLobattoQuadratureWeights(myMesh.getOrder(), weights);
  // get basis function and corresponding derivatives
  myQkBasis.getDerivativeBasisFunction1D(myMesh.getOrder(), quadraturePoints,
                                         derivativeBasisFunction1D);
#endif // USE_SEMCLASSIC
}

//************************************************************************
//  Allocate arrays for the solver
//  This function allocates all arrays needed for the solver
//  It allocates arrays for global nodes, global coordinates, and sponge
//  It also allocates arrays for the mass matrix and the global pressure field
//************************************************************************
void SEMsolver::allocateFEarrays() {
  int nbQuadraturePoints = (myMesh.getOrder() + 1) * (myMesh.getOrder() + 1) *
                           (myMesh.getOrder() + 1);

  cout << "Allocate host memory for arrays in the solver ..." << endl;

#ifdef USE_SEMCLASSIC
  quadraturePoints =
      allocateVector<VECTOR_REAL_VIEW>(nbQuadraturePoints, "quadraturePoints");
  weights = allocateVector<VECTOR_REAL_VIEW>(nbQuadraturePoints, "weights");
  derivativeBasisFunction1D = allocateArray2D<ARRAY_REAL_VIEW>(
      myMesh.getNumberOfNodes(), nbQuadraturePoints,
      "derivativeBasisFunction1D");
#endif // USE_SEMCLASSIC

  // shared arrays
  massMatrixGlobal = allocateVector<VECTOR_REAL_VIEW>(myMesh.getNumberOfNodes(),
                                                      "massMatrixGlobal");
  yGlobal =
      allocateVector<VECTOR_REAL_VIEW>(myMesh.getNumberOfNodes(), "yGlobal");

  // sponge allocation
  spongeTaperCoeff = allocateVector<VECTOR_REAL_VIEW>(myMesh.getNumberOfNodes(),
                                                      "spongeTaperCoeff");
}
