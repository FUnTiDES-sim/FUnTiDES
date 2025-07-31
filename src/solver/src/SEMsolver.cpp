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

void SEMsolver::computeOneStep(const int &timeSample, const float dt,
                               const int &i1, const int &i2,
                               const ARRAY_REAL_VIEW &rhsTerm,
                               const ARRAY_REAL_VIEW &pnGlobal,
                               const VECTOR_INT_VIEW &rhsElement,
                               const ARRAY_REAL_VIEW &rhsWeights) {
  resetGlobalVectors(myMesh.getNumberOfNodes());
  FENCE
  applyRHSTerm(timeSample, dt, i2, rhsTerm, rhsElement, pnGlobal, rhsWeights);
  FENCE
  computeElementContributions(i2, pnGlobal);
  FENCE
  updatePressureField(dt, i1, i2, pnGlobal);
  FENCE
}

void SEMsolver::resetGlobalVectors(int numNodes) {
  LOOPHEAD(numNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND
}

void SEMsolver::applyRHSTerm(int timeSample, float dt, int i2,
                             const ARRAY_REAL_VIEW &rhsTerm,
                             const VECTOR_INT_VIEW &rhsElement,
                             const ARRAY_REAL_VIEW &pnGlobal,
                             const ARRAY_REAL_VIEW &rhsWeights) {
  float const dt2 = dt * dt;
  int nb_rhs_element = rhsElement.extent(0);
  LOOPHEAD(nb_rhs_element, i)
    for(int z = 0; z < myMesh.getOrder() + 1; z++)
    {
      for(int y = 0; y < myMesh.getOrder() + 1; y++)
      {
        for (int x = 0; x < myMesh.getOrder() + 1; x++)
        {
          int localNodeId = x + y * (myMesh.getOrder()+1) + z * (myMesh.getOrder() + 1) * (myMesh.getOrder() + 1);
          int nodeRHS = myMesh.globalNodeIndex(rhsElement[i], x, y, z);
          float scale =  dt2 * myMesh.getModelVpOnElement(rhsElement[i]) * myMesh.getModelVpOnElement(rhsElement[i]);
          float source = scale * rhsTerm(i, timeSample) * rhsWeights(i, localNodeId);
          yGlobal(nodeRHS) -= source;
        }
      }
    }
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
  float nodeCoords[myMesh.getNumberOfPointsPerElement()][3];
  int I = 0;
  for (int k = 0; k< myMesh.getOrder() + 1; k++) {
    for (int j = 0; j < myMesh.getOrder() + 1; j++) {
      for (int i = 0; i < myMesh.getOrder() + 1; i++) {
        int nodeIdx = myMesh.globalNodeIndex(elementNumber, i, j, k);
        nodeCoords[I][0] = myMesh.nodeCoord(nodeIdx, 0);
        nodeCoords[I][2] = myMesh.nodeCoord(nodeIdx, 2);
        nodeCoords[I][1] = myMesh.nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }
  myQkIntegrals.computeMassMatrixAndStiffnessVector(
      myMesh.getOrder(), myMesh.getNumberOfPointsPerElement(),
      nodeCoords, weights,
      derivativeBasisFunction1D, massMatrixLocal, pnLocal, Y);
#else
  // Init coordinates of element's corner
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, myMesh.getOrder()};
  for (int k : nodes_corner) {
    for (int j : nodes_corner) {
      for (int i : nodes_corner) {
        int nodeIdx = myMesh.globalNodeIndex(elementNumber, i, j, k);
        cornerCoords[I][0] = myMesh.nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = myMesh.nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = myMesh.nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }
  myQkIntegrals.computeMassMatrixAndStiffnessVector(
      elementNumber, myMesh.getNumberOfPointsPerElement(), cornerCoords, massMatrixLocal,
      pnLocal, Y);
#endif

  auto const inv_model2 =
      1.0f / (myMesh.getModelVpOnElement(elementNumber) * myMesh.getModelVpOnElement(elementNumber));
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

void SEMsolver::updatePressureField(float dt, int i1, int i2,
                                    const ARRAY_REAL_VIEW &pnGlobal) {

  float const dt2 = dt * dt;
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
  float sum = 0.0;
  Kokkos::parallel_reduce(
      "Sum pnGlobal", myMesh.getNumberOfNodes(),
      KOKKOS_LAMBDA(int i, float &local_sum) { local_sum += pnGlobal(i, i1); },
      sum);
  if (indexTimeStep % 50 == 0) {
    cout << "TimeStep=" << indexTimeStep
         << ";  pnGlobal @ elementSource location " << myElementSource
         << " after computeOneStep = "
         << pnGlobal(myMesh.globalNodeIndex(myElementSource, 0, 0, 0), i1)
         << " and sum pnGlobal is " << sum << endl;
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

  initSpongeValues();
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


void SEMsolver::initSpongeValues() {
  // Init all taper to 1 (default value)
  double alpha = -0.0001;
  LOOPHEAD(myMesh.getNumberOfNodes(), n)
    float x = myMesh.nodeCoord(n, 0);
    float y = myMesh.nodeCoord(n, 1);
    float z = myMesh.nodeCoord(n, 2);

    float distToFrontierX = min(myMesh.domainSize(0) - x, x);
    float distToFrontierY = min(myMesh.domainSize(1) - y, y);  // Fixed: use index 1
    float distToFrontierZ = min(myMesh.domainSize(2) - z, z);  // Fixed: use index 2

    // Find closest distance to domain's end in x y z coordinate
    float minDistToFrontier;
    if (!isSurface) minDistToFrontier = min(distToFrontierX, min(distToFrontierY, distToFrontierZ));
    else minDistToFrontier = min(distToFrontierY, distToFrontierZ);

    // Compute taper Coeff with proper Gaussian decay
    if (minDistToFrontier < m_spongeSize) {
        // Normalize distance to avoid underflow - alpha is already negative
        double normalizedDist = minDistToFrontier / m_spongeSize;
        spongeTaperCoeff(n) = std::exp(alpha * normalizedDist * normalizedDist);
    } else {
        spongeTaperCoeff(n) = 1.0;
    }
  LOOPEND

  FENCE
}
