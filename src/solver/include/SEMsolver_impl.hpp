//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "finiteElement/classic/SEMQkGLIntegralsClassic.hpp"
#include "finiteElement/optim/SEMQkGLIntegralsOptim.hpp"
#include "finiteElement/geos/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "finiteElement/shiva/SEMQkGLIntegralsShiva.hpp"

#include "SEMsolver.hpp"
#include "dataType.hpp"
#ifdef USE_EZV
#include "ezvLauncher.hpp"
#include <cstdlib>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // USE_KOKKOS
#endif // USE_EZV

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::computeFEInit( SEMinfo const & myInfo_in, 
                                                Mesh const & mesh ) 
{
  this->myInfo = &myInfo_in;
  myMesh = mesh;
  order = myInfo_in.myOrderNumber;
  allocateFEarrays(myInfo_in);
  initFEarrays(myInfo_in, mesh);
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
computeOneStep( const int &timeSample, 
                const int &nPointsPerElement, 
                const int &i1,
                const int &i2, 
                SEMinfo const & myInfo,
                const ARRAY_REAL_VIEW &rhsTerm, 
                const ARRAY_REAL_VIEW &pnGlobal,
                const VECTOR_INT_VIEW &rhsElement )
{
  resetGlobalVectors(myInfo.numberOfNodes);
  applyRHSTerm(timeSample, i2, rhsTerm, rhsElement, myInfo, pnGlobal);
  FENCE
  computeElementContributions(order, nPointsPerElement, myInfo, i2, pnGlobal);
  FENCE
  updatePressureField(i1, i2, myInfo, pnGlobal);
  FENCE
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
resetGlobalVectors(int numNodes) {
  LOOPHEAD(numNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
applyRHSTerm(int timeSample, int i2, const ARRAY_REAL_VIEW &rhsTerm,
                             const VECTOR_INT_VIEW &rhsElement, SEMinfo const &myInfo,
                             const ARRAY_REAL_VIEW &pnGlobal) 
{
  float const dt2 = myInfo.myTimeStep * myInfo.myTimeStep;
  LOOPHEAD(myInfo.myNumberOfRHS, i)
    int nodeRHS = globalNodesList(rhsElement[i], 0);
    float scale = dt2 * model[rhsElement[i]] * model[rhsElement[i]];
    pnGlobal(nodeRHS, i2) += scale * rhsTerm(i, timeSample);
  LOOPEND
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
computeElementContributions(int order, int nPointsPerElement,
                                            SEMinfo const &myInfo, int i2,
                                            const ARRAY_REAL_VIEW &pnGlobal) {
  MAINLOOPHEAD(myInfo.numberOfElements, elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= myInfo.numberOfElements)
    return;

  float massMatrixLocal[SEMinfo::nPointsPerElement] = {0};
  float pnLocal[SEMinfo::nPointsPerElement] = {0};
  float Y[SEMinfo::nPointsPerElement] = {0};

  for (int i = 0; i < nPointsPerElement; ++i) 
  {
    int const globalIdx = globalNodesList(elementNumber, i);
    pnLocal[i] = pnGlobal(globalIdx, i2);
  }

#ifdef USE_SEMCLASSIC
  myQkIntegrals.computeMassMatrixAndStiffnessVector(elementNumber, 
                                                    order,
                                                    nPointsPerElement, 
                                                    globalNodesCoordsX,
                                                    globalNodesCoordsY, 
                                                    globalNodesCoordsZ, 
                                                    weights,
                                                    derivativeBasisFunction1D, 
                                                    massMatrixLocal, 
                                                    pnLocal, 
                                                    Y);
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

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
updatePressureField(int i1, int i2, SEMinfo const &myInfo,
                    const ARRAY_REAL_VIEW &pnGlobal) 
{

  float const dt2 = myInfo.myTimeStep * myInfo.myTimeStep;
  LOOPHEAD(myInfo.numberOfNodes, I)
    pnGlobal(I, i1) = 2 * pnGlobal(I, i2) - pnGlobal(I, i1) - dt2 * yGlobal[I] / massMatrixGlobal[I];
    pnGlobal(I, i1) *= spongeTaperCoeff(I);
    pnGlobal(I, i2) *= spongeTaperCoeff(I);
  LOOPEND
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
outputPnValues(Mesh mesh, const int &indexTimeStep, int &i1,
                               int &myElementSource,
                               const ARRAY_REAL_VIEW &pnGlobal) {
  // writes debugging ascii file.
  if (indexTimeStep % 50== 0) {
    cout << "TimeStep=" << indexTimeStep
         << ";  pnGlobal @ elementSource location " << myElementSource
         << " after computeOneStep = "
         << pnGlobal(globalNodesList(myElementSource, 0), i1) << endl;
#ifdef SEM_SAVE_SNAPSHOTS
    mesh.saveSnapShot(indexTimeStep, i1, pnGlobal);
#endif // SEM_SAVE_SNAPSHOTS
  }
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
initFEarrays( SEMinfo const & myInfo, Mesh mesh) {
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
  myQkBasis.gaussLobattoQuadraturePoints(order,quadraturePoints);
  // get gauss-lobatto weights
  myQkBasis.gaussLobattoQuadratureWeights(order,weights);
  // get basis function and corresponding derivatives
  myQkBasis.getDerivativeBasisFunction1D(order,quadraturePoints,
                                         derivativeBasisFunction1D);
#endif // USE_SEMCLASSIC

  // Sponge boundaries
  initSpongeValues(mesh, myInfo);
  FENCE;
}

//************************************************************************
//  Allocate arrays for the solver
//  This function allocates all arrays needed for the solver
//  It allocates arrays for global nodes, global coordinates, and sponge
//  It also allocates arrays for the mass matrix and the global pressure field
//************************************************************************
template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
allocateFEarrays( SEMinfo const & myInfo) 
{
  int nbQuadraturePoints = (order + 1) * (order + 1) * (order + 1);
  // interior elements
  cout << "Allocate host memory for arrays in the solver ..." << endl;
  globalNodesList = allocateArray2D<ARRAY_INT_VIEW>(myInfo.numberOfElements,
                                              myInfo.numberOfPointsPerElement,
                                              "globalNodesList");
  listOfInteriorNodes = allocateVector<VECTOR_INT_VIEW>(myInfo.numberOfInteriorNodes,
                                                  "listOfInteriorNodes");
  listOfDampingNodes = allocateVector<VECTOR_INT_VIEW>(myInfo.numberOfDampingNodes,
                                                 "listOfDampingNodes");

  // global coordinates
  globalNodesCoordsX = allocateArray2D<ARRAY_REAL_VIEW>(
      myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsX");
  globalNodesCoordsY = allocateArray2D<ARRAY_REAL_VIEW>(
      myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsY");
  globalNodesCoordsZ = allocateArray2D<ARRAY_REAL_VIEW>(
      myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsZ");

  model = allocateVector<VECTOR_REAL_VIEW>(myInfo.numberOfElements, "model");

  cout << "Allocate model ..." << endl;
  
  #ifdef USE_SEMCLASSIC
  quadraturePoints = allocateVector<VECTOR_REAL_VIEW>(nbQuadraturePoints, "quadraturePoints");
  weights = allocateVector<VECTOR_REAL_VIEW>(nbQuadraturePoints, "weights");
  derivativeBasisFunction1D = allocateArray2D<ARRAY_REAL_VIEW>(myInfo.numberOfNodes,
                                                           nbQuadraturePoints, "derivativeBasisFunction1D");
  #endif // USE_SEMCLASSIC

  // shared arrays
  massMatrixGlobal =
      allocateVector<VECTOR_REAL_VIEW>(myInfo.numberOfNodes, "massMatrixGlobal");
  yGlobal = allocateVector<VECTOR_REAL_VIEW>(myInfo.numberOfNodes, "yGlobal");

  // sponge allocation
  spongeTaperCoeff =
      allocateVector<VECTOR_REAL_VIEW>(myInfo.numberOfNodes, "spongeTaperCoeff");
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
initSpongeValues(Mesh &mesh, SEMinfo const &myInfo) 
{
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

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
spongeUpdate( const ARRAY_REAL_VIEW &pnGlobal, 
              const int i1,
              const int i2) 
{
  // for (int i = 0; i < myInfo->numberOfNodes; i++) {
  LOOPHEAD(myInfo->numberOfNodes, i)
  pnGlobal(i, i1) *= spongeTaperCoeff(i);
  pnGlobal(i, i2) *= spongeTaperCoeff(i);
  LOOPEND
}


template class SEMsolver< 2, SEMQkGLIntegralsOptim< 2, float, float > >;
template class SEMsolver< 3, SEMQkGLIntegralsOptim< 3, float, float > >;
template class SEMsolver< 4, SEMQkGLIntegralsOptim< 4, float, float > >;
