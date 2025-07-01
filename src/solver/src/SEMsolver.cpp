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

void SEMsolver::allocateSolverDIVA(SEMinfo &myInfo_in) {
  myInfo = &myInfo_in;
  order = myInfo_in.myOrderNumber;
  allocateFEarraysWithoutMesh(myInfo_in);
}

void SEMsolver::FEInitDIVA(SEMinfo &myInfo,
                            const array4DInt &elemsToNodesDIVA,
                            const arrayReal &nodeCoordsDIVA,
                            const vectorInt &interiorNodes,
                            const vectorReal &rhomodelOnNodes,
                            const vectorReal &vpmodelOnNodes) {
  initFEarraysDIVA(myInfo, elemsToNodesDIVA, nodeCoordsDIVA, interiorNodes, rhomodelOnNodes, vpmodelOnNodes);
}

void SEMsolver::FEInitDIVA_(SEMinfo &myInfo,
                              Kokkos::Experimental::python_view_type_t<Kokkos::View<int64_t ****, Layout, MemSpace>> elemsToNodesDIVA,
                              Kokkos::Experimental::python_view_type_t<Kokkos::View<float **, Layout, MemSpace>> nodeCoordsDIVA,
                              Kokkos::Experimental::python_view_type_t<Kokkos::View<int *, Layout, MemSpace>> interiorNodes,
                              Kokkos::Experimental::python_view_type_t<Kokkos::View<float *, Layout, MemSpace>> rhomodelOnNodes,
                              Kokkos::Experimental::python_view_type_t<Kokkos::View<float *, Layout, MemSpace>> vpmodelOnNodes) {
   FEInitDIVA(myInfo,
              elemsToNodesDIVA,
              nodeCoordsDIVA,
              interiorNodes,
              rhomodelOnNodes,
              vpmodelOnNodes);
}

void SEMsolver::allocateFEarraysWithoutMesh(SEMinfo &myInfo) {
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

  // TODO quick and dirty fix to avoid uninitialized values in spongeTaperCoeff
  LOOPHEAD(myInfo.numberOfNodes, i)
  spongeTaperCoeff(i) = 1;
  LOOPEND
}

void SEMsolver::initFEarraysDIVA(SEMinfo &myInfo,
                                  const array4DInt &elemsToNodesDIVA,
                                  const arrayReal &nodeCoordsDIVA,
                                  const vectorInt &interiorNodes,
                                  const vectorReal &rhomodelOnNodes,
                                  const vectorReal &vpmodelOnNodes) {

  globalelemsToNodesDIVA = elemsToNodesDIVA;
  globalnodeCoordsDIVA   = nodeCoordsDIVA;
  listOfInteriorNodes    = interiorNodes;
  vpmodel                = vpmodelOnNodes;
  rhomodel               = rhomodelOnNodes;

  // get minimal wavespeed
  double min;
  auto model_ = this->vpmodel; // Avoid implicit capture
#ifdef USE_KOKKOS
  Kokkos::parallel_reduce(
        "vMinFind", myInfo.numberOfNodes,
        KOKKOS_LAMBDA(const int &n, double &lmin) {
        double val = model_[n];
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
  //initSpongeValues(mesh, myInfo); TODO not yet available from mesher
  FENCE
}

// void SEMsolver::initFEarraysDIVA_(SEMinfo &myInfo,
//                                          Kokkos::Experimental::python_view_type_t<Kokkos::View<int ****, Layout, MemSpace>> elemsToNodesDIVA,
//                                          Kokkos::Experimental::python_view_type_t<Kokkos::View<float **, Layout, MemSpace>> nodeCoordsDIVA,
//                                          Kokkos::Experimental::python_view_type_t<Kokkos::View<int *, Layout, MemSpace>> interiorNodes,
//                                          Kokkos::Experimental::python_view_type_t<Kokkos::View<float *, Layout, MemSpace>> rhomodelOnNodes,
//                                          Kokkos::Experimental::python_view_type_t<Kokkos::View<float *, Layout, MemSpace>> vpmodelOnNodes) {
//    initFEarraysDIVA(myInfo,
//                       elemsToNodesDIVA,
//                       nodeCoordsDIVA,
//                       interiorNodes,
//                       rhomodelOnNodes,
//                       vpmodelOnNodes);
//    }

void SEMsolver::computeOneStepDIVA(const int &timeSample, const int &order,
                               const int &nPointsPerElement, const int &i1,
                               const int &i2, SEMinfo &myInfo,
                               const arrayReal &rhsTerm,
                               const arrayReal &pnGlobal,
                               const vectorInt &rhsElement) {
  resetGlobalVectors(myInfo.numberOfNodes);
  FENCE
  computeElementContributionsDIVA(order, nPointsPerElement, myInfo, i2, pnGlobal);
  FENCE
  applyRHSTerm(timeSample, rhsTerm, rhsElement, myInfo);
  FENCE
  updatePressureField(i1, i2, myInfo, pnGlobal);
  FENCE
  spongeUpdate(pnGlobal, i1, i2);
  FENCE
}

#ifdef ENABLE_PYWRAP
void SEMsolver::computeOneStepDIVA_wrapper(
    int t, int order, int npts, int i1, int i2, SEMinfo info,
    Kokkos::Experimental::python_view_type_t<
        Kokkos::View<float **, Layout, MemSpace>>
        rhsTerm,
    Kokkos::Experimental::python_view_type_t<
        Kokkos::View<float **, Layout, MemSpace>>
        pnGlobal,
    Kokkos::Experimental::python_view_type_t<
        Kokkos::View<int *, Layout, MemSpace>>
        rhsElement) {
  // arrayReal pnGlobal_raw(pnGlobal);
  computeOneStepDIVA(t, order, npts, i1, i2, info, rhsTerm, pnGlobal, rhsElement);
}
#endif // ENABLE_PYWRAP

void SEMsolver::resetGlobalVectors(int numNodes) {
  LOOPHEAD(numNodes, i)
  massMatrixGlobal[i] = 0;
  yGlobal[i] = 0;
  LOOPEND
}

void SEMsolver::applyRHSTerm(int timeSample, const arrayReal &rhsTerm,
                             const vectorInt &rhsElement, SEMinfo &myInfo) {
  LOOPHEAD(myInfo.myNumberOfRHS, i)
  int nodeRHS = globalelemsToNodesDIVA(0,0,0,rhsElement[i]); // TODO modify when we have the coeffs
  yGlobal[nodeRHS] -= rhsTerm(i, timeSample);
  LOOPEND
}


void SEMsolver::computeElementContributionsDIVA(int order, int nPointsPerElement,
                                            SEMinfo &myInfo, int i2,
                                            const arrayReal &pnGlobal) {
  MAINLOOPHEAD(myInfo.numberOfElements, elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= myInfo.numberOfElements)
    return;

  float massMatrixLocal[ROW] = {0};
  float pnLocal[ROW] = {0};
  float Y[ROW] = {0};

  int qindex = 0;
  for (int k = 0; k < order + 1; ++k) {
    for (int j = 0; j < order + 1; ++j) {
      for (int i = 0; i < order + 1; ++i) {
           int globalIdx = globalelemsToNodesDIVA(i, j, k, elementNumber);
           pnLocal[qindex] = pnGlobal(globalIdx, i2);
           qindex++;
      }
    }
  }

#if defined(USE_SEMOPTIM) || defined(USE_SHIVA)
  constexpr int ORDER = SEMinfo::myOrderNumber;
  myQkIntegrals.computeMassMatrixAndStiffnessVectorDIVA<ORDER>(
      elementNumber, nPointsPerElement, globalelemsToNodesDIVA, globalnodeCoordsDIVA,
      massMatrixLocal, pnLocal, Y, rhomodel, vpmodel);
#endif

  qindex = 0;  
  for (int k = 0; k < order + 1; ++k) {
    for (int j = 0; j < order + 1; ++j) {
      for (int i = 0; i < order + 1; ++i) {
           int gIndex = globalelemsToNodesDIVA(i, j, k, elementNumber);
           float massValue = massMatrixLocal[qindex];
           ATOMICADD(massMatrixGlobal[gIndex], massValue);
           ATOMICADD(yGlobal[gIndex], Y[qindex]);
           qindex++;
       }
    }
  }

  MAINLOOPEND
}

void SEMsolver::updatePressureField(int i1, int i2, SEMinfo &myInfo,
                                    const arrayReal &pnGlobal) {
  LOOPHEAD(myInfo.numberOfInteriorNodes, i)
  int I = listOfInteriorNodes[i];
  // if (I >7550000 & I < 7570000) {
  // // Debugging output  
  // printf("RHSTerm[%d] = %f\n", I, yGlobal[I]);
  // printf("pnGlobal[%d][%d] = %f\n", I, i2, pnGlobal(I, i2));
  // printf("pnGlobal[%d][%d] = %f\n", I, i1, pnGlobal(I, i1));
  // printf("massMatrixGlobal[%d] = %f\n", I, massMatrixGlobal[I]);
  // printf("rho[%d] = %f\n", I, rhomodel[I]);
  // printf("vp[%d] = %f\n", I, vpmodel[I]);
  // }
  // Update pressure field
  pnGlobal(I, i1) =
      2 * pnGlobal(I, i2) - pnGlobal(I, i1) -
      myInfo.myTimeStep * myInfo.myTimeStep * yGlobal[I] / massMatrixGlobal[I];
  LOOPEND
}

// void SEMsolver::initSpongeValues(Mesh &mesh, SEMinfo &myInfo) {
//   // Init all taper to 1 (default value)
//   LOOPHEAD(myInfo.numberOfNodes, i)
//   spongeTaperCoeff(i) = 1;
//   LOOPEND

//   // int n = 0;
//   double alpha = -0.0001;
//   int spongeSize = mesh.getSpongeSize();
//   int nx = mesh.getNx();
//   int ny = mesh.getNy();
//   int nz = mesh.getNz();

//   // Update X boundaries
//   LOOPHEAD(nz, k)
//   for (int j = 0; j < ny; j++) {
//     // lower x
//     for (int i = 0; i <= spongeSize; i++) {
//       int n = mesh.ijktoI(i, j, k);
//       double value = spongeSize - i;
//       spongeTaperCoeff(n) = std::exp(alpha * static_cast<float>(value * value));
//     }
//     // upper x
//     for (int i = nx - spongeSize - 1; i < nx; i++) {
//       int n = mesh.ijktoI(i, j, k);
//       double value = spongeSize - (nx - i);
//       spongeTaperCoeff(n) = std::exp(alpha * static_cast<float>(value * value));
//     }
//   }
//   LOOPEND

//   // Update Y boundaries
//   // for (int k = 0; k < nz; k++) {
//   LOOPHEAD(nz, k)
//   int n;
//   for (int i = 0; i < nx; i++) {
//     // lower y
//     for (int j = 0; j <= spongeSize; j++) {
//       n = mesh.ijktoI(i, j, k);
//       double value = spongeSize - j;
//       spongeTaperCoeff(n) =
//           std::exp(alpha * static_cast<double>(value * value));
//     }
//     // upper y
//     for (int j = ny - spongeSize - 1; j < ny; j++) {
//       n = mesh.ijktoI(i, j, k);
//       double value = spongeSize - (ny - j);
//       spongeTaperCoeff(n) =
//           std::exp(alpha * static_cast<double>(value * value));
//     }
//   }
//   LOOPEND

//   // Update Z boundaries
//   LOOPHEAD(ny, j)
//   int n;
//   for (int i = 0; i < nx; i++) {
//     // lower z
//     for (int k = 0; k <= spongeSize; k++) {
//       n = mesh.ijktoI(i, j, k);
//       double value = spongeSize - k;
//       spongeTaperCoeff(n) =
//           std::exp(alpha * static_cast<double>(value * value));
//     }
//     // upper z
//     for (int k = nz - spongeSize - 1; k < nz; k++) {
//       n = mesh.ijktoI(i, j, k);
//       double value = spongeSize - (nz - k);
//       spongeTaperCoeff(n) =
//           std::exp(alpha * static_cast<double>(value * value));
//     }
//   }
//   LOOPEND

//   FENCE
// }

void SEMsolver::spongeUpdate(const arrayReal &pnGlobal, const int i1,
                             const int i2) {
  // for (int i = 0; i < myInfo->numberOfNodes; i++) {
  LOOPHEAD(myInfo->numberOfNodes, i)
  pnGlobal(i, i1) *= spongeTaperCoeff(i);
  pnGlobal(i, i2) *= spongeTaperCoeff(i);
  LOOPEND
}
// }
