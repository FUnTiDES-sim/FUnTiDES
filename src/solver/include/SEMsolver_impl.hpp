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
#include "fe/Integrals.hpp"
#ifdef USE_EZV
#include "ezvLauncher.hpp"
#include <cstdlib>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // USE_KOKKOS
#endif // USE_EZV

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::computeFEInit( Mesh const & mesh ) 
{
  m_mesh = mesh;
  allocateFEarrays( );
  initFEarrays( );
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
computeOneStep( const float &dt, 
                const int &timeSample,
                SolverBase::DataStruct & data )
{
  // Cast to the specific DataStruct type
  auto & myData = dynamic_cast<SEMsolverData &>(data);

  int const & i1 = myData.m_i1;
  int const & i2 = myData.m_i2;
  ARRAY_REAL_VIEW const & rhsTerm = myData.m_rhsTerm;
  ARRAY_REAL_VIEW const & pnGlobal = myData.m_pnGlobal;
  VECTOR_INT_VIEW const & rhsElement = myData.m_rhsElement;

  resetGlobalVectors();
  applyRHSTerm( dt, timeSample, i2, rhsTerm, rhsElement, pnGlobal);
  FENCE
  computeElementContributions( numPointsPerElem, i2, pnGlobal);
  FENCE
  updatePressureField( dt, rhsTerm, i1, i2, pnGlobal);
  FENCE
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
resetGlobalVectors() 
{
  
  LOOPHEAD( massMatrixGlobal.extent(0), i )
  {
    massMatrixGlobal[i] = 0;
    yGlobal[i] = 0;
  }
  LOOPEND
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
applyRHSTerm( float const dt,
              int timeSample, 
              int i2, 
              const ARRAY_REAL_VIEW &rhsTerm,
              const VECTOR_INT_VIEW &rhsElement, 
              const ARRAY_REAL_VIEW &pnGlobal) 
{
  float const dt2 = dt * dt;
  LOOPHEAD(rhsElement.extent(0), i)
  {
    int nodeRHS = globalNodesList(rhsElement[i], 0);
    float scale = dt2 * model[rhsElement[i]] * model[rhsElement[i]];
    pnGlobal(nodeRHS, i2) += scale * rhsTerm(i, timeSample);
  }
  LOOPEND
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
computeElementContributions( int nPointsPerElement,
                             int i2,
                             const ARRAY_REAL_VIEW &pnGlobal )
{
  int const numberOfElements = m_mesh.getNumberOfElements();
  MAINLOOPHEAD( numberOfElements, elementNumber)

  // Guard for extra threads (Kokkos might launch more than needed)
  if (elementNumber >= numberOfElements)
    return;

  float massMatrixLocal[ numPointsPerElem ] = {0};
  float pnLocal[ numPointsPerElem ] = {0};
  float Y[ numPointsPerElem ] = {0};

  for (int i = 0; i < nPointsPerElement; ++i) 
  {
    int const globalIdx = globalNodesList(elementNumber, i);
    pnLocal[i] = pnGlobal(globalIdx, i2);
  }


  INTEGRAL_TYPE::computeMassMatrixAndStiffnessVector( elementNumber, 
                                                      nPointsPerElement, 
                                                      globalNodesCoordsX, 
                                                      globalNodesCoordsY,
                                                      globalNodesCoordsZ, 
                                                      m_precomputedIntegralData,
                                                      massMatrixLocal, 
                                                      pnLocal, 
                                                      Y );

  auto const inv_model2 = 1.0f / (model[elementNumber] * model[elementNumber]);
  for (int i = 0; i < numPointsPerElem; ++i) 
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
updatePressureField( float const dt,
                     const ARRAY_REAL_VIEW &rhsTerm,
                     int i1, 
                     int i2, 
                     const ARRAY_REAL_VIEW &pnGlobal )
{
  float const dt2 = dt * dt;
  LOOPHEAD( spongeTaperCoeff.extent(0), I)
    pnGlobal(I, i1) = 2 * pnGlobal(I, i2) - pnGlobal(I, i1) - dt2 * yGlobal[I] / massMatrixGlobal[I];
    pnGlobal(I, i1) *= spongeTaperCoeff(I);
    pnGlobal(I, i2) *= spongeTaperCoeff(I);
  LOOPEND
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
outputPnValues( const int &indexTimeStep, int &i1,
                               int &myElementSource,
                               const ARRAY_REAL_VIEW &pnGlobal) {
  // writes debugging ascii file.
  if (indexTimeStep % 50== 0) {
    cout << "TimeStep=" << indexTimeStep
         << ";  pnGlobal @ elementSource location " << myElementSource
         << " after computeOneStep = "
         << pnGlobal(globalNodesList(myElementSource, 0), i1) << endl;
#ifdef SEM_SAVE_SNAPSHOTS
    m_mesh.saveSnapShot(indexTimeStep, i1, pnGlobal);
#endif // SEM_SAVE_SNAPSHOTS
  }
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
initFEarrays() 
{
  int const numberOfElements = m_mesh.getNumberOfElements();
  int const numberOfInteriorNodes = m_mesh.getNumberOfInteriorNodes();
  // interior elements
  m_mesh.globalNodesList( numberOfElements, globalNodesList);
  m_mesh.getListOfInteriorNodes(numberOfInteriorNodes,
                              listOfInteriorNodes);
  // mesh coordinates
  m_mesh.nodesCoordinates(globalNodesCoordsX, globalNodesCoordsZ,
                        globalNodesCoordsY);

  // get model
  m_mesh.getModel( numberOfElements, model);

  // get minimal wavespeed
  double min;
  auto model_ = this->model; // Avoid implicit capture
#ifdef USE_KOKKOS
  Kokkos::parallel_reduce( "vMinFind", 
                           numberOfElements,
                           KOKKOS_LAMBDA(const int &e, double &lmin) 
                           {
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

  INTEGRAL_TYPE::init( m_precomputedIntegralData);

  // Sponge boundaries
  initSpongeValues(m_mesh);
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
allocateFEarrays() 
{
  int const numberOfElements = m_mesh.getNumberOfElements();
  int const numberOfNodes = m_mesh.getNumberOfNodes();
  int nbQuadraturePoints = numPointsPerElem;
  int const numberOfInteriorNodes = m_mesh.getNumberOfInteriorNodes();
  int const numberOfDampingNodes = m_mesh.getNumberOfDampingNodes();
  // interior elements
  cout << "Allocate host memory for arrays in the solver ..." << endl;
  globalNodesList = allocateArray2D<ARRAY_INT_VIEW>(numberOfElements,
                                                    numPointsPerElem,
                                              "globalNodesList");
  listOfInteriorNodes = allocateVector<VECTOR_INT_VIEW>(numberOfInteriorNodes,
                                                  "listOfInteriorNodes");
  listOfDampingNodes = allocateVector<VECTOR_INT_VIEW>(numberOfDampingNodes,
                                                 "listOfDampingNodes");

  // global coordinates
  globalNodesCoordsX = allocateArray2D<ARRAY_REAL_VIEW>( numberOfElements, nbQuadraturePoints, "globalNodesCoordsX");
  globalNodesCoordsY = allocateArray2D<ARRAY_REAL_VIEW>( numberOfElements, nbQuadraturePoints, "globalNodesCoordsY");
  globalNodesCoordsZ = allocateArray2D<ARRAY_REAL_VIEW>( numberOfElements, nbQuadraturePoints, "globalNodesCoordsZ");

  model = allocateVector<VECTOR_REAL_VIEW>(numberOfElements, "model");

  cout << "Allocate model ..." << endl;
  
  // shared arrays
  massMatrixGlobal =
      allocateVector<VECTOR_REAL_VIEW>( numberOfNodes, "massMatrixGlobal");
  yGlobal = allocateVector<VECTOR_REAL_VIEW>( numberOfNodes, "yGlobal");

  // sponge allocation
  spongeTaperCoeff =
      allocateVector<VECTOR_REAL_VIEW>( numberOfNodes, "spongeTaperCoeff");
}

template< int ORDER, typename INTEGRAL_TYPE >
void 
SEMsolver<ORDER, INTEGRAL_TYPE>::
initSpongeValues(Mesh &mesh) 
{
  int const numberOfNodes = mesh.getNumberOfNodes();
  // Init all taper to 1 (default value)
  LOOPHEAD( numberOfNodes, i)
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
  LOOPHEAD(spongeTaperCoeff.extent(0), i)
  pnGlobal(i, i1) *= spongeTaperCoeff(i);
  pnGlobal(i, i2) *= spongeTaperCoeff(i);
  LOOPEND
}


extern template class SEMsolver< 2, IntegralTypeSelector< 2, 0 >::type >;
extern template class SEMsolver< 2, IntegralTypeSelector< 2, 1 >::type >;
extern template class SEMsolver< 2, IntegralTypeSelector< 2, 2 >::type >;
extern template class SEMsolver< 2, IntegralTypeSelector< 2, 3 >::type >;
