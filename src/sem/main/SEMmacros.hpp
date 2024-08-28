#ifndef SEM_MACROS_HPP_
#define SEM_MACROS_HPP_

#include "commonConfig.hpp"

// define Macros test SEM proxy 2D case
#if defined (SEM2D)
#define DIMENSION  2
#define ROW 36
#define COL 4
#define ZEROED2D 0
// test SEM proxy 3D case
#else
#define DIMENSION  3
#define ROW 64
#define COL 6
#define ZEROED2D 1
#endif

// create views only for RAJA
#if defined (USE_RAJA)
#define CREATEVIEWS \
        arrayRealView rhsTerm=RHS_Term.toView(); \
        arrayRealView pnGlobal=PN_Global.toView(); \
        vectorIntView rhsElement=RHS_Element.toView(); \
        vectorRealView massMatrixGlobal=this->massMatrixGlobal.toView(); \
        vectorRealView yGlobal=this->yGlobal.toView(); \
        vectorRealView ShGlobal=this->ShGlobal.toView(); \
        arrayIntView globalNodesList=this->globalNodesList.toView(); \
        vectorRealView model=this->model.toView(); \
        vectorIntView listOfInteriorNodes=this->listOfInteriorNodes.toView(); \
        vectorIntView listOfBoundaryNodes=this->listOfBoundaryNodes.toView(); \
        arrayIntView faceInfos=this->faceInfos.toView(); \
        arrayIntView localFaceNodeToGlobalFaceNode=this->localFaceNodeToGlobalFaceNode.toView(); \
        vectorDoubleView weights=this->weights.toView(); \
        arrayRealView globalNodesCoords=this->globalNodesCoords.toView(); \
        arrayDoubleView derivativeBasisFunction1D=this->derivativeBasisFunction1D.toView(); \
        arrayIntView listOfElementsByColor=this->listOfElementsByColor.toView();
#else
#define CREATEVIEWS
#define RHS_Term rhsTerm
#define PN_Global pnGlobal
#define RHS_Element rhsElement
#endif

// define atomic add operation
#if defined(USE_RAJA) && !defined(SEM_MESHCOLOR)
  #define ATOMICADD(ADD1,ADD2) RAJA::atomicAdd< deviceAtomicPolicy >(&ADD1,ADD2);
#elif defined (USE_KOKKOS) && !defined(SEM_MESHCOLOR)
  #define ATOMICADD(ADD1,ADD2) Kokkos::atomic_add(&ADD1,ADD2)
#else
  #define ATOMICADD(ADD1,ADD2) ADD1+=ADD2
#endif

// define fence
#if defined (USE_RAJA)
  #define FENCE\
  if(timeSample%50==0)\
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, myInfo.numberOfNodes ), [pnGlobal] LVARRAY_HOST_DEVICE ( int i ){});
#elif defined (USE_KOKKOS)
  #define FENCE\
  if(timeSample%50==0) Kokkos::fence();
#else
  #define FENCE 
#endif


#endif //SEM_MACROS_HPP_
