#include "commonConfig.hpp"

// define Macros for function type
#if defined (USE_RAJA)
  #define PROXY_HOST_DEVICE LVARRAY_HOST_DEVICE LVARRAY_FORCE_INLINE
#elif defined (USE_KOKKOS)
  #define PROXY_HOST_DEVICE KOKKOS_INLINE_FUNCTION 
#else
  #define PROXY_HOST_DEVICE 
#endif

#if defined (USE_KOKKOS)
  #define LOOPHEAD(Range, Iterator)\
    Kokkos::parallel_for( Range, KOKKOS_CLASS_LAMBDA ( const int Iterator ){
  #define LOOPEND   });

#elif defined (USE_RAJA)
  #define LOOPHEAD(Range, Iterator)\
    RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, Range),  [=] LVARRAY_HOST_DEVICE  ( int Iterator ) {
  #define LOOPEND   });

#elif defined (USE_OMP)
  #define LOOPHEAD(Range, Iterator)\
    _Pragma("omp parallel for")\
    for( int Iterator=0; Iterator<Range; Iterator++ ){
  #define LOOPEND   }

#else 
  // the sequential case
  #define LOOPHEAD(Range, Iterator)\
    for( int Iterator=0; Iterator<Range; Iterator++ ){
  #define LOOPEND   }
#endif


#if defined (USE_KOKKOS) && defined (USE_KOKKOS_TEAMS)
  #define LaunchMaxThreadsPerBlock 64 
  #define LaunchMinBlocksPerSM 1
  #define nthreads 64
  #define MAINLOOPHEAD(Range, Iterator)\
    const int leagueSize=(Range-1)/nthreads; \
    const Kokkos::TeamPolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>> teamPolicy(leagueSize, nthreads); \
    Kokkos::parallel_for("Loop", teamPolicy, KOKKOS_CLASS_LAMBDA ( const Kokkos::TeamPolicy<>::member_type & thread ) { \
      Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, nthreads), [=] (const int index) { \
        int Iterator = thread.league_rank()*nthreads+index;
  #define MAINLOOPEND }); });

#elif defined (USE_KOKKOS) && !defined(SEM_MESHCOLOR)
  #define LaunchMaxThreadsPerBlock 64
  #define LaunchMinBlocksPerSM 1
  #define MAINLOOPHEAD(Range, Iterator)\
    Kokkos::parallel_for( Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, Range), \
                          KOKKOS_CLASS_LAMBDA ( const int Iterator ){
  #define MAINLOOPEND   });

#elif defined (SEM_MESHCOLOR)
  #define MAINLOOPHEAD(Range, Iterator)\
    for (int color=0; color<myInfo.numberOfColors;color++) {\
    LOOPHEAD( myInfo.numberOfElementsByColor[color], eColor) \
    int Iterator=listOfElementsByColor(color,eColor);
  #define MAINLOOPEND   }); }

#else
  #define MAINLOOPHEAD LOOPHEAD
  #define MAINLOOPEND LOOPEND
#endif


#if defined (USE_RAJA)
  #define ARRAY_DOUBLE_VIEW arrayDoubleView 
  #define ARRAY_REAL_VIEW arrayRealView
  #define ARRAY_INT_VIEW arrayIntView
  #define VECTOR_DOUBLE_VIEW vectorDoubleView
  #define VECTOR_REAL_VIEW vectorRealView
  #define VECTOR_INT_VIEW vectorIntView
#else
  #define ARRAY_DOUBLE_VIEW arrayDouble
  #define ARRAY_REAL_VIEW arrayReal
  #define ARRAY_INT_VIEW arrayInt
  #define VECTOR_DOUBLE_VIEW vectorDouble
  #define VECTOR_REAL_VIEW vectorReal
  #define VECTOR_INT_VIEW vectorInt
#endif

#if defined (USE_KOKKOS)
  #define KOKKOSNAME "v",
#else
  #define KOKKOSNAME 
#endif
