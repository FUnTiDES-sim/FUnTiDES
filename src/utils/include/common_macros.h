#include "common_config.h"

// define Macros for function type
#if defined(USE_KOKKOS)
#define PROXY_HOST_DEVICE KOKKOS_FORCEINLINE_FUNCTION
#else
#define PROXY_HOST_DEVICE
#endif

// LOOPHEAD
#if defined(USE_KOKKOS)
#define LOOPHEAD(Range, Iterator) \
    Kokkos::parallel_for( Range, KOKKOS_CLASS_LAMBDA ( const int Iterator ){
#define LOOPEND \
  });

#else
// the sequential case
#define LOOPHEAD(Range, Iterator)                      \
  for (int Iterator = 0; Iterator < Range; Iterator++) \
  {
#define LOOPEND }
#endif

// MAINLOOP
#if defined(USE_KOKKOS) && defined(USE_KOKKOS_TEAMS)
#define MAINLOOPHEAD(Range, Iterator)                                \
  const int nthreads = 128;                                          \
  const int leagueSize = (Range - 1) / nthreads;                     \
  const Kokkos::TeamPolicy<> teamPolicy(leagueSize, Kokkos::AUTO()); \
    Kokkos::parallel_for("teamLoop", teamPolicy,KOKKOS_CLASS_LAMBDA(const Kokkos::TeamPolicy<>::member_type & team_member ){\
		int Iterator=team_member.league_rank () * team_member.team_size () + team_member.team_rank ();
#define MAINLOOPEND \
  });

#elif defined(USE_KOKKOS) && !defined(SEM_MESHCOLOR)
#define LaunchMaxThreadsPerBlock 64
#define LaunchMinBlocksPerSM 1
#define MAINLOOPHEAD(Range, Iterator)                                                                                          \
    Kokkos::parallel_for( Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, Range), \
                          KOKKOS_CLASS_LAMBDA ( const int Iterator ){
#define MAINLOOPEND \
  });

#elif defined(SEM_MESHCOLOR)
#define MAINLOOPHEAD(Range, Iterator)                         \
  for (int color = 0; color < myInfo.numberOfColors; color++) \
  {                                                           \
    LOOPHEAD(myInfo.numberOfElementsByColor[color], eColor)   \
    int Iterator = listOfElementsByColor(color, eColor);
#define MAINLOOPEND \
  });               \
  }

#else
#define MAINLOOPHEAD LOOPHEAD
#define MAINLOOPEND LOOPEND
#endif

// FIND_MAX FOR 1D ARRAY
// The kokkos case
#if defined(USE_KOKKOS)
#define FIND_MAX_1D(Array, Range, Result)                                  \
  if (Array.extent(0) == 0)                                                \
    throw std::runtime_error("Error in FIND_MAX_1D: Array has zero size"); \
  Kokkos::parallel_reduce(                                                 \
      "FindMax1D", Range,                                                  \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_max) {      \
        if (Array[i] > local_max) local_max = Array[i];                    \
      },                                                                   \
      Kokkos::Max<decltype(Result)>(Result));
#else
// The sequential case
#define FIND_MAX_1D(Array, Range, Result)                                  \
  if (Array.extent(0) == 0)                                                \
    throw std::runtime_error("Error in FIND_MAX_1D: Array has zero size"); \
  Result = Array[0];                                                       \
  for (int i = 1; i < Range; i++)                                          \
  {                                                                        \
    if (Array[i] > Result) Result = Array[i];                              \
  }
#endif

// FIND_MIN
#if defined(USE_KOKKOS)
#define FIND_MIN(Array, Range, Result)                                \
  Kokkos::parallel_reduce(                                            \
      Range,                                                          \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_min) { \
        local_min = Array[i];                                         \
      },                                                              \
      Kokkos::Min<decltype(Result)>(Result));
#else
// The sequential case
#define FIND_MIN(Array, Range, Result)        \
  Result = Array[0];                          \
  for (int i = 1; i < Range; i++)             \
  {                                           \
    if (Array[i] < Result) Result = Array[i]; \
  }
#endif

// SUM
#if defined(USE_KOKKOS)
#define SUM(Array, Range, Result)                                     \
  Kokkos::parallel_reduce(                                            \
      Range,                                                          \
      KOKKOS_CLASS_LAMBDA(const int i, decltype(Result)& local_sum) { \
        local_sum = Array[i];                                         \
      },                                                              \
      Kokkos::Sum<decltype(Result)>(Result));
#else
// The sequential case
#define SUM(Array, Range, Result) \
  Result = decltype(Result){0};   \
  for (int i = 0; i < Range; i++) \
  {                               \
    Result += Array[i];           \
  }
#endif

#define ARRAY_DOUBLE_VIEW arrayReal
#define ARRAY_REAL_VIEW arrayReal
#define ARRAY_INT_VIEW arrayInt
#define VECTOR_DOUBLE_VIEW vectorReal
#define VECTOR_REAL_VIEW vectorReal
#define VECTOR_INT_VIEW vectorInt

#define ARRAY_TYPE_VIEW arrayReal
#define VECTOR_TYPE_VIEW vectorReal

#if defined(USE_KOKKOS)
#define KOKKOSNAME "v",
#else
#define KOKKOSNAME
#endif
