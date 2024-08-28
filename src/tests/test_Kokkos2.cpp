#include <Kokkos_Core.hpp>
#include <cstdio>

int main( int argc, char * argv[] )
{
  struct Foo {
    KOKKOS_INLINE_FUNCTION
    // Each team handles a slice of the data
    void operator() (const Kokkos::TeamPolicy<>::member_type& thread) const {
        printf("Loop 2: Greetings from thread %i of team %i out of (thread_size %d, team_size %d)\n",
                thread.team_rank(),thread.league_rank(), thread.team_size(), thread.league_size() );
    }
};

  Kokkos::initialize( argc, argv );

  int N = (argc > 1)? atoi( argv[1] ) : 10;

  Kokkos::parallel_for( "Loop1", N, KOKKOS_LAMBDA ( const int i )
  {
    printf( "Loop 1: Greeting from iteration %i\n", i );
  } );


  // Set up TeamPolicy with N teams with maximum number of threads per team
  // and 16 vector lanes. Kokkos::AUTO will determine the number of threads
  const Kokkos::TeamPolicy<> policy1(N, Kokkos::AUTO, 16);
  // The maximum vector length is hardware dependent but can always be smaller
  // than the hardware allows. The vector length must be a power of 2.
  Foo foo;
  Kokkos::parallel_for("Loop2", policy1, foo);


  const Kokkos::TeamPolicy<> policy2(N, 12);
  Kokkos::parallel_for("Loop3", policy2, KOKKOS_LAMBDA ( const Kokkos::TeamPolicy<>::member_type & thread ) {
        printf("Loop3: Greetings from thread %i of team %i out of (thread_size %d, team_size %d) \n",
                thread.team_rank(),thread.league_rank(), thread.team_size(), thread.league_size() );
  });

  const Kokkos::TeamPolicy<> policy3(N, 3);
  Kokkos::parallel_for("Loop4outer", policy3, KOKKOS_LAMBDA ( const Kokkos::TeamPolicy<>::member_type & thread ) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, 31), KOKKOS_LAMBDA (const int i) {
        printf("Loop3: Greetings from thread %i of team %i out of (thread_size %d, team_size %d) executed loop %d \n",
                thread.team_rank(),thread.league_rank(), thread.team_size(), thread.league_size(), i);
    });
  });
 
  Kokkos::finalize();
}
