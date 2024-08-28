#include <Kokkos_Core.hpp>
#include <cstdio>

int main( int argc, char * argv[] )
{
  Kokkos::initialize( argc, argv );

  int N = atoi( argv[1] );

  Kokkos::parallel_for( "Loop1", N, KOKKOS_LAMBDA ( const int i )
  {
    printf( "Greeting from iteration %i\n", i );
  } );

  Kokkos::finalize();
}
