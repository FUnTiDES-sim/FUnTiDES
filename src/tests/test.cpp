/* Copyright (c) 2019 CSC Training */
/* Copyright (c) 2021 ENCCS */
#include <stdio.h>
#include <math.h>
//#define NX 1024000

int main( void )
{
  constexpr int NX=102400;
  double vecA[NX], vecB[NX], vecC[NX];
  double r=0.2;
  printf( "ici\n" );

/* Initialization of vectors */
  for( int i = 0; i < NX; i++ )
  {
    vecA[i] = pow( r, i );
    vecB[i] = 1.0;
  }
  printf( "ici\n" );

/* dot product of two vectors */
  //#pragma omp target teams distribute map(from:vecC[0:NX]) map(to:vecA[0:NX],vecB[0:NX])
  #pragma omp target teams distribute map(from:vecC) map(to:vecA,vecB)
  for( int i = 0; i < NX; i++ )
  {
    vecC[i] = vecA[i] * vecB[i];
  }

  double sum = 0.0;
  /* calculate the sum */
  #pragma omp target map(tofrom:sum)
  for( int i = 0; i < NX; i++ )
  {
    sum += vecC[i];
  }
  printf( "The sum is: %8.6f \n", sum );
  return 0;
}
