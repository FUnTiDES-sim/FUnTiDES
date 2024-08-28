//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <chrono>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <omp.h>

#define MemSpace Kokkos::SharedSpace
using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy< ExecSpace >;

//#define h_MemSpace Kokkos::HostSpace
//using h_ExecSpace = h_MemSpace::execution_space;
//using h_range_policy = Kokkos::RangePolicy<h_ExecSpace>;

typedef Kokkos::LayoutLeft layout;
typedef Kokkos::View< int *, layout, MemSpace >   vectorInt;
typedef Kokkos::View< float *, layout, MemSpace >   vectorReal;
typedef Kokkos::View< int * *, layout, MemSpace >   arrayInt;
typedef Kokkos::View< float * *, layout, MemSpace >   arrayReal;
typedef Kokkos::View< double * * *, layout, MemSpace >   array3DReal;
typedef Kokkos::View< double * *, layout, MemSpace >   arrayDouble;
typedef Kokkos::View< double * * *, layout, MemSpace >   array3DDouble;


template< class T >
T allocateVector( int n1 )
{
  std::cout<<"allocate vector of size "<<n1<<std::endl;
  T vect( "v", n1 );
  return vect;
}
template< class T >
T allocateArray2D( int n1, int n2 )
{
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
  T array( "a", n1, n2 );
  return array;
}
template< class T >
T allocateArray3D( int n1, int n2, int n3 )
{
  std::cout<<"allocate array of size "<<n1<<" ,"<<n2<<" ,"<<n3<<std::endl;
  T array( "a", n1, n2, n3 );
  return array;
}

struct model{
             vectorReal coef;
             array3DReal pnp1;
             array3DReal pn;
             array3DReal pnm1;
            };
int main( int argc, char *argv[] )
{
  Kokkos::initialize( argc, argv );
  {
    const int n1=400;
    const int n2=400;
    const int n3=400;
    const int ncoefs=4;
    float dt=0.001;

    model m;

    //vectorReal coef;
    m.coef=allocateVector<vectorReal>(ncoefs);
    //array3DReal pnp1;
    m.pnp1=allocateArray3D<array3DReal>(n1,n2,n3);
    //array3DReal pn;
    m.pn=allocateArray3D<array3DReal>(n1,n2,n3);
    //array3DReal pnm1;
    m.pnm1=allocateArray3D<array3DReal>(n1,n2,n3);
    

 #pragma omp parallel for
    for( int k=0; k<n3; k++ )
    {
      for( int j=0; j<n2; j++ )
      {
        for( int i=0; i<n1; i++ )
        {
          m.pnp1(i,j,k)=1;
          m.pn(i,j,k)=1;
          m.pnm1(i,j,k)=1;
	}
      }
    }

    m.coef(0)=4*dt*dt;
    m.coef(1)=1*dt*dt;
    m.coef(2)=1*dt*dt;
    m.coef(3)=1*dt*dt;

    Kokkos::Timer timer1;
    for( int k=3; k<n3-3; k++ )
    {
      for( int j=3; j<n2-3; j++ )
      {
        for( int i=3; i<n1-3; i++ )
        {
		m.pnp1(i,j,k)=(2. +8*m.coef(0))*m.pn(i,j,k)+m.coef(1)*(m.pn(i+1,j,k)+m.pn(i-1,j,k)+m.pn(i,j+1,k)+m.pn(i,j-1,k)+m.pn(i,j,k+1)+m.pn(i,j,k-1))
                             +m.coef(2)*(m.pn(i+2,j,k)+m.pn(i-2,j,k)+m.pn(i,j+2,k)+m.pn(i,j-2,k)+m.pn(i,j,k+2)+m.pn(i,j,k-2))
                             +m.coef(3)*(m.pn(i+3,j,k)+m.pn(i-3,j,k)+m.pn(i,j+3,k)+m.pn(i,j-3,k)+m.pn(i,j,k+3)+m.pn(i,j,k-3))
			     -m.pnm1(i,j,k);
        }
      }
    }
    double time1=timer1.seconds();
    printf("result 1 %f\n",m.pnp1(n1/2,n2/2,n3/2));
    std::cout << "Elapsed Time sequential  loop : "<<time1 <<" seconds.\n\n";
    Kokkos::fence();

 #pragma omp parallel for
    for( int k=0; k<n3; k++ )
    {
      for( int j=0; j<n2; j++ )
      {
        for( int i=0; i<n1; i++ )
        {
          m.pnp1(i,j,k)=1;
          m.pn(i,j,k)=1;
          m.pnm1(i,j,k)=1;
	}
      }
    }

    Kokkos::Timer timer2;
 #pragma omp parallel for
    for( int k=3; k<n3-3; k++ )
    {
      for( int j=3; j<n2-3; j++ )
      {
        for( int i=3; i<n1-3; i++ )
        {
		m.pnp1(i,j,k)=(2. +8*m.coef(0))*m.pn(i,j,k)+m.coef(1)*(m.pn(i+1,j,k)+m.pn(i-1,j,k)+m.pn(i,j+1,k)+m.pn(i,j-1,k)+m.pn(i,j,k+1)+m.pn(i,j,k-1))
                             +m.coef(2)*(m.pn(i+2,j,k)+m.pn(i-2,j,k)+m.pn(i,j+2,k)+m.pn(i,j-2,k)+m.pn(i,j,k+2)+m.pn(i,j,k-2))
                             +m.coef(3)*(m.pn(i+3,j,k)+m.pn(i-3,j,k)+m.pn(i,j+3,k)+m.pn(i,j-3,k)+m.pn(i,j,k+3)+m.pn(i,j,k-3))
			     -m.pnm1(i,j,k);
        }
      }
    }
    double time2=timer2.seconds();
    printf("result 2 %f\n",m.pnp1(n1/2,n2/2,n3/2));
    std::cout << "Elapsed Time omp  loop : "<<time2 <<" seconds.\n\n";
    Kokkos::fence();

 #pragma omp parallel for
    for( int k=0; k<n3; k++ )
    {
      for( int j=0; j<n2; j++ )
      {
        for( int i=0; i<n1; i++ )
        {
          m.pnp1(i,j,k)=1;
          m.pn(i,j,k)=1;
          m.pnm1(i,j,k)=1;
	}
      }
    }
    Kokkos::Timer timer3;
 #pragma omp parallel for collapse(3)
    for( int k=3; k<n3-3; k++ )
    {
      for( int j=3; j<n2-3; j++ )
      {
        for( int i=3; i<n1-3; i++ )
        {
		m.pnp1(i,j,k)=(2. +8*m.coef(0))*m.pn(i,j,k)+m.coef(1)*(m.pn(i+1,j,k)+m.pn(i-1,j,k)+m.pn(i,j+1,k)+m.pn(i,j-1,k)+m.pn(i,j,k+1)+m.pn(i,j,k-1))
                             +m.coef(2)*(m.pn(i+2,j,k)+m.pn(i-2,j,k)+m.pn(i,j+2,k)+m.pn(i,j-2,k)+m.pn(i,j,k+2)+m.pn(i,j,k-2))
                             +m.coef(3)*(m.pn(i+3,j,k)+m.pn(i-3,j,k)+m.pn(i,j+3,k)+m.pn(i,j-3,k)+m.pn(i,j,k+3)+m.pn(i,j,k-3))
			     -m.pnm1(i,j,k);
    printf("result 3 %f\n",m.pnp1(n1/2,n2/2,n3/2));
    std::cout << "Elapsed Time omp collapse  loop : "<<time3 <<" seconds.\n\n";
    Kokkos::fence();

    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>> ({3,3,3},{n1,n2,n3}) , KOKKOS_LAMBDA ( int i,int j, int k)
        {
          m.pnp1(i,j,k)=1;
          m.pn(i,j,k)=1;
          m.pnm1(i,j,k)=1;
	});

    Kokkos::Timer timer4;
    Kokkos::parallel_for( n1, KOKKOS_LAMBDA ( int i )
    {
      for( int j=3; j<n2-3; j++ )
      {
        for( int k=3; k<n3-3; k++ )
        {
		m.pnp1(i,j,k)=(2. +8*m.coef(0))*m.pn(i,j,k)+m.coef(1)*(m.pn(i+1,j,k)+m.pn(i-1,j,k)+m.pn(i,j+1,k)+m.pn(i,j-1,k)+m.pn(i,j,k+1)+m.pn(i,j,k-1))
                             +m.coef(2)*(m.pn(i+2,j,k)+m.pn(i-2,j,k)+m.pn(i,j+2,k)+m.pn(i,j-2,k)+m.pn(i,j,k+2)+m.pn(i,j,k-2))
                             +m.coef(3)*(m.pn(i+3,j,k)+m.pn(i-3,j,k)+m.pn(i,j+3,k)+m.pn(i,j-3,k)+m.pn(i,j,k+3)+m.pn(i,j,k-3))
			     -m.pnm1(i,j,k);
        }
      }
    } );
    Kokkos::fence();
    double time4=timer4.seconds();
    printf("result 4 %f\n",m.pnp1(n1/2,n2/2,n3/2));
    std::cout << "Elapsed Time parallel_for kokkos  loop : "<<time4 <<" seconds.\n\n";
    
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>> ({0,0,0},{n1,n2,n3}) , KOKKOS_LAMBDA ( int i,int j, int k)
        {
          m.pnp1(i,j,k)=1;
          m.pn(i,j,k)=1;
          m.pnm1(i,j,k)=1;
	});

    Kokkos::Timer timer5;
    Kokkos::parallel_for( Kokkos::MDRangePolicy< Kokkos::Rank< 3 > >( {3, 3, 3}, {n1-3, n2-3, n3-3} ), KOKKOS_LAMBDA ( int i, int j, int k )
    {
		m.pnp1(i,j,k)=(2. +8*m.coef(0))*m.pn(i,j,k)+m.coef(1)*(m.pn(i+1,j,k)+m.pn(i-1,j,k)+m.pn(i,j+1,k)+m.pn(i,j-1,k)+m.pn(i,j,k+1)+m.pn(i,j,k-1))
                             +m.coef(2)*(m.pn(i+2,j,k)+m.pn(i-2,j,k)+m.pn(i,j+2,k)+m.pn(i,j-2,k)+m.pn(i,j,k+2)+m.pn(i,j,k-2))
                             +m.coef(3)*(m.pn(i+3,j,k)+m.pn(i-3,j,k)+m.pn(i,j+3,k)+m.pn(i,j-3,k)+m.pn(i,j,k+3)+m.pn(i,j,k-3))
			     -m.pnm1(i,j,k);
    });
    Kokkos::fence();
    double time5=timer5.seconds();
    printf("result 5 %f\n",m.pnp1(n1/2,n2/2,n3/2));
    std::cout << "Elapsed Time parallel_for MDRANGE kokkos  loop : "<<time5 <<" seconds.\n\n";

  }
  Kokkos::finalize();
  return (0);
}
