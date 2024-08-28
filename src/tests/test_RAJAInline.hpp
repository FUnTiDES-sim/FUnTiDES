#ifndef TEST_RAJAINLINE_HPP_
#define TEST_RAJAINLINE_HPP_
#include <vector>
#include "Array.hpp"
#include "RAJA/RAJA.hpp"
#include "Macros.hpp"
#include "ChaiBuffer.hpp"

// test lvArray
using vectorInt=LvArray::Array< int,
                                1,
                                camp::idx_seq< 0 >,
                                std::ptrdiff_t,
                                LvArray::ChaiBuffer >;
using vectorIntView=LvArray::ArrayView< int,
                                        1,
                                        0,
                                        std::ptrdiff_t,
                                        LvArray::ChaiBuffer >;
using arrayInt=LvArray::Array< int,
                               2,
                               camp::idx_seq< 1, 0 >,
                               std::ptrdiff_t,
                               LvArray::ChaiBuffer >;
using arrayIntView=LvArray::ArrayView< int,
                                       2,
                                       0,
                                       std::ptrdiff_t,
                                       LvArray::ChaiBuffer >;
using arrayReal=LvArray::Array< float,
                                2,
                                camp::idx_seq< 1, 0 >,
                                std::ptrdiff_t,
                                LvArray::ChaiBuffer >;
using arrayRealView=LvArray::ArrayView< float,
                                        2,
                                        0,
                                        std::ptrdiff_t,
                                        LvArray::ChaiBuffer >;
using array2DDouble=LvArray::Array< double,
                                    2,
                                    camp::idx_seq< 1, 0 >,
                                    std::ptrdiff_t,
                                    LvArray::ChaiBuffer >;
using array2DDoubleView=LvArray::ArrayView< double,
                                            2,
                                            0,
                                            std::ptrdiff_t,
                                            LvArray::ChaiBuffer >;
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

class mesh
{

public:
  LVARRAY_HOST_DEVICE mesh();
  LVARRAY_HOST_DEVICE  ~mesh();

  LVARRAY_HOST_DEVICE int localToGlobalNodes( const int & e,
                                              const int & nPointsPerElement, arrayIntView const & nodesList,
                                              int localToGlobal[] )const
  {
    //vectorInt localToGlobal( nPointsPerElement );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      localToGlobal[i]=nodesList( e, i );
      //if(e==0)printf("i=%d %d\n",i,localToGlobal(e,i));
    }
    //return localToGlobal;
    return 0;
  }

  LVARRAY_HOST_DEVICE int getXi( const int & numberOfPointsPerElement, arrayRealView const & globalNodesCoords,
                                 int const localToGlobal[], double Xi[][2] )const
  {
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      Xi[i][0]=globalNodesCoords( localToGlobal[i], 0 );
      Xi[i][1]=globalNodesCoords( localToGlobal[i], 1 );
    }
    return 0;
  }
};


#endif
