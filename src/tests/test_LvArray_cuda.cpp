// Source includes
#include "Array.hpp"
#include "MallocBuffer.hpp"
#if defined(LVARRAY_USE_CHAI)
  #include "ChaiBuffer.hpp"
#endif

// TPL includes
#include <RAJA/RAJA.hpp>

// system includes
#include <string>

int main( int argc, char *argv[] )
{
  LvArray::Array< int,
                  2,
                  camp::idx_seq< 1, 0 >,
                  std::ptrdiff_t,
                  LvArray::ChaiBuffer > array( 5, 6 );

  // Move the array to the device.
  array.move( LvArray::MemorySpace::cuda );
  int * const devicePointer = array.data();

  RAJA::forall< RAJA::cuda_exec< 32 > >(
    RAJA::TypedRangeSegment< std::ptrdiff_t >( 0, array.size()),
    [devicePointer] __device__ ( std::ptrdiff_t const i )
  {
    devicePointer[ i ] = i;
  }
    );

  LvArray::ArrayView< int,
                      2,
                      0,
                      std::ptrdiff_t,
                      LvArray::ChaiBuffer > const & view = array;

  // Capture the view in a host kernel which moves the data back to the host.
  RAJA::forall< RAJA::loop_exec >(
    RAJA::TypedRangeSegment< std::ptrdiff_t >( 0, view.size()),
    [view] ( std::ptrdiff_t const i )
  {
    //EXPECT_EQ( view.data()[ i ], i );
  }
    );
}
