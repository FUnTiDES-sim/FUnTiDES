#include "SEMdata.hpp"

using tfloat = float;
using gfloat = float;

#ifdef USE_SEMCLASSIC
    #include <fe/SEMKernels/src/finiteElement/classic/SEMQkGLIntegralsClassic.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic<SEMinfo::myOrderNumber> ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <fe/SEMKernels/src/finiteElement/optim/SEMQkGLIntegralsOptim.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsOptim<SEMinfo::myOrderNumber, tfloat, gfloat>;
#endif
#ifdef USE_SHIVA
    #include <fe/SEMKernels/src/finiteElement/shiva/SEMQkGLIntegralsShiva.hpp>
    using TransformType =
    LinearTransform< tfloat,
                     InterpolatedShape< tfloat,
                                        Cube< tfloat >,
                                        LagrangeBasis< tfloat, 1, EqualSpacing >,
                                        LagrangeBasis< tfloat, 1, EqualSpacing >,
                                        LagrangeBasis< tfloat, 1, EqualSpacing > > >;

  using ParentElementType =
    ParentElement< gfloat,
                   Cube< gfloat >,
                   LagrangeBasis< gfloat, SEMinfo::myOrderNumber, EqualSpacing >,
                   LagrangeBasis< gfloat, SEMinfo::myOrderNumber, EqualSpacing >,
                   LagrangeBasis< gfloat, SEMinfo::myOrderNumber, EqualSpacing > >;

  using SEMQkGLIntegrals = SEMQkGLIntegralsShiva< SEMinfo::myOrderNumber, TransformType, ParentElementType >;    
#endif
