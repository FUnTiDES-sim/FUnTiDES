#include "SEMdata.hpp"


#ifdef USE_SEMCLASSIC
    #include <fe/SEMKernels/src/finiteElement/SEMQkGLIntegralsClassic.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic<SEMinfo::myOrderNumber> ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <fe/SEMKernels/src/finiteElement/SEMQkGLIntegralsOptim.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsOptim<SEMinfo::myOrderNumber>;
#endif
#ifdef USE_SHIVA
    #include <fe/SEMKernels/src/finiteElement/SEMQkGLIntegralsShiva.hpp>
    using TransformType =
    LinearTransform< double,
                     InterpolatedShape< double,
                                        Cube< double >,
                                        LagrangeBasis< double, 1, EqualSpacing >,
                                        LagrangeBasis< double, 1, EqualSpacing >,
                                        LagrangeBasis< double, 1, EqualSpacing > > >;

  using ParentElementType =
    ParentElement< double,
                   Cube< double >,
                   LagrangeBasis< double, SEMinfo::myOrderNumber, EqualSpacing >,
                   LagrangeBasis< double, SEMinfo::myOrderNumber, EqualSpacing >,
                   LagrangeBasis< double, SEMinfo::myOrderNumber, EqualSpacing > >;

  using SEMQkGLIntegrals = SEMQkGLIntegralsShiva< double, SEMinfo::myOrderNumber, TransformType, ParentElementType >;    
#endif
