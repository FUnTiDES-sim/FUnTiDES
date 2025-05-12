#include "SEMdata.hpp"


#ifdef USE_SEMCLASSIC
    #include <fe/SEMQkGLIntegralsClassic.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic ;
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

  constexpr int order = 1;
  using ParentElementType =
    ParentElement< double,
                   Cube< double >,
                   LagrangeBasis< double, order, EqualSpacing >,
                   LagrangeBasis< double, order, EqualSpacing >,
                   LagrangeBasis< double, order, EqualSpacing > >;

  using SEMQkGLIntegrals = SEMQkGLIntegralsShiva< double, order, TransformType, ParentElementType >;    
#endif
