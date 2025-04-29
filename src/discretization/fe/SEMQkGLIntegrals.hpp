#ifdef USE_SEMCLASSIC
    #include <fe/SEMQkGLIntegralsClassic.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <fe/SEMQkGLIntegralsOptim.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsOptim;
#endif
#ifdef USE_SHIVA
    #include <fe/SEMQkGLIntegralsShiva.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsShiva ;
#endif
