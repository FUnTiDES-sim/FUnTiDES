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
#ifdef USE_SEMGEOS
    #include <fe/SEMQkGLIntegralsGeos.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsGeos;
#endif // USE_SEMGEOS
