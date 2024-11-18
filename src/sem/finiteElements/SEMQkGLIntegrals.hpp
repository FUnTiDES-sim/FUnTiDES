#ifdef USE_SEMCLASSIC
    #include <SEMQkGLIntegralsClassic.hpp> 
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <SEMQkGLIntegralsOptim.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsOptim;
#endif
#ifdef USE_SHIVA
    #include <SEMQkGLIntegralsShiva.hpp> 
    using SEMQkGLIntegrals = SEMQkGLIntegralsShiva ;
#endif
