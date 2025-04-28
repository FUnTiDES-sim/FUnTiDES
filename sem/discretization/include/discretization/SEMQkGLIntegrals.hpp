#ifdef USE_SEMCLASSIC
    #include <discretization/SEMQkGLIntegralsClassic.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <discretization/SEMQkGLIntegralsOptim.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsOptim;
#endif
#ifdef USE_SHIVA
    #include <discretization/SEMQkGLIntegralsShiva.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsShiva ;
#endif
