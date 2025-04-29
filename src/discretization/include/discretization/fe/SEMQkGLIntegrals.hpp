#ifdef USE_SEMCLASSIC
    #include <discretization/fe/SEMQkGLIntegralsClassic.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <discretization/fe/SEMQkGLIntegralsOptim.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsOptim;
#endif
#ifdef USE_SHIVA
    #include <discretization/fe/SEMQkGLIntegralsShiva.hpp>
    using SEMQkGLIntegrals = SEMQkGLIntegralsShiva ;
#endif
