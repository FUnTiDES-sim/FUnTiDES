#ifndef SEMQKGLBASISFUNCTIONS_HPP_
#define SEMQKGLBASISFUNCTIONS_HPP_
#include <utils/dataType.hpp>
#include <utils/SEMmacros.hpp>
#include <utils/SEMdata.hpp>
using namespace std;
#ifdef USE_SEMCLASSIC
    #include <discretization/SEMQkGLBasisFunctionsClassic.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <discretization/SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#ifdef  USE_SHIVA
    #include <discretization/SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#endif //SEMQKGLBASISFUNCTIONS_HPP_
