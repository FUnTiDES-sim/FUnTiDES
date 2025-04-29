#ifndef SEMQKGLBASISFUNCTIONS_HPP_
#define SEMQKGLBASISFUNCTIONS_HPP_
#include <utils/dataType.hpp>
#include <utils/SEMmacros.hpp>
#include <utils/SEMdata.hpp>
using namespace std;
#ifdef USE_SEMCLASSIC
    #include <discretization/fe/SEMQkGLBasisFunctionsClassic.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <discretization/fe/SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#ifdef  USE_SHIVA
    #include <discretization/fe/SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#endif //SEMQKGLBASISFUNCTIONS_HPP_
