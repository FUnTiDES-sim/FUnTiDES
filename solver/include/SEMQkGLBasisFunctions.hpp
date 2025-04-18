#ifndef SEMQKGLBASISFUNCTIONS_HPP_
#define SEMQKGLBASISFUNCTIONS_HPP_
#include <utils/dataType.hpp>
#include <utils/SEMmacros.hpp>
#include <mesh/SEMdata.hpp>
using namespace std;
#ifdef USE_SEMCLASSIC
    #include <SEMQkGLBasisFunctionsClassic.hpp> 
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#ifdef  USE_SHIVA
    #include <SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#endif //SEMQKGLBASISFUNCTIONS_HPP_
