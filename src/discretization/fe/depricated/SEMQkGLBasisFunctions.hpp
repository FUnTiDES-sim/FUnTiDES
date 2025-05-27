#ifndef SEMQKGLBASISFUNCTIONS_HPP_
#define SEMQKGLBASISFUNCTIONS_HPP_
#include <dataType.hpp>
#include <SEMmacros.hpp>
#include <SEMdata.hpp>
using namespace std;
#ifdef USE_SEMCLASSIC
    #include <fe/SEMQkGLBasisFunctionsClassic.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsClassic ;
#endif
#ifdef  USE_SEMOPTIM 
    #include <fe/SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#ifdef  USE_SHIVA
    #include <fe/SEMQkGLBasisFunctionsOptim.hpp>
    using SEMQkGLBasisFunctions = SEMQkGLBasisFunctionsOptim;
#endif
#endif //SEMQKGLBASISFUNCTIONS_HPP_
