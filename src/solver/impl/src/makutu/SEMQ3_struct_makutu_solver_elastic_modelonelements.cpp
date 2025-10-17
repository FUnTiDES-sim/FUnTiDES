#include <model_struct.h>

#include "sem_solver_elastic_impl.h"

constexpr int ORDER = 3;
template class SEMsolverElastic<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type,
    model::ModelStruct<float, int, ORDER>, false>;
