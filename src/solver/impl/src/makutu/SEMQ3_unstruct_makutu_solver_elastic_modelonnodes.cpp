#include <model_unstruct.h>

#include "sem_solver_elastic_impl.h"

constexpr int ORDER = 3;
template class SEMsolverElastic<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type,
    model::ModelUnstruct<float, int>, true>;
