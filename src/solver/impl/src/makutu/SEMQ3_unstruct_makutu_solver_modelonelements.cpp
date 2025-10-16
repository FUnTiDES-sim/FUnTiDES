#include <model_unstruct.h>

#include "sem_solver_impl.h"

constexpr int ORDER = 3;
template class SEMsolver<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type,
    model::ModelUnstruct<float, int>, false>;
