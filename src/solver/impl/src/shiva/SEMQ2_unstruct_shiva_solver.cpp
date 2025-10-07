#include <model_unstruct.h>

#include "sem_solver_impl.h"

constexpr int ORDER = 2;
template class SEMsolver<ORDER,
                         IntegralTypeSelector<ORDER, false, IntegralType::SHIVA>::type,
                         model::ModelUnstruct<float, int> >;
