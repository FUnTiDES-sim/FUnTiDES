#include <model_struct.h>

#include "sem_solver_impl.h"

constexpr int ORDER = 3;
template class SEMsolver<ORDER,
                         IntegralTypeSelector<ORDER, IntegralType::GEOS>::type,
                         model::ModelStruct<float, int, ORDER> >;
