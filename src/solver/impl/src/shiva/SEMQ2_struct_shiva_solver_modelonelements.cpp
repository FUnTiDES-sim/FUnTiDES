#include <model_struct.h>

#include "sem_solver_impl.h"

constexpr int ORDER = 2;
template class SEMsolver<ORDER,
                         IntegralTypeSelector<ORDER, IntegralType::SHIVA>::type,
                         model::ModelStruct<float, int, ORDER>, false >;
