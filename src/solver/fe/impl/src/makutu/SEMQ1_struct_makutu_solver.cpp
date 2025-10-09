#include <model_struct.h>

#include "sem_solver_impl.h"

constexpr int ORDER = 1;
template class SEMsolver<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type,
    model::mesh::ModelStruct<float, int, ORDER> >;
