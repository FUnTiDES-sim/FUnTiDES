#include <model_struct.h>

#include "sem_solver_acoustic_impl.h"

constexpr int ORDER = 2;
template class SEMsolverAcoustic<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type,
    model::ModelStruct<float, int, ORDER>, true>;
