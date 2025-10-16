#include <model_unstruct.h>

#include "sem_solver_acoustic_impl.h"

constexpr int ORDER = 3;
template class SEMsolverAcoustic<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type,
    model::ModelUnstruct<float, int>, true>;
