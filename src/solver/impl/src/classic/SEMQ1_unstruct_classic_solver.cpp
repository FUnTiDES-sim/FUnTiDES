#include <model_unstruct.h>

#include "SEMsolver_impl.hpp"

constexpr int ORDER = 1;
template class SEMsolver<
    ORDER, IntegralTypeSelector<ORDER, IntegralType::CLASSIC>::type,
    model::ModelUnstruct<float, int> >;
