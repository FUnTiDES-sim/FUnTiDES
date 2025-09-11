#include <model_unstruct.h>

#include "SEMsolver_impl.hpp"

constexpr int ORDER = 3;
template class SEMsolver<ORDER,
                         IntegralTypeSelector<ORDER, IntegralType::GEOS>::type,
                         model::ModelUnstruct<float, int> >;
