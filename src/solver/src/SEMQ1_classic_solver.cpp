#include "SEMsolver_impl.hpp"

constexpr int ORDER = 1;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type, model::ModelStruct<float, int> >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type, model::ModelUnstruct<float, int> >;
