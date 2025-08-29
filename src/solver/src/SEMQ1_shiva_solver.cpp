#include "SEMsolver_impl.hpp"

constexpr int ORDER = 1;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::SHIVA>::type, model::ModelStruct<float, int, ORDER> >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::SHIVA>::type, model::ModelUnstruct<float, int> >;
