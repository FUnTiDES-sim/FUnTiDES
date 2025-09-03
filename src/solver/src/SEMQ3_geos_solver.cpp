#include "SEMsolver_impl.hpp"

constexpr int ORDER = 3;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::GEOS >::type, model::ModelStruct< float, int, ORDER > >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::GEOS >::type, model::ModelUnstruct< float, int > >;
