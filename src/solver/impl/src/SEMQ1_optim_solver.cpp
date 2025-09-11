#include <model_struct.h>
#include <model_unstruct.h>
#include "SEMsolver_impl.hpp"

constexpr int ORDER = 1;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM>::type, model::ModelStruct<float, int, ORDER> >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM>::type, model::ModelUnstruct<float, int> >;
