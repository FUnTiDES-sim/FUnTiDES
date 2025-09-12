#include <model_unstruct.h>
#include "SEMsolver_impl.hpp"

constexpr int ORDER = 2;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM>::type, model::ModelUnstruct<float, int> >;
