#include <model_struct.h>
#include "SEMsolver_impl.hpp"

constexpr int ORDER = 2;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type, model::ModelStruct<float, int, ORDER> >;
