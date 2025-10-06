#include "model_unstruct.h"

namespace model
{
namespace mesh
{
template class ModelUnstruct<float, int>;
template class ModelUnstruct<float, long>;
template class ModelUnstruct<double, int>;
template class ModelUnstruct<double, long>;
}  // namespace mesh
}  // namespace model
