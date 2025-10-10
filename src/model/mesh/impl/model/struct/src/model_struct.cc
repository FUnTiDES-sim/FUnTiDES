#include "model_struct.h"

namespace model
{
namespace mesh
{
template class ModelStruct<float, int, 1>;
template class ModelStruct<float, int, 2>;
template class ModelStruct<float, int, 3>;
template class ModelStruct<float, int, 4>;
}  // namespace mesh
}  // namespace model
