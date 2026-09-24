#include "model_struct.h"

/// Explicit instantiations of ModelStruct<float, int, Order> for Order = 1 to 9.
namespace model {
template class ModelStruct<float, int, 1>;
template class ModelStruct<float, int, 2>;
template class ModelStruct<float, int, 3>;
template class ModelStruct<float, int, 4>;
template class ModelStruct<float, int, 5>;
template class ModelStruct<float, int, 6>;
template class ModelStruct<float, int, 7>;
template class ModelStruct<float, int, 8>;
template class ModelStruct<float, int, 9>;
}  // namespace model
