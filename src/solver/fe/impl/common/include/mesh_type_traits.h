#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MESH_TYPE_TRAITS_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MESH_TYPE_TRAITS_H_
#include <sem_enums.h>

#include <type_traits>

#include "model_struct.h"
#include "model_unstruct.h"

namespace solver {
namespace fe {

/// @brief Maps a mesh type to its utils::enums::meshType value.
/// Default: kUnstruct. Specialise for each concrete mesh class.
template <typename T>
struct MeshTypeTraits : std::integral_constant<utils::enums::meshType, utils::enums::meshType::kUnstruct> {};

template <typename F, typename S, int Order>
struct MeshTypeTraits<model::ModelStruct<F, S, Order>>
    : std::integral_constant<utils::enums::meshType, utils::enums::meshType::kStruct> {};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MESH_TYPE_TRAITS_H_
