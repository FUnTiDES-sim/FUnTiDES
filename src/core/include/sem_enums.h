#ifndef FUNTIDES_CORE_INCLUDE_SEM_ENUMS_H_
#define FUNTIDES_CORE_INCLUDE_SEM_ENUMS_H_
#include <string>

namespace utils {
namespace enums {

/// @brief Discretization method used by the solver.
enum class methodType { kSem, kDg, kDgSem, kDgPAdaptive };

/// @brief Runtime selector of the finite-element back-end.
enum class implemType { kMakutu };

/// @brief Kind of mesh: structured or unstructured.
enum class meshType { kStruct, kUnstruct };

/// @brief Where the model properties are stored.
enum class modelLocationType { kOnNodes, kOnElements };

/// @brief Physics solved by the simulation.
enum class physicType : int { kAcoustic, kElastic, kAcoustoElastic };

/// @brief How the mesh builder assigned material properties to the nodes that
/// sit on an acoustic/elastic interface.
enum class interfacePropertyConvention {
  /// Interface nodes carry the fluid state; the solid state is rebuilt from an
  /// adjacent elastic element.
  kFluidOnInterfaceNodes,
  /// Interface nodes carry a single state, meant to be used by both sides.
  kSharedOnInterfaceNodes
};

/// @brief Human-readable name of a method type ("Unknown" if out of range).
inline std::string to_string(methodType m) {
  switch (m) {
    case methodType::kSem:
      return "SEM";
    case methodType::kDg:
      return "DG";
    case methodType::kDgSem:
      return "DG-SEM";
    case methodType::kDgPAdaptive:
      return "DGPAdaptive";
    default:
      return "Unknown";
  }
}

/// @brief Human-readable name of an implementation type ("Unknown" if out of range).
inline std::string to_string(implemType i) {
  switch (i) {
    case implemType::kMakutu:
      return "MAKUTU";
    default:
      return "Unknown";
  }
}

/// @brief Human-readable name of a mesh type ("Unknown" if out of range).
inline std::string to_string(meshType m) {
  switch (m) {
    case meshType::kStruct:
      return "Struct";
    case meshType::kUnstruct:
      return "Unstruct";
    default:
      return "Unknown";
  }
}

/// @brief Human-readable name of a model location ("Unknown" if out of range).
inline std::string to_string(modelLocationType loc) {
  switch (loc) {
    case modelLocationType::kOnNodes:
      return "OnNodes";
    case modelLocationType::kOnElements:
      return "OnElements";
    default:
      return "Unknown";
  }
}

/// @brief Human-readable name of a physics type ("Unknown" if out of range).
inline std::string to_string(physicType p) {
  switch (p) {
    case physicType::kAcoustic:
      return "Acoustic";
    case physicType::kElastic:
      return "Elastic";
    case physicType::kAcoustoElastic:
      return "AcoustoElastic";
    default:
      return "Unknown";
  }
}

/// @brief Human-readable name of an interface property convention ("Unknown" if out of range).
inline std::string to_string(interfacePropertyConvention c) {
  switch (c) {
    case interfacePropertyConvention::kFluidOnInterfaceNodes:
      return "FluidOnInterfaceNodes";
    case interfacePropertyConvention::kSharedOnInterfaceNodes:
      return "SharedOnInterfaceNodes";
    default:
      return "Unknown";
  }
}

}  // namespace enums
}  // namespace utils
#endif  // FUNTIDES_CORE_INCLUDE_SEM_ENUMS_H_
