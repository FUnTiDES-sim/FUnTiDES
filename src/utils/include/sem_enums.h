#ifndef FUNTIDES_UTILS_INCLUDE_SEM_ENUMS_H_
#define FUNTIDES_UTILS_INCLUDE_SEM_ENUMS_H_
#include <string>

namespace utils {
namespace enums {

enum class methodType { kSem, kDg };
enum class implemType { kMakutu };
enum class meshType { kStruct, kUnstruct };
enum class modelLocationType { kOnNodes, kOnElements };
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

inline std::string to_string(methodType m) {
  switch (m) {
    case methodType::kSem:
      return "SEM";
    case methodType::kDg:
      return "DG";
    default:
      return "Unknown";
  }
}

inline std::string to_string(implemType i) {
  switch (i) {
    case implemType::kMakutu:
      return "MAKUTU";
    default:
      return "Unknown";
  }
}

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
#endif  // FUNTIDES_UTILS_INCLUDE_SEM_ENUMS_H_
