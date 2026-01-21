#ifndef SOLVER_FE_IMPL_INCLUDE_SEM_ENUMS_H_
#define SOLVER_FE_IMPL_INCLUDE_SEM_ENUMS_H_

#include <string>

namespace solver
{
namespace fe
{
namespace enums
{

enum class methodType
{
  kSem,
  kDg
};
enum class implemType
{
  kMakutu,
  kShiva
};
enum class meshType
{
  kStruct,
  kUnstruct
};
enum class modelLocationType
{
  kOnNodes,
  kOnElements
};
enum class physicType : int
{
  kAcoustic,
  kElastic
};

inline std::string to_string(methodType m)
{
  switch (m)
  {
    case methodType::kSem:
      return "SEM";
    case methodType::kDg:
      return "DG";
    default:
      return "Unknown";
  }
}

inline std::string to_string(implemType i)
{
  switch (i)
  {
    case implemType::kMakutu:
      return "MAKUTU";
    case implemType::kShiva:
      return "SHIVA";
    default:
      return "Unknown";
  }
}

inline std::string to_string(meshType m)
{
  switch (m)
  {
    case meshType::kStruct:
      return "Struct";
    case meshType::kUnstruct:
      return "Unstruct";
    default:
      return "Unknown";
  }
}

inline std::string to_string(modelLocationType loc)
{
  switch (loc)
  {
    case modelLocationType::kOnNodes:
      return "OnNodes";
    case modelLocationType::kOnElements:
      return "OnElements";
    default:
      return "Unknown";
  }
}

inline std::string to_string(physicType p)
{
  switch (p)
  {
    case physicType::kAcoustic:
      return "Acoustic";
    case physicType::kElastic:
      return "Elastic";
    default:
      return "Unknown";
  }
}

}  // namespace enums
}  // namespace fe
}  // namespace solver
#endif