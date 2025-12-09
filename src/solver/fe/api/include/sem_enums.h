#pragma once

#include <string>

namespace solver
{
namespace fe
{

enum methodType
{
  kSem,
  kDg
};
enum implemType
{
  kMakutu,
  kShiva
};
enum meshType
{
  kStruct,
  kUnstruct
};
enum modelLocationType
{
  kOnNodes,
  kOnElements
};
enum physicType
{
  kAcoustic,
  kElastic
};

inline std::string to_string(methodType m)
{
  switch (m)
  {
    case kSem:
      return "SEM";
    case kDg:
      return "DG";
    default:
      return "Unknown";
  }
}

inline std::string to_string(implemType i)
{
  switch (i)
  {
    case kMakutu:
      return "MAKUTU";
    case kShiva:
      return "SHIVA";
    default:
      return "Unknown";
  }
}

inline std::string to_string(meshType m)
{
  switch (m)
  {
    case kStruct:
      return "Struct";
    case kUnstruct:
      return "Unstruct";
    default:
      return "Unknown";
  }
}

inline std::string to_string(modelLocationType loc)
{
  switch (loc)
  {
    case kOnNodes:
      return "OnNodes";
    case kOnElements:
      return "OnElements";
    default:
      return "Unknown";
  }
}

inline std::string to_string(physicType p)
{
  switch (p)
  {
    case kAcoustic:
      return "Acoustic";
    case kElastic:
      return "Elastic";
    default:
      return "Unknown";
  }
}

}  // namespace fe
}  // namespace solver
