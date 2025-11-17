#pragma once

#include <string>

namespace model
{

// helper to make a readable suffix for class names (specialize as needed)
template <typename FloatType, typename ScalarType>
constexpr const char* type_suffix();
template <>
constexpr const char* type_suffix<float, int>()
{
  return "f32_i32";
}
template <>
constexpr const char* type_suffix<double, int>()
{
  return "f64_i32";
}
template <>
constexpr const char* type_suffix<float, long>()
{
  return "f32_i64";
}
template <>
constexpr const char* type_suffix<double, long>()
{
  return "f64_i64";
}

// helper to make a readable suffix for order
constexpr const char* order_suffix(int order)
{
  switch (order)
  {
    case 1:
      return "O1";
    case 2:
      return "O2";
    case 3:
      return "O3";
    default:
      throw std::runtime_error("Unsupported order for binding: " +
                               std::to_string(order));
  }
}

template <typename FUNC>
auto orderDispatch(int const order, FUNC&& func)
{
  switch (order)
  {
    case 1:
      return func(std::integral_constant<int, 1>{});
    case 2:
      return func(std::integral_constant<int, 2>{});
    case 3:
      return func(std::integral_constant<int, 3>{});
    default:
      throw std::invalid_argument("Unsupported order for binding: " +
                                  std::to_string(order));
  }
}

// helper to generate class name
template <typename FloatType, typename ScalarType, int Order>
std::string model_class_name(std::string basename)
{
  return basename + "_" + type_suffix<FloatType, ScalarType>() + "_" +
         order_suffix(Order);
}

// helper to generate class name
template <typename FloatType, typename ScalarType>
std::string model_class_name(std::string basename)
{
  return basename + "_" + type_suffix<FloatType, ScalarType>();
}

}  // namespace model