#ifndef FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_UTILS_H_
#define FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_UTILS_H_

#pragma once

#include <string>

namespace model {

/**
 * @brief Short suffix naming a (FloatType, ScalarType) pair in Python class names.
 *
 * The primary template is only declared: a pair without a specialization
 * (float/double with int/long are provided) does not compile.
 * @tparam FloatType Floating-point type of the bound class.
 * @tparam ScalarType Integer type of the bound class.
 * @return Suffix such as "f32_i32", a string literal with static storage.
 */
template <typename FloatType, typename ScalarType>
constexpr const char* type_suffix();
template <>
constexpr const char* type_suffix<float, int>() {
  return "f32_i32";
}
template <>
constexpr const char* type_suffix<double, int>() {
  return "f64_i32";
}
template <>
constexpr const char* type_suffix<float, long>() {
  return "f32_i64";
}
template <>
constexpr const char* type_suffix<double, long>() {
  return "f64_i64";
}

/**
 * @brief Short suffix naming a polynomial order in Python class names.
 * @param[in] order Polynomial order, from 1 to 9.
 * @return Suffix such as "O3", a string literal with static storage.
 * @throws std::runtime_error If order is outside 1..9.
 */
constexpr const char* order_suffix(int order) {
  switch (order) {
    case 1:
      return "O1";
    case 2:
      return "O2";
    case 3:
      return "O3";
    case 4:
      return "O4";
    case 5:
      return "O5";
    case 6:
      return "O6";
    case 7:
      return "O7";
    case 8:
      return "O8";
    case 9:
      return "O9";
    default:
      throw std::runtime_error("Unsupported order for binding: " + std::to_string(order));
  }
}

/**
 * @brief Calls a functor with the run-time order turned into a compile-time constant.
 * @tparam FUNC Callable taking a std::integral_constant<int, Order>.
 * @param[in] order Polynomial order, from 1 to 9.
 * @param[in] func Callable invoked once. Its return type must be the same for every order.
 * @return The value returned by func.
 * @throws std::invalid_argument If order is outside 1..9.
 */
template <typename FUNC>
auto orderDispatch(int const order, FUNC&& func) {
  switch (order) {
    case 1:
      return func(std::integral_constant<int, 1>{});
    case 2:
      return func(std::integral_constant<int, 2>{});
    case 3:
      return func(std::integral_constant<int, 3>{});
    case 4:
      return func(std::integral_constant<int, 4>{});
    case 5:
      return func(std::integral_constant<int, 5>{});
    case 6:
      return func(std::integral_constant<int, 6>{});
    case 7:
      return func(std::integral_constant<int, 7>{});
    case 8:
      return func(std::integral_constant<int, 8>{});
    case 9:
      return func(std::integral_constant<int, 9>{});
    default:
      throw std::invalid_argument("Unsupported order for binding: " + std::to_string(order));
  }
}

/**
 * @brief Python class name for a class specialized on types and order.
 * @tparam FloatType Floating-point type of the bound class.
 * @tparam ScalarType Integer type of the bound class.
 * @tparam Order Polynomial order, from 1 to 9.
 * @param[in] basename Class name without suffix.
 * @return basename followed by the type suffix and the order suffix, joined by "_".
 * @throws std::runtime_error If Order is outside 1..9.
 */
template <typename FloatType, typename ScalarType, int Order>
std::string model_class_name(std::string basename) {
  return basename + "_" + type_suffix<FloatType, ScalarType>() + "_" + order_suffix(Order);
}

/**
 * @brief Python class name for a class specialized on types only.
 * @tparam FloatType Floating-point type of the bound class.
 * @tparam ScalarType Integer type of the bound class.
 * @param[in] basename Class name without suffix.
 * @return basename followed by the type suffix, joined by "_".
 */
template <typename FloatType, typename ScalarType>
std::string model_class_name(std::string basename) {
  return basename + "_" + type_suffix<FloatType, ScalarType>();
}

}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_UTILS_H_
