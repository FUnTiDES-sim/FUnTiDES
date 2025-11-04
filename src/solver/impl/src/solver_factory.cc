#include "solver_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver.h"

namespace SolverFactory
{

template <typename FUNC>
std::unique_ptr<SolverBase> orderDispatch(int const order, FUNC&& func)
{
  if (order == 1)
  {
    return func(integral_constant<int, 1>{});
  }
  if (order == 2)
  {
    return func(integral_constant<int, 2>{});
  }
  else if (order == 3)
  {
    return func(integral_constant<int, 3>{});
  }
  // else if( order == 4 )
  // {
  //   return func( integral_constant< int, 4>{} );
  // }
  abort();
}

template <auto ImplTag>
static std::unique_ptr<SolverBase> make_sem_solver(
    int order, meshType mesh, modelLocationType modelLocation)
{
  bool isModelOnNodes = (modelLocation == OnNodes);

  switch (mesh)
  {
    case Struct:
      return orderDispatch(
          order, [isModelOnNodes](auto orderIC) -> std::unique_ptr<SolverBase> {
            constexpr int ORDER = decltype(orderIC)::value;
            using SelectedIntegral =
                typename IntegralTypeSelector<ORDER, ImplTag>::type;
            using MeshT = model::ModelStruct<float, int, ORDER>;

            if (isModelOnNodes)
            {
              return std::make_unique<
                  SEMsolver<ORDER, SelectedIntegral, MeshT, true>>();
            }
            else
            {
              return std::make_unique<
                  SEMsolver<ORDER, SelectedIntegral, MeshT, false>>();
            }
          });
    case Unstruct:
      return orderDispatch(
          order, [isModelOnNodes](auto orderIC) -> std::unique_ptr<SolverBase> {
            constexpr int ORDER = decltype(orderIC)::value;
            using SelectedIntegral =
                typename IntegralTypeSelector<ORDER, ImplTag>::type;
            using MeshT = model::ModelUnstruct<float, int>;

            if (isModelOnNodes)
            {
              return std::make_unique<
                  SEMsolver<ORDER, SelectedIntegral, MeshT, true>>();
            }
            else
            {
              return std::make_unique<
                  SEMsolver<ORDER, SelectedIntegral, MeshT, false>>();
            }
          });
  }
  throw std::runtime_error("Unknown mesh type");
}

std::unique_ptr<SolverBase> createSolver(methodType const methodType,
                                         implemType const implemType,
                                         meshType const mesh,
                                         modelLocationType const modelLocation,
                                         int const order)
{
  if (methodType == SEM)
  {
    switch (implemType)
    {
      case MAKUTU:
        return make_sem_solver<IntegralType::MAKUTU>(order, mesh,
                                                     modelLocation);
      case SHIVA:
        return make_sem_solver<IntegralType::SHIVA>(order, mesh, modelLocation);
    }
  }

  // Add DG or other methods as needed
  throw std::runtime_error(
      "Unsupported solver configuration: methodType=" + to_string(methodType) +
      ", implemType=" + to_string(implemType));
}
}  // namespace SolverFactory
