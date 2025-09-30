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

// ImplTag is one of IntegralType::MAKUTU/SHIVA (compile-time)
template <auto ImplTag>
static std::unique_ptr<SolverBase> make_sem_solver(int order, meshType mesh)
{
  switch (mesh)
  {
    case Struct:
      // ORDER-dependent mesh type
      return orderDispatch(
          order, [](auto orderIC) -> std::unique_ptr<SolverBase> {
            constexpr int ORDER = decltype(orderIC)::value;
            using SelectedIntegral =
                typename IntegralTypeSelector<ORDER, ImplTag>::type;
            using MeshT = model::ModelStruct<float, int, ORDER>;
            return std::make_unique<
                SEMsolver<ORDER, SelectedIntegral, MeshT>>();
          });
    case Unstruct:
      return orderDispatch(
          order, [](auto orderIC) -> std::unique_ptr<SolverBase> {
            constexpr int ORDER = decltype(orderIC)::value;
            using SelectedIntegral =
                typename IntegralTypeSelector<ORDER, ImplTag>::type;
            using MeshT = model::ModelUnstruct<float, int>;
            return std::make_unique<
                SEMsolver<ORDER, SelectedIntegral, MeshT>>();
          });
  }

  throw std::runtime_error("Unknown mesh type");
}

std::unique_ptr<SolverBase> createSolver(methodType const methodType,
                                         implemType const implemType,
                                         meshType const mesh, int const order)
{
  if (methodType == SEM)
  {
    switch (implemType)
    {
      case MAKUTU:
        return make_sem_solver<IntegralType::MAKUTU>(order, mesh);
        // case SHIVA:
        //   return make_sem_solver<IntegralType::SHIVA>(order, mesh);
    }
  }

  // Add DG or other methods as needed
  throw std::runtime_error("Unsupported solver configuration");
}
}  // namespace SolverFactory
