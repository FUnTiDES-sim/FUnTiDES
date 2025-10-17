#include "solver_factory.h"

#include <model_struct.h>
#include <model_unstruct.h>

#include "sem_solver_acoustic.h"
#include "sem_solver_elastic.h"
namespace SolverFactory
{

template <typename FUNC>
std::unique_ptr<SEMSolverBase> orderDispatch(int const order, FUNC&& func)
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
static std::unique_ptr<SEMSolverBase> make_sem_solver(int order, meshType mesh, modelLocationType modelLocation, physicType physic)
{
  bool isModelOnNodes = (modelLocation == OnNodes);

  switch (mesh)
  {
    case Struct:
      switch(physic)
      {
        case Acoustic:
          return orderDispatch(
            order, [isModelOnNodes](auto orderIC) -> std::unique_ptr<SEMSolverBase> {
              constexpr int ORDER = decltype(orderIC)::value;
              using SelectedIntegral =
                  typename IntegralTypeSelector<ORDER, ImplTag>::type;
              using MeshT = model::ModelStruct<float, int, ORDER>;
              if (isModelOnNodes)
              {
                return std::make_unique<
                    SEMsolverAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
              }
              else
              {
                return std::make_unique<
                    SEMsolverAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
              }
            });
        case Elastic:
          return orderDispatch(
            order, [isModelOnNodes](auto orderIC) -> std::unique_ptr<SEMSolverBase> {
              constexpr int ORDER = decltype(orderIC)::value;
              using SelectedIntegral =
                  typename IntegralTypeSelector<ORDER, ImplTag>::type;
              using MeshT = model::ModelStruct<float, int, ORDER>;
              if (isModelOnNodes)
              {
                return std::make_unique<
                    SEMsolverElastic<ORDER, SelectedIntegral, MeshT, true>>();
              }
              else
              {
                return std::make_unique<
                    SEMsolverElastic<ORDER, SelectedIntegral, MeshT, false>>();
              }
            });
      }
    case Unstruct:
      switch (physic)
      {
        case Acoustic:
          return orderDispatch(
            order, [isModelOnNodes](auto orderIC) -> std::unique_ptr<SEMSolverBase> {
              constexpr int ORDER = decltype(orderIC)::value;
              using SelectedIntegral =
                  typename IntegralTypeSelector<ORDER, ImplTag>::type;
              using MeshT = model::ModelUnstruct<float, int>;
              if (isModelOnNodes)
              {
                return std::make_unique<
                    SEMsolverAcoustic<ORDER, SelectedIntegral, MeshT, true>>();
              }
              else
              {
                return std::make_unique<
                    SEMsolverAcoustic<ORDER, SelectedIntegral, MeshT, false>>();
              }
            });
        case Elastic:
          return orderDispatch(
            order, [isModelOnNodes](auto orderIC) -> std::unique_ptr<SEMSolverBase> {
              constexpr int ORDER = decltype(orderIC)::value;
              using SelectedIntegral =
                  typename IntegralTypeSelector<ORDER, ImplTag>::type;
              using MeshT = model::ModelUnstruct<float, int>;
              if (isModelOnNodes)
              {
                return std::make_unique<
                    SEMsolverElastic<ORDER, SelectedIntegral, MeshT, true>>();
              }
              else
              {
                return std::make_unique<
                    SEMsolverElastic<ORDER, SelectedIntegral, MeshT, false>>();
              }
            });

      }
      default:
      break;


  }
  throw std::runtime_error("Unknown mesh type");
}

std::unique_ptr<SEMSolverBase> createSolver(methodType const methodType,
                                         implemType const implemType,
                                         meshType const mesh, 
                                         modelLocationType const modelLocation,
                                         physicType const physicType,
                                         int const order)
{
  if (methodType == SEM)
  {
    switch (implemType)
    {
      case MAKUTU:
        return make_sem_solver<IntegralType::MAKUTU>(order, mesh,
                                                     modelLocation, physicType);
        // case SHIVA:
        //   return make_sem_solver<IntegralType::SHIVA>(order, mesh,
        //   modelLocation);
    }
  }

  // Add DG or other methods as needed
  throw std::runtime_error(
      "Unsupported solver configuration: methodType=" + to_string(methodType) +
      ", implemType=" + to_string(implemType) + ", physicType=" + to_string(physicType));
}
}  // namespace SolverFactory
