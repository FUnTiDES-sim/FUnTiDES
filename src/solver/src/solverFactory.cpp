#include "solverFactory.hpp"
#include "SEMsolver.hpp"

namespace SolverFactory
{

template< typename FUNC >
std::unique_ptr<SolverBase> orderDispatch( int const order,
                                           FUNC && func )
{
  if( order == 1 )
  {
    return func( integral_constant< int, 1>{} );
  }
  if( order == 2 )
  {
    return func( integral_constant< int, 2>{} );
  }
  else if( order == 3 )
  {
    return func( integral_constant< int, 3>{} );
  }
  // else if( order == 4 )
  // {
  //   return func( integral_constant< int, 4>{} );
  // }
  abort();
}

// ImplTag is one of IntegralType::CLASSIC/OPTIM/GEOS/SHIVA (compile-time)
template <auto ImplTag>
static std::unique_ptr<SolverBase>
make_sem_solver(int order, meshType mesh)
{
  switch (mesh) {
    case CARTESIAN:
      // ORDER-dependent mesh type
      return orderDispatch(order, [] (auto orderIC) -> std::unique_ptr<SolverBase> {
        constexpr int ORDER = decltype(orderIC)::value;
        using SelectedIntegral = typename IntegralTypeSelector<ORDER, ImplTag>::type;
        using MeshT = CartesianSEMmesh<float, int, ORDER>;
        return std::make_unique<SEMsolver<ORDER, SelectedIntegral, MeshT>>();
      });

    case UNSTRUCT_CARTESIAN:
    case DIVA:
      throw std::runtime_error("SEM mesh type not implemented yet");
  }

  throw std::runtime_error("Unknown mesh type");
}

std::unique_ptr<SolverBase>
createSolver(methodType const methodType,
              implemType const implemType,
              meshType   const mesh,
              int        const order)
{
  if (methodType == SEM) {
    switch (implemType) {
      case CLASSIC: return make_sem_solver<IntegralType::CLASSIC>(order, mesh);
      case GEOS:    return make_sem_solver<IntegralType::GEOS>(order, mesh);
      case OPTIM:   return make_sem_solver<IntegralType::OPTIM>(order, mesh);
      case SHIVA:   return make_sem_solver<IntegralType::SHIVA>(order, mesh);
    }
  }

  // Add DG or other methods as needed
  throw std::runtime_error("Unsupported solver configuration");
}
} // namespace SolverFactory

// std::unique_ptr<SolverBase> createSolver( methodType const methodType,
//                                           implemType const implemType,
//                                           meshType const meshType,
//                                           int const order )
// {
//     if (methodType == SEM)
//     {
//         if (implemType == CLASSIC)
//         {
//           return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
//           {
//             constexpr int ORDER = decltype(orderIC)::value;
//             using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type;
//             return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
//           });
//         }
//         else
//         if (implemType == OPTIM)
//         {
//           return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
//           {
//             constexpr int ORDER = decltype(orderIC)::value;
//             using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::OPTIM >::type;
//             return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
//           });
//         }
//         else if (implemType == GEOS)
//         {
//           return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
//           {
//             constexpr int ORDER = decltype(orderIC)::value;
//             using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::GEOS >::type;
//             return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
//           });
//         }
//         else if (implemType == SHIVA)
//         {
//           return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
//           {
//             constexpr int ORDER = decltype(orderIC)::value;
//             using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::SHIVA >::type;
//             return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
//           });
//         }

//     }
//     // Add more physics types as needed

//     throw std::runtime_error("Unsupported solver configuration");

// }
// } // namespace SolverFactory
