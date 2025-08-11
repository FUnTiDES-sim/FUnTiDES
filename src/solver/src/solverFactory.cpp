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

std::unique_ptr<SolverBase> createSolver( methodType const methodType,
                                          implemType const implemType,
                                          int const order )
{
    if (methodType == SEM)
    {
        if (implemType == CLASSIC)
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }
        else
        if (implemType == OPTIM)
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::OPTIM >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }
        else if (implemType == GEOS)
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::GEOS >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }
        else if (implemType == SHIVA)
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, IntegralType::SHIVA >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }

    }
    // Add more physics types as needed

    throw std::runtime_error("Unsupported solver configuration");

}
} // namespace SolverFactory
