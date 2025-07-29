
#include "solverFactory.hpp"
#include "SEMsolver.hpp"



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

std::unique_ptr<SolverBase> createSolver( int const physicsType,
                                          int const methodType,
                                          int const order )
{
    if (physicsType == 0) 
    {
        if (methodType == 0) 
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, 0 >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }
        else 
        if (methodType == 1) 
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, 1 >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }
        else if (methodType == 2) 
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, 2 >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }
        else if (methodType == 3) 
        {
          return orderDispatch( order, []( auto orderIC ) -> std::unique_ptr<SolverBase>
          {
            constexpr int ORDER = decltype(orderIC)::value;
            using IntegralType = typename IntegralTypeSelector< ORDER, 3 >::type;
            return std::unique_ptr<SolverBase>( new SEMsolver<ORDER, IntegralType>() );
          });
        }

    }
    // Add more physics types as needed

    throw std::runtime_error("Unsupported solver configuration");

}