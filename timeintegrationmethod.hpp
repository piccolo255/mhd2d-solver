#ifndef TIMEINTEGRATIONMETHOD_HPP
#define TIMEINTEGRATIONMETHOD_HPP

#include "mhd2d.hpp"

class TimeIntegrationMethod
{
   public:
      explicit TimeIntegrationMethod( size_t       nx
                                    , size_t       ny
                                    , size_t       bufferWidth
                                    );
      virtual ~TimeIntegrationMethod();

      virtual t_status step( t_matrices    U
                           , double       &current_dt
                           ) = 0;

   protected:
      const size_t      nxProblem;
      const size_t      nyProblem;
      const size_t      bufferWidth;

      const size_t nxFirst;
      const size_t nyFirst;
      const size_t nxLast;
      const size_t nyLast;
      const size_t nxTotal;
      const size_t nyTotal;

   private:
};

#endif // TIMEINTEGRATIONMETHOD_HPP
