#ifndef TIMEINTEGRATIONEULER_HPP
#define TIMEINTEGRATIONEULER_HPP

#include "timeintegrationmethod.hpp"

class TimeIntegrationEuler : public TimeIntegrationMethod
{
   public:
      explicit TimeIntegrationEuler( size_t nx
                                   , size_t ny
                                   , size_t bufferWidth
                                   , double dtMin
                                   , double dtMax
                                   , double cflNumber
                                   , std::unique_ptr<SpatialIntegrationMethod> method
                                   );
      virtual ~TimeIntegrationEuler();

      t_status step( t_matrices     U
                   , t_matrices     cx
                   , t_matrices     cy
                   , t_matrices      LUx
                   , t_matrices      LUy
                   , borderVectors  borderFlux
                   , double        &dtCurrent
                   ) override;

   protected:

   private:
      t_matrices UL;
};

#endif // TIMEINTEGRATIONEULER_HPP
