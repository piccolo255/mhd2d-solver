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

      t_status step( t_matrices   U
                   , t_vectors    borderFluxLRUD
                   , double      &dtCurrent
                   ) override;

   protected:

   private:
      t_matrices UL;
};

#endif // TIMEINTEGRATIONEULER_HPP
