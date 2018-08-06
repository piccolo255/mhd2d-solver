#ifndef TIMEINTEGRATIONRK3_HPP
#define TIMEINTEGRATIONRK3_HPP

#include "timeintegrationmethod.hpp"

class TimeIntegrationRK3 : public TimeIntegrationMethod
{
   public:
      explicit TimeIntegrationRK3( size_t nx
                                 , size_t ny
                                 , size_t bufferWidth
                                 , double dtMin
                                 , double dtMax
                                 , double cflNumber
                                 , std::unique_ptr<SpatialIntegrationMethod> method
                                 );
      virtual ~TimeIntegrationRK3();

      t_status step( t_matrices   U
                   , double      &dtCurrent
                   ) override;

   protected:

   private:
      t_matrices U1;
      t_matrices U2;
      t_matrices UL;
};

#endif // TIMEINTEGRATIONRK3_HPP