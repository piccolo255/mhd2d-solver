#ifndef TIMEINTEGRATIONMETHOD_HPP
#define TIMEINTEGRATIONMETHOD_HPP

#include "mhd2d.hpp"
#include "spatialintegrationmethod.hpp"

class TimeIntegrationMethod
{
   public:
      explicit TimeIntegrationMethod( size_t nx
                                    , size_t ny
                                    , size_t bufferWidth
                                    , double dtMin
                                    , double dtMax
                                    , double cflNumber
                                    , std::unique_ptr<SpatialIntegrationMethod> method
                                    );
      virtual ~TimeIntegrationMethod();

      virtual t_status step( t_matrices      U
                           , borderVectors   borderFlux
                           , double         &dtCurrent
                           ) = 0;

   protected:
      const size_t   nxProblem;
      const size_t   nyProblem;
      const size_t   bufferWidth;

      const size_t   nxFirst;
      const size_t   nyFirst;
      const size_t   nxLast;
      const size_t   nyLast;
      const size_t   nxTotal;
      const size_t   nyTotal;

      double   dtMin;
      double   dtMax;
      double   cflNumber;
      bool     variableTime;

      std::unique_ptr<SpatialIntegrationMethod> method;

      virtual t_status updateDt( double  &dtCurrent
                               , double   dtIdeal
                               );

   private:
};

#endif // TIMEINTEGRATIONMETHOD_HPP
