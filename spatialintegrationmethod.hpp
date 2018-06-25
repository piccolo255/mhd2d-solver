#ifndef SPATIALINTEGRATIONMETHOD_H
#define SPATIALINTEGRATIONMETHOD_H

#include "mhd2d.hpp"

class SpatialIntegrationMethod
{
   public:
      explicit SpatialIntegrationMethod( int          nx
                                       , int          ny
                                       , int          bufferWidth
                                       , double       dx
                                       , double       dy
                                       , t_boundary   boundary
                                       , double       gamma
                                       );
      virtual ~SpatialIntegrationMethod();

      virtual t_status integrate( t_matrices    U
                                , t_matrices    UL
                                , double       &ideal_dt
                                ) = 0;

   protected:
      const int         nxProblem;
      const int         nyProblem;
      const int         bufferWidth;
      const double      dx;
      const double      dy;
      const t_boundary  boundary;
      const double      gamma;

      const int nxFirst;
      const int nyFirst;
      const int nxLast;
      const int nyLast;
      const int nxTotal;
      const int nyTotal;

      virtual t_status applyBoundaryConditions( t_matrices U
                                          );

   private:
};

#endif // SPATIALINTEGRATIONMETHOD_H
