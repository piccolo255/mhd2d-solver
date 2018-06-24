#ifndef SPATIALINTEGRATIONMETHOD_H
#define SPATIALINTEGRATIONMETHOD_H

#include "mhd2d.hpp"

class SpatialIntegrationMethod
{
   public:
      explicit SpatialIntegrationMethod( int          nx
                                       , int          ny
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

   private:
      int         nx;
      int         ny;
      double      dx;
      double      dy;
      t_boundary  boundary;
      double      gamma;
};

#endif // SPATIALINTEGRATIONMETHOD_H
