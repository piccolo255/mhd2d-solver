#ifndef SPATIALMETHODCENTRALFD_H
#define SPATIALMETHODCENTRALFD_H

#include "spatialintegrationmethod.hpp"

class SpatialMethodCentralFD2 : public SpatialIntegrationMethod
{
   public:
      explicit SpatialMethodCentralFD2( int          nx
                                     , int          ny
                                     , int          bufferWidth
                                     , double       dx
                                     , double       dy
                                     , t_boundary   boundary
                                     , double       gamma
                                     );
      virtual ~SpatialMethodCentralFD2();

      t_status integrate( t_matrices    U
                        , t_matrices    UL
                        , double       &ideal_dt
                        ) override;

   protected:

   private:
      const int requiredBufferWidth = 1;

      t_matrices F;
      t_matrices G;
};

#endif // SPATIALMETHODCENTRALFD_H
