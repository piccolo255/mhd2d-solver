#ifndef SPATIALMETHODENO_HPP
#define SPATIALMETHODENO_HPP

#include "spatialintegrationmethod.hpp"

class SpatialMethodEno : public SpatialIntegrationMethod
{
   public:
      explicit SpatialMethodEno( size_t      nx
                               , size_t      ny
                               , size_t      bufferWidth
                               , double      dx
                               , double      dy
                               , t_boundary  boundary
                               , double      gamma
                               );
      virtual ~SpatialMethodEno();

      t_status integrate( t_matrices    U
                        , t_matrices    UL
                        , double       &ideal_dt
                        ) override;

      static size_t requiredBufferWidth();

   protected:

   private:
      static const size_t minimumBufferWidth = 1;

      t_matrices F;
      t_matrices G;
};

#endif // SPATIALMETHODCENTRALFD_HPP
