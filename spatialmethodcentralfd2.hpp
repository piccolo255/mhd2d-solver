#ifndef SPATIALMETHODCENTRALFD_H
#define SPATIALMETHODCENTRALFD_H

#include "spatialintegrationmethod.hpp"

class SpatialMethodCentralFD2 : public SpatialIntegrationMethod
{
   public:
      explicit SpatialMethodCentralFD2( size_t       nx
                                      , size_t       ny
                                      , size_t       bufferWidth
                                      , double       dx
                                      , double       dy
                                      , t_boundary   boundary
                                      , double       gamma
                                      );
      virtual ~SpatialMethodCentralFD2();

      t_status integrate( t_matrices      U
                        , t_matrices      UL
                        , borderVectors   borderFlux
                        , double         &dtIdeal
                        ) override;

      static size_t requiredBufferWidth();

   protected:

   private:
      static const size_t minimumBufferWidth = 1;

      t_matrices F;
      t_matrices G;
};

#endif // SPATIALMETHODCENTRALFD_H
