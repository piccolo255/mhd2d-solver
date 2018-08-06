#include "spatialmethodeno.hpp"

#include <stdexcept>

SpatialMethodEno::SpatialMethodEno
   ( size_t       nx
   , size_t       ny
   , size_t       bufferWidth
   , double       dx
   , double       dy
   , t_boundary   boundary
   , double       gamma
)
   : SpatialIntegrationMethod{ nx, ny, bufferWidth, dx, dy, boundary, gamma }
{
   if( bufferWidth < minimumBufferWidth ){
      throw std::domain_error( "Buffer is too small for ENO;"
                               " need at least 2 extra grids" );
   }

   F = createMatrices( PRB_DIM, nxTotal, nyTotal );
   G = createMatrices( PRB_DIM, nxTotal, nyTotal );
}

SpatialMethodEno::~SpatialMethodEno
   (
){
   freeMatrices( F );
   freeMatrices( G );
}
