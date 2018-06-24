#include "spatialintegrationmethod.hpp"

SpatialIntegrationMethod::SpatialIntegrationMethod
   ( int          nx
   , int          ny
   , double       dx
   , double       dy
   , t_boundary   boundary
   , double       gamma
)
   : nx( nx )
   , ny( ny )
   , dx( dx )
   , dy( dy )
   , boundary( boundary )
   , gamma( gamma )
{
   //ctor
}

SpatialIntegrationMethod::~SpatialIntegrationMethod
   (
){
   //dtor
}
