#include "spatialintegrationmethod.hpp"

SpatialIntegrationMethod::SpatialIntegrationMethod
   ( int          nx
   , int          ny
   , int          bufferWidth
   , double       dx
   , double       dy
   , t_boundary   boundary
   , double       gamma
)
   : nxProblem{ nx } , nyProblem{ ny }
   , bufferWidth{ bufferWidth }
   , dx{ dx }, dy{ dy }
   , boundary{ boundary }
   , gamma{ gamma }
   , nxFirst{ bufferWidth      }, nyFirst{ bufferWidth      }
   , nxLast { nx+bufferWidth   }, nyLast { ny+bufferWidth   }
   , nxTotal{ nx+2*bufferWidth }, nyTotal{ ny+2*bufferWidth }
{
   //ctor
}

SpatialIntegrationMethod::~SpatialIntegrationMethod
   (
){
   //dtor
}
