#include "spatialmethodenoroe.hpp"

#include <stdexcept>

SpatialMethodEnoRoe::SpatialMethodEnoRoe
   ( size_t       nx
   , size_t       ny
   , size_t       bufferWidth
   , double       dx
   , double       dy
   , t_boundary   boundary
   , double       gamma
)
   : SpatialMethodEno{ nx, ny, bufferWidth, dx, dy, boundary, gamma }
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SpatialMethodEnoRoe::~SpatialMethodEnoRoe
   (
){
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoRoe::integrate
   ( t_matrices   U
   , t_matrices   UL
   , double      &dtIdeal
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO Roe method not yet implemented" };
}
