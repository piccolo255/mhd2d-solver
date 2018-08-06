#include "spatialmethodenolf.hpp"

#include <stdexcept>

SpatialMethodEnoLF::SpatialMethodEnoLF
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
SpatialMethodEnoLF::~SpatialMethodEnoLF
   (
){
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoLF::integrate
   ( t_matrices   U
   , t_matrices   UL
   , double      &dtIdeal
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO-LF is not yet implemented" };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoLF::getNumericalFluxF
   ( double    U1[PRB_DIM]
   , double    U2[PRB_DIM]
   , double    F[PRB_DIM][8]
   , double    alpha[PRB_DIM]
   , double    F_[PRB_DIM]
   , double   &maxWaveSpeed
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO-LF is not yet implemented" };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoLF::getNumericalFluxG
   ( double    U1[PRB_DIM]
   , double    U2[PRB_DIM]
   , double    G[PRB_DIM][8]
   , double    alpha[PRB_DIM]
   , double    G_[PRB_DIM]
   , double   &maxWaveSpeed
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO-LF is not yet implemented" };
}
