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
   VUF = createMatrices( PRB_DIM, nxTotal, nyTotal );
   VUG = createMatrices( PRB_DIM, nxTotal, nyTotal );
   F_  = createMatrices( PRB_DIM, nxTotal, nyTotal );
   G_  = createMatrices( PRB_DIM, nxTotal, nyTotal );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SpatialMethodEnoRoe::~SpatialMethodEnoRoe
   (
){
   freeMatrices( VUF );
   freeMatrices( VUG );
   freeMatrices( F_ );
   freeMatrices( G_ );

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoRoe::integrate
   ( t_matrices   U
   , t_matrices   UL
   , double      &dtIdeal
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO Roe method not yet implemented" };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoRoe::getNumericalFluxF
   ( double    U1[8]
   , double    U2[8]
   , double    F[8][8]
   , double    VUF[8][7]
   , double    F_[8]
   , double   &maxWaveSpeed
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO Roe method not yet implemented" };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoRoe::getNumericalFluxG
   ( double    U1[8]
   , double    U2[8]
   , double    G[8][8]
   , double    VUG[8][7]
   , double    G_[8]
   , double   &maxWaveSpeed
){
   return { true, ReturnStatus::ErrorNotImplemented, "ENO Roe method not yet implemented" };
}
