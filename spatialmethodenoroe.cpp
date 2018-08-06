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
   // for processing function return values
   auto status = t_status{};
   // track maximum wave speeds (i.e., eigenvalues)
   double maxWaveSpeedX = 0.0;
   double maxWaveSpeedY = 0.0;

   // temporary storage for function calls
   const int TN      = 8;
   const int TNFIRST = 2;
   //const int TNLAST  = 6;
   double tU1[PRB_DIM], tU2[PRB_DIM];
   double tF[PRB_DIM][TN], tVUF[PRB_DIM][TN-1], tF_[PRB_DIM];
   double tMaxWaveSpeed;

   // first, boundary conditions
   status = applyBoundaryConditions( U );
   if( status.isError ){
      status.message += "\n! SpatialMethodEnoRoe::integrate";
      return status;
   }

   // data is ready, calculate fluxes
   updateFluxes( U );

   // undivided differences V_UF, V_UG
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      for( auto i = nxFirst-2; i < nxLast+1; i++ ){
         for( auto j = nyFirst-2; j < nyLast+1; j++ ){
            VUF[k][i][j] = F[k][i+1][j] - F[k][i][j];
            VUG[k][i][j] = G[k][i][j+1] - G[k][i][j];
         }
      }
   }

   /******************************
    *           F flux           *
    ******************************/
   for( auto i = nxFirst-1; i < nxLast; i++ ){
      for( auto j = nyFirst; j < nyLast; j++ ){
         // prepare data
         for( auto k = size_t{0}; k < PRB_DIM; k++ ){
            if( boundary.left == BoundaryCondition::Open && i < nxFirst ){
               tU1[k] = U[k][i+1][j];       // left point
            } else {
               tU1[k] = U[k][i][j];         // left point
            }
            if( boundary.right == BoundaryCondition::Open && i > nxLast-2 ){
               tU2[k] = U[k][i][j];         // right point
            } else {
               tU2[k] = U[k][i+1][j];       // right point
            }
            for( auto l = size_t{0}; l < 4; l++ ){
               tF[k][TNFIRST+l] = F[k][i-1+l][j]; // physical flux
            }
            for( auto l = size_t{0}; l < 3; l++ ){
               tVUF[k][TNFIRST+l] = VUF[k][i-1+l][j]; // undivided differences
            }
         }

         // find numerical flux
         status = getNumericalFluxF( tU1, tU2, tF, tVUF, tF_, tMaxWaveSpeed );
         if( status.isError ){
            status.message += "\n! SpatialMethodEnoRoe::integrate: F flux";
            return status;
         }

         //process results
         for( auto k = size_t{0}; k < PRB_DIM; k++ ){
            F_[k][i][j] = tF_[k];
         }
         if( tMaxWaveSpeed > maxWaveSpeedX )
            maxWaveSpeedX = tMaxWaveSpeed;
      }
   }

   /******************************
    *           G flux           *
    ******************************/
   for( auto i = nxFirst; i < nxLast; i++ ){
      for( auto j = nyFirst-1; j < nyLast; j++ ){
         // prepare data
         for( auto k = size_t{0}; k < PRB_DIM; k++ ){
            if( boundary.bottom == BoundaryCondition::Open && j < nyFirst ){
               tU1[k] = U[k][i][j+1];       // down point
            } else {
               tU1[k] = U[k][i][j];         // down point
            }
            if( boundary.top == BoundaryCondition::Open && j > nyLast-2 ){
               tU2[k] = U[k][i][j];         // up point
            } else {
               tU2[k] = U[k][i][j+1];       // up point
            }
            for( auto l = size_t{0}; l < 4; l++ ){
               tF[k][TNFIRST+l] = G[k][i][j-1+l]; // physical flux
            }
            for( auto l = size_t{0}; l < 3; l++ ){
               tVUF[k][TNFIRST+l] = VUG[k][i][j-1+l]; // undivided differences
            }
         }

         // find numerical flux
         status = getNumericalFluxG( tU1, tU2, tF, tVUF, tF_, tMaxWaveSpeed );
         if( status.isError ){
            status.message += "\n! SpatialMethodEnoRoe::integrate: G flux";
            return status;
         }

         //process results
         for( auto k = size_t{0}; k < PRB_DIM; k++ ){
            G_[k][i][j] = tF_[k];
         }
         if( tMaxWaveSpeed > maxWaveSpeedY )
            maxWaveSpeedY = tMaxWaveSpeed;
      }
   }

   #ifdef DEBUG_MAX_VELOCITY
      OUT << "*** DEBUG: maxWaveSpeedX = " << maxWaveSpeedX
          << ", maxWaveSpeedY = " << maxWaveSpeedY << "\n";
   #endif // DEBUG_MAX_VELOCITY

   // update local dt from max wave speeds
   dtIdeal = std::min( dx/maxWaveSpeedX, dy/maxWaveSpeedY );

   // use numerical fluxes to calculate dU/dt
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      for( auto i = nxFirst; i < nxLast; i++ ){
         for( auto j = nyFirst; j < nyLast; j++ ){
            UL[k][i][j] = -(1.0/dx)*( F_[k][i][j] - F_[k][i-1][j] )
                          -(1.0/dy)*( G_[k][i][j] - G_[k][i][j-1] );
         }
      }
   }

   return status;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoRoe::getNumericalFluxF
   ( double    U1[PRB_DIM]
   , double    U2[PRB_DIM]
   , double    F[PRB_DIM][8]
   , double    VUF[PRB_DIM][7]
   , double    F_[PRB_DIM]
   , double   &maxWaveSpeed
){
   const auto TNFIRST = size_t{2};
   // current index
   const auto i = TNFIRST+1;

   // eigenvalues & eigenvectors @ u_(i+1/2)
   double a[PRB_DIM], lv[PRB_DIM][PRB_DIM], rv[PRB_DIM][PRB_DIM];
   // characteristics and differences
   double w[PRB_DIM][8], Vw[PRB_DIM][7];
   // fluxes on minus edge, plus edge, and final
   double wm[PRB_DIM], wp[PRB_DIM], fw[PRB_DIM];
   // maximum wave speed (i.e., eigenvalue)
   maxWaveSpeed = 0.0;

   // eigensystem
   auto status = getEigensF( U1, U2, a, lv, rv );
   if( status.status != ReturnStatus::OK ){
      status.message += "\n! SpatialMethodEnoRoe::getNumericalFluxF";
   }

   // update maximum wave speed
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      if( auto aabs = fabs(a[k]) > maxWaveSpeed )
         maxWaveSpeed = aabs;
   }

   // local characteristics and undivided differences
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      for( auto l = i-1; l <= i+2; l++ ){
         w[k][l] = 0.0;
         for( auto m = size_t{0}; m < PRB_DIM; m++ )
            w[k][l] += lv[k][m]*F[m][l];
      }
   }
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      for( auto l = i-1; l <= i+1; l++ ){
         Vw[k][l] = 0.0;
         for( auto m = size_t{0}; m < PRB_DIM; m++ )
            Vw[k][l] += lv[k][m]*VUF[m][l];
      }
   }

   // TODO: set up boundaries for open boundary conditions

   // minus edge
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      if( ABS(Vw[k][i-1]) < ABS(Vw[k][i]) ) // left (stencil -3/2, -1/2, +1/2)
         wm[k] = -(1.0/2.0)*w[k][i-1] + (3.0/2.0)*w[k][i];
      else // right (stencil -1/2, +1/2, +3/2)
         wm[k] = (1.0/2.0)*w[k][i] + (1.0/2.0)*w[k][i+1];
   }

   // plus edge
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      if( ABS(Vw[k][i]) < ABS(Vw[k][i+1]) ) // left (stencil -3/2, -1/2, +1/2)
         wp[k] = (1.0/2.0)*w[k][i] + (1.0/2.0)*w[k][i+1];
      else // right (stencil -1/2, +1/2, +3/2)
         wp[k] = (3.0/2.0)*w[k][i+1] - (1.0/2.0)*w[k][i+2];
   }

   // fluxes
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      if( a[k] >= 0 && wm[k] < 1.0e10 )
         fw[k] = wm[k];
      else if( a[k] <= 0 && wp[k] < 1.0e10 )
         fw[k] = wp[k];
      else
         fw[k] = 0.0;
   }

   // return to physical space fluxes
   for( auto k = size_t{0}; k < PRB_DIM; k++ ){
      F_[k] = 0.0;
      for( auto l = size_t{0}; l < PRB_DIM; l++ )
         F_[k] += rv[l][k]*fw[l];
   }

   return status;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialMethodEnoRoe::getNumericalFluxG
   ( double    U1[PRB_DIM]
   , double    U2[PRB_DIM]
   , double    G[PRB_DIM][8]
   , double    VUG[PRB_DIM][7]
   , double    G_[PRB_DIM]
   , double   &maxWaveSpeed
){
   // invert x and y axis
   std::swap( U1[1], U1[2] ); std::swap( U1[4], U1[5] );
   std::swap( U2[1], U2[2] ); std::swap( U2[4], U2[5] );
   for( auto i = size_t{0}; i < 8; i++ ){
      std::swap( G[1][i], G[2][i] ); std::swap( G[4][i], G[5][i] );
   }
   for( auto i = size_t{0}; i < 7; i++ ){
      std::swap( VUG[1][i], VUG[2][i] ); std::swap( VUG[4][i], VUG[5][i] );
   }

   // find flux for inverted x and y
   auto status = getNumericalFluxF( U1, U2, G, VUG, G_, maxWaveSpeed );
   if( status.status != ReturnStatus::OK ){
      status.message += "\n! SpatialMethodEnoRoe::getNumericalFluxG";
   }

   // return found flux to proper axes
   std::swap( G_[1], G_[2] ); std::swap( G_[4], G_[5] );

   return status;
}
