#include "spatialintegrationmethod.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
   , nxFirst{ 0 +bufferWidth   }, nyFirst{ 0 +bufferWidth   }
   , nxLast { nx+bufferWidth   }, nyLast { ny+bufferWidth   }
   , nxTotal{ nx+bufferWidth*2 }, nyTotal{ ny+bufferWidth*2 }
{
   //ctor
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SpatialIntegrationMethod::~SpatialIntegrationMethod
   (
){
   //dtor
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
t_status SpatialIntegrationMethod::applyBoundaryConditions
   ( t_matrices U
){
   // left boundary
   switch( boundary.left ){
   case BoundaryCondition::Undefined:
      return { true, ReturnStatus::ErrorWrongParameter, "Unknown left boundary condition." };
      break;
   case BoundaryCondition::Periodic:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               U[k][nxFirst-i][j] = U[k][nxLast-i][j];
            }
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               /// TODO
               //U[k][nxFirst-i][j] = boundary_dirichlet_U[params.b_left][k][j];
               return { true, ReturnStatus::ErrorNotImplemented, "Not implemented: Dirichlet BC." };
            }
         }
      }
      break;
   case BoundaryCondition::Neumann:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               U[k][nxFirst-i][j] = U[k][nxFirst+i][j];
            }
         }
      }
      // use div B = 0 for dBx/dx and dBy/dy boundary conditions
      {
         double ratio_xy = dx/dy;
         for( int i = nxFirst; i > 1; i-- ){
            for( int j = nyFirst+1; j < nyLast-1; j++ ){
               U[4][i-1][j] = U[4][i+1][j] - ratio_xy*( U[5][i][j-1] - U[5][i][j+1] );
            }
         }
      }
      break;
   case BoundaryCondition::Open:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               U[k][nxFirst-i][j] = U[k][nxFirst+i][j];
            }
         }
      }
      break;
   }

   // right boundary
   switch( boundary.right ){
   case BoundaryCondition::Undefined:
      return { true, ReturnStatus::ErrorWrongParameter, "Unknown right boundary condition." };
      break;
   case BoundaryCondition::Periodic:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               U[k][nxLast-1+i][j] = U[k][nxFirst-1+i][j];
            }
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               /// TODO
               //U[k][nxLast-1+i][j] = params.boundary_dirichlet_U[params.b_right][k][j];
               return { true, ReturnStatus::ErrorNotImplemented, "Not implemented: Dirichlet BC." };
            }
         }
      }
      break;
   case BoundaryCondition::Neumann:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               U[k][nxLast-1+i][j] = U[k][nxLast-1-i][j];
            }
         }
      }
      // use div B = 0 for dBx/dx and dBy/dy boundary conditions
      {
         double ratio_xy = dx/dy;
         for( int i = nxLast-1; i < nxTotal-1; i++ ){
            for( int j = nyFirst+1; j < nyLast-1; j++ ){
               U[4][i+1][j] = U[4][i-1][j] - ratio_xy*( U[5][i][j+1] - U[5][i][j-1] );
            }
         }
      }
      break;
   case BoundaryCondition::Open:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = bufferWidth; i > 0; i-- ){
            for( int j = nyFirst; j < nyLast; j++ ){
               U[k][nxLast-1+i][j] = U[k][nxLast-1-i][j];
            }
         }
      }
      break;
   }

   // bottom boundary
   switch( boundary.bottom ){
   case BoundaryCondition::Undefined:
      return { true, ReturnStatus::ErrorWrongParameter, "Unknown bottom boundary condition." };
      break;
   case BoundaryCondition::Periodic:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               U[k][i][nyFirst-j] = U[k][i][nyLast-j];
            }
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               /// TODO
               //U[k][i][nyFirst-j] = params.boundary_dirichlet_U[params.b_bottom][k][i];
               return { true, ReturnStatus::ErrorNotImplemented, "Not implemented: Dirichlet BC." };
            }
         }
      }
      break;
   case BoundaryCondition::Neumann:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               U[k][i][nyFirst-j] = U[k][i][nyFirst+j];
            }
         }
      }
      // use div B = 0 for dBx/dx and dBy/dy boundary conditions
      {
         double ratio_yx = dy/dx;
         for( int i = nxFirst+1; i < nxLast-1; i++ ){
            for( int j = nyFirst; j > 1; j-- ){
               U[5][i][j-1] = U[5][i][j+1] - ratio_yx*( U[4][i-1][j] - U[4][i+1][j] );
            }
         }
      }
      break;
   case BoundaryCondition::Open:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               U[k][i][nyFirst-j] = U[k][i][nyFirst+j];
            }
         }
      }
      break;
   }

   // top boundary
   switch( boundary.top ){
   case BoundaryCondition::Undefined:
      return { true, ReturnStatus::ErrorWrongParameter, "Unknown top boundary condition." };
      break;
   case BoundaryCondition::Periodic:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               U[k][i][nyLast-1+j] = U[k][i][nyFirst-1+j];
            }
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               /// TODO
               //U[k][i][nyLast-1+j] = params.boundary_dirichlet_U[params.b_top][k][i];
               return { true, ReturnStatus::ErrorNotImplemented, "Not implemented: Dirichlet BC." };
            }
         }
      }
      break;
   case BoundaryCondition::Neumann:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               U[k][i][nyLast-1+j] = U[k][i][nyLast-1-j];
            }
         }
      }
      // use div B = 0 for dBx/dx and dBy/dy boundary conditions
      {
         double ratio_yx = dy/dx;
         for( int i = nxFirst+1; i < nxLast-1; i++ ){
            for( int j = nyLast-1; j < nyTotal-1; j++ ){
               U[5][i][j+1] = U[5][i][j-1] - ratio_yx*( U[4][i+1][j] - U[4][i-1][j] );
            }
         }
      }
      break;
   case BoundaryCondition::Open:
      for( int k = 0; k < PRB_DIM; k++ ){
         for( int i = nxFirst; i < nxLast; i++ ){
            for( int j = bufferWidth; j > 0; j-- ){
               U[k][i][nyLast-1+j] = U[k][i][nyLast-1-j];
            }
         }
      }
      break;
   }

   return { false, ReturnStatus::OK, "" };
}
