#include "timeintegrationrk3.hpp"

TimeIntegrationRK3::TimeIntegrationRK3
   ( size_t nx
   , size_t ny
   , size_t bufferWidth
   , double dtMin
   , double dtMax
   , double cflNumber
   , std::unique_ptr<SpatialIntegrationMethod> method
)
   : TimeIntegrationMethod{ nx, ny, bufferWidth, dtMin, dtMax, cflNumber, std::move( method ) }
{
   UL = createMatrices( PRB_DIM, nxTotal, nyTotal );
}

TimeIntegrationRK3::~TimeIntegrationRK3
   (
){
   freeMatrices( UL );
}

t_status TimeIntegrationRK3::step
   ( t_matrices   U
   , double      &dtCurrent
){
   return { true, ReturnStatus::ErrorNotImplemented, "not implemented: 3rd order Runge-Kutta time stepper" };
}
