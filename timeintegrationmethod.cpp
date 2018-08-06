#include "timeintegrationmethod.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TimeIntegrationMethod::TimeIntegrationMethod
   ( size_t nx
   , size_t ny
   , size_t bufferWidth
   , double dtMin
   , double dtMax
   , double cflNumber
   , std::unique_ptr<SpatialIntegrationMethod> method
)
   : nxProblem{ nx } , nyProblem{ ny }
   , bufferWidth{ bufferWidth }
   , nxFirst{ 0 +bufferWidth   }, nyFirst{ 0 +bufferWidth   }
   , nxLast { nx+bufferWidth   }, nyLast { ny+bufferWidth   }
   , nxTotal{ nx+bufferWidth*2 }, nyTotal{ ny+bufferWidth*2 }
   , dtMin{ dtMin }
   , dtMax{ dtMax }
   , cflNumber{ cflNumber }
   , method{ std::move( method ) }
{
   //ctor
}

TimeIntegrationMethod::~TimeIntegrationMethod
   (
){
   //dtor
}

t_status TimeIntegrationMethod::updateDt
   ( double   &dtCurrent
   , double    dtIdeal
){
   auto dt = dtCurrent;

   // resize dt to fit between 0.5*cfl and cfl
   while( dt < 0.5*cflNumber*dtIdeal && dt <= dtMax ){
      dt *= 2.0;

      #ifdef DEBUG_DT_UPDATE
         OUT << "*** DEBUG: raised dt to " << dt << LF;
      #endif // DEBUG_DT_UPDATE
   }

   while( dt > cflNumber*dtIdeal && dt >= dtMin ){
      dt /= 2.0;

      #ifdef DEBUG_DT_UPDATE
         OUT << "*** DEBUG: lowered dt to " << dt << LF;
      #endif // DEBUG_DT_UPDATE
   }

   // clamp from above at dtMax
   if( dt > dtMax ){
      dt = dtMax;
   }

   // error if < dtMin
   if( dt < dtMin ){
      return { true, ReturnStatus::ErrorTimeUnderflow
             , "smallest allowed dt does not satisfy CFL condition\n! updateDt" };
   }

   if( dt == dtCurrent ){
      #ifdef DEBUG_DT_UPDATE
         OUT << "*** DEBUG: dt unchanged, dt = " << dt << LF;
      #endif // DEBUG_DT_UPDATE

      return { false, ReturnStatus::NoChange, "" };
   } else {
      #ifdef DEBUG_DT_UPDATE
         OUT << "*** DEBUG: dt updated, dt = " << dt << LF;
      #endif // DEBUG_DT_UPDATE
      dtCurrent = dt;

      return { false, ReturnStatus::Updated, "" };
   }
}
