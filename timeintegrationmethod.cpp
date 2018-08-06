#include "timeintegrationmethod.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TimeIntegrationMethod::TimeIntegrationMethod
   ( size_t       nx
   , size_t       ny
   , size_t       bufferWidth
)
   : nxProblem{ nx } , nyProblem{ ny }
   , bufferWidth{ bufferWidth }
   , nxFirst{ 0 +bufferWidth   }, nyFirst{ 0 +bufferWidth   }
   , nxLast { nx+bufferWidth   }, nyLast { ny+bufferWidth   }
   , nxTotal{ nx+bufferWidth*2 }, nyTotal{ ny+bufferWidth*2 }
{
   //ctor
}

TimeIntegrationMethod::~TimeIntegrationMethod
   (
){
   //dtor
}
