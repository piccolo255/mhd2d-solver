#ifndef ENUMS_HPP_INCLUDED
#define ENUMS_HPP_INCLUDED

// C++ headers
#include <string>
#include <algorithm>

// Function return status
enum class ReturnStatus
   { OK
   , NoChange
   , Updated

   , ErrorNotImplemented
   , ErrorTimeUnderflow
   , ErrorNotConverged
   , ErrorNegativeDensity
   , ErrorNegativePressure

   , ErrorFileNotFound
   , ErrorFailedToParse
   , ErrorWrongParameter
};

// Output skipping modes
enum class WriteSkipMode
   { Undefined
   , Step
   , Time
};

// Output file modes
enum class WriteFileCountMode
   { Undefined
   , Single
   , Multiple
};

// Time stepping modes
enum class TimeStepMode
   { Undefined
   , Constant
   , Variable
};

// Time stepping methods
enum class TimeStepMethod
   { Undefined
   , Euler
   , RungeKutta3_TVD
};

// Spatial integration methods
enum class IntegrationMethod
   { Undefined
   , CentralFD
   , ENO_Roe
   , ENO_LF
};

// Boundary conditions
enum class BoundaryCondition
   { Undefined
   , Periodic
   , Dirichlet
   , Neumann
   , Open
};

enum class DivBCorrectionMethod
   { Undefined
   , SOR
};

template <typename T>
T fromString
   ( std::string name );

template <typename T>
std::string toString
   ( T member );

#endif // ENUMS_HPP_INCLUDED
