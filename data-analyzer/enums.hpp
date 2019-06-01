#ifndef ENUMS_HPP_INCLUDED
#define ENUMS_HPP_INCLUDED

// C++ headers
#include <string>

// Function return status
enum class ReturnStatus
   { OK
   , NoChange
   , Updated

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

template <typename T>
T fromString
   ( std::string name );

template <typename T>
std::string toString
   ( T member );

#endif // ENUMS_HPP_INCLUDED
