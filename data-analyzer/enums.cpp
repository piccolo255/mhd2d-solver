/******************************************************************************
 *                                                                            *
 *               Enums and enum <> string conversion functions                *
 *                                                                            *
 ******************************************************************************/

// Main header
#include "enums.hpp"

// C++ headers
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
WriteSkipMode fromString<WriteSkipMode>
   ( std::string name
){
   WriteSkipMode member;
   std::transform( name.begin(), name.end(), name.begin(), ::tolower );

   if( false ){ // just to align the ifs below
   } else if( name == "step" ){
      member = WriteSkipMode::Step;
   } else if( name == "time" ){
      member = WriteSkipMode::Time;
   } else {
      member = WriteSkipMode::Undefined;
   }

   return member;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
std::string toString<WriteSkipMode>
   ( WriteSkipMode member
){
   std::string name;

   switch( member ){
   case WriteSkipMode::Undefined:
      name = "undefined";
      break;
   case WriteSkipMode::Step:
      name = "step";
      break;
   case WriteSkipMode::Time:
      name = "time";
      break;
   }

   return name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
WriteFileCountMode fromString<WriteFileCountMode>
   ( std::string name
){
   WriteFileCountMode member;
   std::transform( name.begin(), name.end(), name.begin(), ::tolower );

   if( false ){ // just to align the ifs below
   } else if( name == "single" ){
      member = WriteFileCountMode::Single;
   } else if( name == "multiple" ){
      member = WriteFileCountMode::Multiple;
   } else {
      member = WriteFileCountMode::Undefined;
   }

   return member;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
std::string toString<WriteFileCountMode>
   ( WriteFileCountMode member
){
   std::string name;

   switch( member ){
   case WriteFileCountMode::Undefined:
      name = "undefined";
      break;
   case WriteFileCountMode::Single:
      name = "single";
      break;
   case WriteFileCountMode::Multiple:
      name = "multiple";
      break;
   }

   return name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
TimeStepMode fromString<TimeStepMode>
   ( std::string name
){
   TimeStepMode member;
   std::transform( name.begin(), name.end(), name.begin(), ::tolower );

   if( false ){ // just to align the ifs below
   } else if( name == "constant" ){
      member = TimeStepMode::Constant;
   } else if( name == "variable" ){
      member = TimeStepMode::Variable;
   } else {
      member = TimeStepMode::Undefined;
   }

   return member;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
std::string toString<TimeStepMode>
   ( TimeStepMode member
){
   std::string name;

   switch( member ){
   case TimeStepMode::Undefined:
      name = "undefined";
      break;
   case TimeStepMode::Constant:
      name = "constant";
      break;
   case TimeStepMode::Variable:
      name = "variable";
      break;
   }

   return name;
}
