// Local headers
#include "data-analyzer.hpp"

// C++ headers
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main
   ( int    argc
   , char **argv
){
   /* Variable declarations */
   // output file information
   t_output output_grid;
   t_output output_non_grid;
   // parameters
   t_params params;
   // variables
   t_data data;

   /* Check command line arguments */
   if( argc < 2 ){
      std::cout << "Usage: " << argv[0] << " config_file\n";
      exit(0);
   }

   /* Initialize */
   //inputData( argv[1], output_grid, output_non_grid, params, data );

   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void criticalError
   ( ReturnStatus       error
   , const std::string  message
){
   std::cerr << "Critical error!\n"
             << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
             << message << "\n"
             << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
             << "Terminating execution.\n";

   auto errcode = int{0};

   switch( error ){
   case ReturnStatus::OK                     : errcode = 0; break;
   case ReturnStatus::NoChange               : errcode = 0; break;
   case ReturnStatus::Updated                : errcode = 0; break;

   case ReturnStatus::ErrorFileNotFound      : errcode = 200; break;
   case ReturnStatus::ErrorFailedToParse     : errcode = 201; break;
   case ReturnStatus::ErrorWrongParameter    : errcode = 202; break;
   }

   exit( errcode );
}
