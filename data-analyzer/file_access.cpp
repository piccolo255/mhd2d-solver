/******************************************************************************
 *                                                                            *
 *                         Reading and writing files                          *
 *                                                                            *
 ******************************************************************************/

// Main header
#include "file_access.hpp"

// Boost headers
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

// Local headers
#include "data-analyzer.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// local functions

// Read value of a key, show warning and default value if the key is missing
template <class T>
T readEntry
   ( boost::property_tree::ptree pt
   , std::string section
   , std::string name
   , T           defaultValue );

// Read settings from property tree: Plasma sheet problem
void readProblemPlasmaSheet
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void openFile
   ( t_output &output
){
   if( output.filename == "stdout" ){
      output.file = stdout;
   } else {
      if( output.binary ){
         output.file = fopen( output.filename.c_str(), "wb" );
      } else {
         output.file = fopen( output.filename.c_str(), "w" );
      }
      if( !output.file ){
         criticalError( ReturnStatus::ErrorFileNotFound, std::string{}
                      + "openFile: Unable to open file '"
                      + output.filename + "' for writing." );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void closeFile
   ( t_output &output
){
   if( output.file != stdout ){
      fclose( output.file );
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Read value of a key, show warning and default value if the key is missing */
template <class T>
T readEntry
   ( boost::property_tree::ptree pt
   , std::string section
   , std::string name
   , T           defaultValue
){
   T value;

   try {
      // get value
      value = pt.get<T>( section+"."+name );
   } catch( boost::property_tree::ptree_error &err ) {
      // show warning if key is missing
      std::cerr << "WARNING: readEntry: Key \"" << name << "\" in section [" << section << "] not found.\n"
                << "         Assuming " << name << " = " << defaultValue << "\n";
      value = defaultValue;
   }

   return value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inputData
   ( const std::string filename
   , t_output &output_grid
   , t_output &output_non_grid
   , t_params &params
   , t_data   &data
){
   namespace nptree = ::boost::property_tree;
   nptree::ptree pt;
   std::string tempstr;

   // Parse configuration file
   try {
      nptree::read_ini( filename, pt );
   } catch( std::exception &e ) {
      criticalError( ReturnStatus::ErrorFailedToParse, std::string{}
                   + "inputData: " + e.what() );
   }

   // Grid size
   params.nx = readEntry<int>( pt, "problem", "Nx", 16 );
   params.ny = readEntry<int>( pt, "problem", "Ny", 16 );

   // Prepare data storage
   data.U = createMatrices( PRB_DIM, NX, NY );
   data.u = createMatrices( VEL_DIM, NX, NY );
   data.p = createMatrix( NX, NY );

   // Time stepping mode
   // - keep a constant dt, or
   // - automatically change dt to satisfy CFL condition
   tempstr = readEntry<std::string>( pt, "time", "mode", "constant" );
   params.time_mode = fromString<TimeStepMode>( tempstr );
   switch( params.time_mode ){
   case TimeStepMode::Undefined:
      criticalError( ReturnStatus::ErrorWrongParameter, std::string{}
                   + "inputData: in section [time], key \"mode\":\n"
                   + "Unknown time stepping mode: " + tempstr );
      break;
   case TimeStepMode::Constant:
      params.steps  = readEntry<int>   ( pt, "time", "steps", 8     );
      params.dt_max = readEntry<double>( pt, "time", "dt",    0.125 );
      break;
   case TimeStepMode::Variable:
      params.dt_max     = readEntry<double>( pt, "time", "dt_max",     1.0 );
      params.dt_min     = readEntry<double>( pt, "time", "dt_min",     params.dt_max/16.0 );
      params.t_max      = readEntry<double>( pt, "time", "t_max",      1.0 );
      break;
   }

   // Initialize time
   data.dt        = params.dt_max;
   data.t_current = 0.0;

   // Set up the problem
   tempstr = readEntry<std::string>( pt, "problem", "type", "shock tube" );
   if( tempstr == "plasma sheet" ){
      readProblemPlasmaSheet( pt, params, data );
   } else {
      criticalError( ReturnStatus::ErrorWrongParameter, std::string{}
                   + "inputData: in section [problem], key \"type\":\n"
                   + "Unknown problem type: " + tempstr );
   }

   // Data output file
   output_grid.filename     = readEntry<std::string>( pt, "output grid", "datafile", "%f.g.dat" );
   output_grid.natural      = readEntry<bool>( pt, "output grid", "natural",      true );
   output_grid.conservation = readEntry<bool>( pt, "output grid", "conservation", true );

   // Binary or text file
   output_grid.binary = readEntry<bool>( pt, "output grid", "binary", false );
   if( output_grid.binary ){
      output_grid.binary_double = readEntry<bool>( pt, "output grid", "binary use double", false );
   }

   // Single file or separate files
   output_grid.single_file = readEntry<bool>( pt, "output grid", "single file", true );

   // How often to output data, can be given as time or as steps
   tempstr = readEntry<std::string>( pt, "output grid", "mode", "step" );
   output_grid.skip_mode = fromString<WriteSkipMode>( tempstr );
   switch( output_grid.skip_mode ){
   case WriteSkipMode::Undefined:
      criticalError( ReturnStatus::ErrorWrongParameter, std::string{}
                   + "inputData: in section [output grid], key \"mode\":\n"
                   + "Unknown output mode: " + tempstr );
      break;
   case WriteSkipMode::Step:
      output_grid.skip_steps = readEntry<int>( pt, "output grid", "skip steps", 1 );
      break;
   case WriteSkipMode::Time:
      output_grid.skip_t = readEntry<double>( pt, "output grid", "skip t", params.dt_max );
      break;
   }

   // Overview data output
   output_non_grid.filename = readEntry<std::string>( pt, "output non grid", "datafile", "%f.ng.dat" );

   tempstr = readEntry<std::string>( pt, "output non grid", "mode", "step" );
   output_non_grid.skip_mode = fromString<WriteSkipMode>( tempstr );
   switch( output_non_grid.skip_mode ){
   case WriteSkipMode::Undefined:
      criticalError( ReturnStatus::ErrorWrongParameter, std::string{}
                   + "inputData: in section [output non grid], key \"mode\":\n"
                   + "Unknown output mode: " + tempstr );
      break;
   case WriteSkipMode::Step:
      output_non_grid.skip_steps = readEntry<int>( pt, "output non grid", "skip steps", 1 );
      break;
   case WriteSkipMode::Time:
      output_non_grid.skip_t = readEntry<double>( pt, "output non grid", "skip t", params.dt_max );
      break;
   }

   // Get filename without the extension
   std::string mainname = filename.substr( 0, filename.find_last_of( "." ) );
   // Replace masks for output filenames
   boost::replace_all( output_grid.filename, "%f", filename );
   boost::replace_all( output_grid.filename, "%m", mainname );
   boost::replace_all( output_non_grid.filename, "%f", filename );
   boost::replace_all( output_non_grid.filename, "%m", mainname );
   if( output_grid.single_file ){
      // remove fields that are only used for separate output
      boost::replace_all( output_grid.filename, "%rs", "" );
      boost::replace_all( output_grid.filename, "%ri", "" );
      boost::replace_all( output_grid.filename, "%rt", "" );
   }

   // Read formatting strings
   if( output_grid.filename.find("%rs") != std::string::npos ){
      output_grid.format_step = readEntry<std::string>( pt, "output grid", "rs format", "%08d" );
   }
   if( output_grid.filename.find("%ri") != std::string::npos ){
      output_grid.format_index = readEntry<std::string>( pt, "output grid", "ri format", "%08d" );
   }
   if( output_grid.filename.find("%rt") != std::string::npos ){
      output_grid.format_simtime = readEntry<std::string>( pt, "output grid", "rt format", "%015.8f" );
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readProblemPlasmaSheet
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data
){
   // Domain size and position
   double Lx = readEntry<double>( pt, "problem", "Lx", 2.0 ); params.dx = Lx/params.nx;
   double Ly = readEntry<double>( pt, "problem", "Ly", 2.0 ); params.dy = Ly/params.ny;
   params.start_x = readEntry<double>( pt, "problem", "start_x", 0.0 );
   params.start_y = readEntry<double>( pt, "problem", "start_y", 0.0 );
}
