#include "mhd2d.hpp"
#include "file_access.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void openFile
   ( t_output &output
){
   if( output.filename == "stdout" )
      output.file = stdout;
   else {
      if( output.binary ){
         output.file = fopen( output.filename.c_str(), "wb" );
      } else {
         output.file = fopen( output.filename.c_str(), "w" );
      }
      if( !output.file ){
         ERROUT << "ERROR: openFile: Unable to open file '" << output.filename << "' for writing." << LF;
         exit( RET_ERR_FILE_NOT_FOUND );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void closeFile
   ( t_output &output
){
   if( output.file != stdout )
      fclose( output.file );
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
      ERROUT << "WARNING: readEntry: Key \"" << name << "\" in section [" << section << "] not found.\n"
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
   namespace nptree = boost::property_tree;
   nptree::ptree pt;
   std::string tempstr;

   // Parse configuration file
   try {
      nptree::read_ini( filename, pt );
   } catch( std::exception &e ) {
      ERROUT << "ERROR: inputData: " << e.what() << LF;
      exit( RET_ERR_PROPERTY_TREE );
   }

   // Global problem parameters
   params.gamma = readEntry<double>( pt, "problem", "gamma", 1.4 );

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
      ERROUT << "ERROR: inputData: in section [time], key \"mode\":\n"
             << "       Unknown time stepping mode: " << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
      break;
   case TimeStepMode::Constant:
      params.steps  = readEntry<int>   ( pt, "time", "steps", 8     );
      params.dt_max = readEntry<double>( pt, "time", "dt",    0.125 );
      break;
   case TimeStepMode::Variable:
      params.cfl_number = readEntry<double>( pt, "time", "cfl number", 0.5 );
      params.dt_max     = readEntry<double>( pt, "time", "dt_max",     1.0 );
      params.dt_min     = readEntry<double>( pt, "time", "dt_min",     params.dt_max/16.0 );
      params.t_max      = readEntry<double>( pt, "time", "t_max",      1.0 );
      break;
   }

   // Initialize time
   data.dt        = params.dt_max;
   data.t_current = 0.0;

   // DivB correction
   readParamsDivBCorrector( pt, params );

   // Set up the problem
   tempstr = readEntry<std::string>( pt, "problem", "type", "shock tube" );
   if( tempstr == "shock tube" ){
      readProblemShockTube( pt, params, data );
   } else if( tempstr == "explosion" ){
      readProblemExplosion( pt, params, data );
   } else if( tempstr == "orszag-tang vortex" ){
      readProblemOTVortex( pt, params, data );
   } else if( tempstr == "plasma sheet" ){
      readProblemPlasmaSheet( pt, params, data );
   } else {
      ERROUT << "ERROR: inputData: in section [problem], key \"type\":\n"
             << "       Unknown problem type: " << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
   }

   // Negative pressure handling
   params.break_on_neg_pressure = readEntry<bool>( pt, "problem", "halt on negative pressure", false );

   // Time integration method
   tempstr = readEntry<std::string>( pt, "problem", "time method", "rk3" );
   params.time_stepping = fromString<TimeStepMethod>( tempstr );
   if( params.time_stepping == TimeStepMethod::Undefined ){
      ERROUT << "ERROR: inputData: in section [problem], key \"time method\":\n"
             << "       Unknown time stepping method: " << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
   }

   // Space integration method
   tempstr = readEntry<std::string>( pt, "problem", "space method", "eno-roe" );
   params.scheme = fromString<IntegrationMethod>( tempstr );
   if( params.scheme == IntegrationMethod::Undefined ){
      ERROUT << "ERROR: inputData: in section [problem], key \"space method\":\n"
             << "       Unknown space integration method: " << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
   }

   // Data output file
   output_grid.filename = readEntry<std::string>( pt, "output grid", "datafile", "%f.g.dat" );
   output_grid.natural      = readEntry<bool>( pt, "output grid", "natural",      true );
   output_grid.conservation = readEntry<bool>( pt, "output grid", "conservation", true );

   // Binary or text file
   output_grid.binary = readEntry<bool>( pt, "output grid", "binary", false );
   std::string binary_hint_file;
   if( output_grid.binary ){
      output_grid.binary_double = readEntry<bool>( pt, "output grid", "binary use double", false );
      binary_hint_file = readEntry<std::string>( pt, "output grid", "binary hint", "%f.hint" );
   } else {
      binary_hint_file = "dummy";
   }

   // How often to output data, can be given as time or as steps
   tempstr = readEntry<std::string>( pt, "output grid", "mode", "step" );
   output_grid.skip_mode = fromString<WriteSkipMode>( tempstr );
   switch( output_grid.skip_mode ){
   case WriteSkipMode::Undefined:
      ERROUT << "ERROR: inputData: in section [output grid], key \"mode\":\n"
             << "       Unknown output mode: " << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
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
      ERROUT << "ERROR: inputData: in section [output non grid], key \"mode\":\n"
             << "       Unknown output mode: " << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
      break;
   case WriteSkipMode::Step:
      output_non_grid.skip_steps = readEntry<int>( pt, "output non grid", "skip steps", 1 );
      break;
   case WriteSkipMode::Time:
      output_non_grid.skip_t = readEntry<double>( pt, "output non grid", "skip t", params.dt_max );
      break;
   }

   // Construct time string
   std::time_t tt = std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() );
   struct std::tm * ptm = std::localtime(&tt);
   char buffer[80];
   strftime( buffer, 80, "%Y-%m-%d-%H-%M-%S", ptm );
   std::string timestr( buffer );
   // Get filename without the extension
   std::string mainname = filename.substr( 0, filename.find_last_of( "." ) );
   // Replace masks for output filenames
   boost::replace_all( output_grid.filename, "%f", filename );
   boost::replace_all( output_grid.filename, "%m", mainname );
   boost::replace_all( output_grid.filename, "%t", timestr  );
   boost::replace_all( output_non_grid.filename, "%f", filename );
   boost::replace_all( output_non_grid.filename, "%m", mainname );
   boost::replace_all( output_non_grid.filename, "%t", timestr  );
   boost::replace_all( binary_hint_file, "%f", filename );
   boost::replace_all( binary_hint_file, "%m", mainname );
   boost::replace_all( binary_hint_file, "%t", timestr  );

   // Logging parameters
   params.log_params.r_step    = readEntry<bool>( pt, "logging", "phi correction", false );
   params.log_params.r_end     = readEntry<bool>( pt, "logging", "phi correction final", false );
   params.log_params.divb_step = readEntry<bool>( pt, "logging", "divb correction", false );

   // Create the hint file that records the binary file structure
   if( output_grid.binary ){
      outputBinaryHintFile( binary_hint_file, params, output_grid );
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readParamsDivBCorrector
   ( boost::property_tree::ptree &pt
   , t_params &params
){
   std::string tempstr;

   tempstr = readEntry<std::string>( pt, "divb corrector", "method", "sor" );
   params.divb_method = fromString<DivBCorrectionMethod>( tempstr );
   switch( params.divb_method ){
   case DivBCorrectionMethod::Undefined:
      ERROUT << "ERROR: readParamsDivBCorrector: in section [divb corrector], key \"method\":\n"
             << "       Unknown correction method:" << tempstr << LF;
      exit( RET_ERR_WRONG_PARAMETER );
      break;
   case DivBCorrectionMethod::SOR:
      params.divb_sor_rsteps = readEntry<double>( pt, "divb corrector", "sor rsteps",  0       );
      params.divb_sor_rmax   = readEntry<double>( pt, "divb corrector", "sor rmax",    1.0e-5  );
      params.divb_sor_steps  = readEntry<int>(    pt, "divb corrector", "sor steps",   0       );
      params.divb_sor_max    = readEntry<double>( pt, "divb corrector", "sor divbmax", 1.0e-10 );
      params.divb_sor_omega_use_static = readEntry<bool>  ( pt, "divb corrector", "sor overrelaxation param use static", false );
      params.divb_sor_omega            = readEntry<double>( pt, "divb corrector", "sor overrelaxation param",            0.0   );
      break;
   }

   // How often to calculate the correction
   params.divb_skip_steps = readEntry<int>( pt, "divb corrector", "stepping rate", 1 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readProblemShockTube
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data
){
   std::string tempstr;

   // Domain size and position
   double Lx = readEntry<double>( pt, "problem", "Lx", 2.0 ); params.dx = Lx/params.nx;
   double Ly = readEntry<double>( pt, "problem", "Ly", 2.0 ); params.dy = Ly/params.ny;
   params.start_x = readEntry<double>( pt, "problem", "start_x", -1.0 );
   params.start_y = readEntry<double>( pt, "problem", "start_y", -1.0 );

   // Shock tube angle to the x-axis
   double angle = readEntry<double>( pt, "problem", "interface angle", 0.0 );

   double x0 = 0.5*Lx-0.5*params.dx;
   double y0 = 0.5*Ly-0.5*params.dy;
   double slope = tan( (90.0-angle) * M_PI/180.0 );

   double left[PRB_DIM];
   double right[PRB_DIM];

   left[0]  = readEntry<double>( pt, "problem", "rho_l", 1.0   );
   left[1]  = readEntry<double>( pt, "problem", "u_l",   0.0   );
   left[2]  = readEntry<double>( pt, "problem", "v_l",   0.0   );
   left[3]  = readEntry<double>( pt, "problem", "w_l",   0.0   );
   left[4]  = readEntry<double>( pt, "problem", "Bx_l",  0.75  );
   left[5]  = readEntry<double>( pt, "problem", "By_l",  1.0   );
   left[6]  = readEntry<double>( pt, "problem", "Bz_l",  0.0   );
   left[7]  = readEntry<double>( pt, "problem", "p_l",   1.0   );

   right[0] = readEntry<double>( pt, "problem", "rho_r", 0.125 );
   right[1] = readEntry<double>( pt, "problem", "u_r",   0.0   );
   right[2] = readEntry<double>( pt, "problem", "v_r",   0.0   );
   right[3] = readEntry<double>( pt, "problem", "w_r",   0.0   );
   right[4] = readEntry<double>( pt, "problem", "Bx_r",  0.75  );
   right[5] = readEntry<double>( pt, "problem", "By_r", -1.0   );
   right[6] = readEntry<double>( pt, "problem", "Bz_r",  0.0   );
   right[7] = readEntry<double>( pt, "problem", "p_r",   0.1   );

   for( int i = NXFIRST; i < NXLAST; i++ ){
      double line_y = y0 - slope*( (i-NXFIRST)*params.dx - x0 );
      double line_j = NYFIRST + (line_y/params.dy+0.5);
      int jump = line_j < NYFIRST ? NYFIRST :
                 line_j > NYLAST  ? NYLAST  :
                                    line_j;
      for( int j = NYFIRST; j < jump; j++ ){
         data.U[0][i][j] = left[0];
         data.u[0][i][j] = left[1];
         data.u[1][i][j] = left[2];
         data.u[2][i][j] = left[3];
         data.U[4][i][j] = left[4];
         data.U[5][i][j] = left[5];
         data.U[6][i][j] = left[6];
         data.p[i][j]    = left[7];
      }
      for( int j = jump; j < NYLAST; j++ ){
         data.U[0][i][j] = right[0];
         data.u[0][i][j] = right[1];
         data.u[1][i][j] = right[2];
         data.u[2][i][j] = right[3];
         data.U[4][i][j] = right[4];
         data.U[5][i][j] = right[5];
         data.U[6][i][j] = right[6];
         data.p[i][j]    = right[7];
      }
   }
   toConservationData( params, data );

   std::vector<std::string> boundary_name(params.b_count);
   boundary_name[params.b_left  ] = "left";
   boundary_name[params.b_top   ] = "top";
   boundary_name[params.b_right ] = "right";
   boundary_name[params.b_bottom] = "bottom";
   for( int b = 0; b < params.b_count; b++ ){
      std::string boundary_key = "boundary "+boundary_name[b];
      tempstr = readEntry<std::string>( pt, "problem", boundary_key, "open" );

      params.boundary[b] = fromString<BoundaryCondition>( tempstr );
      switch( params.boundary[b] ){
      case BoundaryCondition::Undefined:
         ERROUT << "ERROR: readProblemShockTube: in section [problem], key \"" << boundary_key << "\":\n"
                << "       Unknown boundary condition: " << tempstr << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case BoundaryCondition::Dirichlet:
         ERROUT << "ERROR: readProblemShockTube: in section [problem], key \"" << boundary_key << "\":\n"
                << "       Boundary condition not yet implemented: " << tempstr << LF;
         exit( RET_ERR_WRONG_PARAMETER );
   //      params.boundary_left_U = createVector( PRB_DIM );
   //      params.boundary_left_U[0] = readEntry<double>( pt, "problem", "boundary rho", data.U[0][NXFIRST] );
         break;
      case BoundaryCondition::Periodic:
      case BoundaryCondition::Neumann:
      case BoundaryCondition::Open:
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readProblemExplosion
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data
){
   std::string tempstr;

   // Domain size and position
   double Lx = readEntry<double>( pt, "problem", "Lx", 2.0 ); params.dx = Lx/params.nx;
   double Ly = readEntry<double>( pt, "problem", "Ly", 2.0 ); params.dy = Ly/params.ny;
   params.start_x = readEntry<double>( pt, "problem", "start_x", -1.0 );
   params.start_y = readEntry<double>( pt, "problem", "start_y", -1.0 );

   double in[PRB_DIM];
   double out[PRB_DIM];
   double r = readEntry<double>( pt, "problem", "r",       0.4   );

   in[0]    = readEntry<double>( pt, "problem", "rho_in",  1.0   );
   in[1]    = readEntry<double>( pt, "problem", "u_in",    0.0   );
   in[2]    = readEntry<double>( pt, "problem", "v_in",    0.0   );
   in[3]    = readEntry<double>( pt, "problem", "w_in",    0.0   );
   in[4]    = readEntry<double>( pt, "problem", "Bx_in",   0.0   );
   in[5]    = readEntry<double>( pt, "problem", "By_in",   0.0   );
   in[6]    = readEntry<double>( pt, "problem", "Bz_in",   0.0   );
   in[7]    = readEntry<double>( pt, "problem", "p_in",    1.0   );

   out[0]   = readEntry<double>( pt, "problem", "rho_out", 0.125 );
   out[1]   = readEntry<double>( pt, "problem", "u_out",   0.0   );
   out[2]   = readEntry<double>( pt, "problem", "v_out",   0.0   );
   out[3]   = readEntry<double>( pt, "problem", "w_out",   0.0   );
   out[4]   = readEntry<double>( pt, "problem", "Bx_out",  0.0   );
   out[5]   = readEntry<double>( pt, "problem", "By_out",  0.0   );
   out[6]   = readEntry<double>( pt, "problem", "Bz_out",  0.0   );
   out[7]   = readEntry<double>( pt, "problem", "p_out",   0.1   );

   for( int i = NXFIRST; i < NXLAST; i++ ){
      for( int j = NYFIRST; j < NYLAST; j++ ){
         if( (i-NX/2)*(i-NX/2)*params.dx*params.dx + (j-NY/2)*(j-NY/2)*params.dy*params.dy < r*r ){
            data.U[0][i][j] = in[0];
            data.u[0][i][j] = in[1];
            data.u[1][i][j] = in[2];
            data.u[2][i][j] = in[3];
            data.U[4][i][j] = in[4];
            data.U[5][i][j] = in[5];
            data.U[6][i][j] = in[6];
            data.p[i][j]    = in[7];
         } else {
            data.U[0][i][j] = out[0];
            data.u[0][i][j] = out[1];
            data.u[1][i][j] = out[2];
            data.u[2][i][j] = out[3];
            data.U[4][i][j] = out[4];
            data.U[5][i][j] = out[5];
            data.U[6][i][j] = out[6];
            data.p[i][j]    = out[7];
         }
      }
   }
   toConservationData( params, data );

   std::vector<std::string> boundary_name(params.b_count);
   boundary_name[params.b_left  ] = "left";
   boundary_name[params.b_top   ] = "top";
   boundary_name[params.b_right ] = "right";
   boundary_name[params.b_bottom] = "bottom";
   for( int b = 0; b < params.b_count; b++ ){
      std::string boundary_key = "boundary "+boundary_name[b];
      tempstr = readEntry<std::string>( pt, "problem", boundary_key, "open" );

      params.boundary[b] = fromString<BoundaryCondition>( tempstr );
      switch( params.boundary[b] ){
      case BoundaryCondition::Undefined:
         ERROUT << "ERROR: readProblemExplosion: in section [problem], key \"" << boundary_key << "\":\n"
                << "       Unknown boundary condition: " << tempstr << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case BoundaryCondition::Dirichlet:
         ERROUT << "ERROR: readProblemExplosion: in section [problem], key \"" << boundary_key << "\":\n"
                << "       Boundary condition not yet implemented: " << tempstr << LF;
         exit( RET_ERR_WRONG_PARAMETER );
   //      params.boundary_left_U = createVector( PRB_DIM );
   //      params.boundary_left_U[0] = readEntry<double>( pt, "problem", "boundary rho", data.U[0][NXFIRST] );
         break;
      case BoundaryCondition::Periodic:
      case BoundaryCondition::Neumann:
      case BoundaryCondition::Open:
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readProblemOTVortex
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data
){
   double x, y;

   // Domain size and position
   double Lx = 2.0*M_PI; params.dx = Lx/params.nx;
   double Ly = 2.0*M_PI; params.dy = Ly/params.ny;
   double x0 = readEntry<double>( pt, "problem", "x0 percent", 0.0 );
   double y0 = readEntry<double>( pt, "problem", "y0 percent", 0.0 );
   params.start_x = x0*2.0*M_PI;
   params.start_y = y0*2.0*M_PI;

   for( int i = NXFIRST; i < NXLAST; i++ ){
      for( int j = NYFIRST; j < NYLAST; j++ ){
         x       = params.start_x + (i-NXFIRST)*params.dx;
         y       = params.start_y + (j-NYFIRST)*params.dy;

         /* rho */ data.U[0][i][j] = params.gamma*params.gamma;
         /*  u  */ data.u[0][i][j] = -sin(y);
         /*  v  */ data.u[1][i][j] =  sin(x);
         /*  w  */ data.u[2][i][j] = 0;
         /*  Bx */ data.U[4][i][j] = -sin(y);
         /*  By */ data.U[5][i][j] =  sin(2*x);
         /*  Bz */ data.U[6][i][j] = 0;
         /*  p  */ data.p   [i][j] = params.gamma;
      }
   }
   toConservationData( params, data );

   for( int b = 0; b < params.b_count; b++ ){
      params.boundary[b] = BoundaryCondition::Periodic;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readProblemPlasmaSheet
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data
){
   std::string tempstr;

   // Domain size and position
   double Lx = readEntry<double>( pt, "problem", "Lx", 2.0 ); params.dx = Lx/params.nx;
   double Ly = readEntry<double>( pt, "problem", "Ly", 2.0 ); params.dy = Ly/params.ny;
   params.start_x = readEntry<double>( pt, "problem", "start_x", 0.0 );
   params.start_y = readEntry<double>( pt, "problem", "start_y", 0.0 );

   double left[PRB_DIM];
   double right[PRB_DIM];
   double up[PRB_DIM];
   double down[PRB_DIM];

   left[0]      = readEntry<double>( pt, "problem", "rho_l", 0.1 );
   left[1]      = readEntry<double>( pt, "problem", "u_l",   0.0 );
   left[2]      = readEntry<double>( pt, "problem", "v_l",   0.0 );
   left[3]      = readEntry<double>( pt, "problem", "w_l",   0.0 );
   left[4]      = readEntry<double>( pt, "problem", "Bx_l",  0.0 );
   left[5]      = readEntry<double>( pt, "problem", "By_l",  0.0 );
   left[6]      = readEntry<double>( pt, "problem", "Bz_l",  0.0 );
   left[7]      = readEntry<double>( pt, "problem", "p_l",   0.1 );

   right[0]     = readEntry<double>( pt, "problem", "rho_r", 1.0 );
   right[1]     = readEntry<double>( pt, "problem", "u_r",   0.0 );
   right[2]     = readEntry<double>( pt, "problem", "v_r",   0.0 );
   right[3]     = readEntry<double>( pt, "problem", "w_r",   0.0 );
   right[4]     = readEntry<double>( pt, "problem", "Bx_r",  0.0 );
   right[5]     = readEntry<double>( pt, "problem", "By_r",  0.0 );
   right[6]     = readEntry<double>( pt, "problem", "Bz_r",  0.0 );
   right[7]     = readEntry<double>( pt, "problem", "p_r",   1.0 );

   up[0]     = readEntry<double>( pt, "problem", "rho_u", 0.1   );
   up[1]     = readEntry<double>( pt, "problem", "u_u",   0.0   );
   up[2]     = readEntry<double>( pt, "problem", "v_u",   0.0   );
   up[3]     = readEntry<double>( pt, "problem", "w_u",   0.0   );
   up[4]     = readEntry<double>( pt, "problem", "Bx_u", -0.5   );
   up[5]     = readEntry<double>( pt, "problem", "By_u",  0.0   );
   up[6]     = readEntry<double>( pt, "problem", "Bz_u",  0.0   );
   up[7]     = readEntry<double>( pt, "problem", "p_u",   0.875 );

   down[0]     = readEntry<double>( pt, "problem", "rho_d", 0.1   );
   down[1]     = readEntry<double>( pt, "problem", "u_d",   0.0   );
   down[2]     = readEntry<double>( pt, "problem", "v_d",   0.0   );
   down[3]     = readEntry<double>( pt, "problem", "w_d",   0.0   );
   down[4]     = readEntry<double>( pt, "problem", "Bx_d",  0.5   );
   down[5]     = readEntry<double>( pt, "problem", "By_d",  0.0   );
   down[6]     = readEntry<double>( pt, "problem", "Bz_d",  0.0   );
   down[7]     = readEntry<double>( pt, "problem", "p_d",   0.875 );

   double sheet_start_head = readEntry<double>( pt, "problem", "sheet start", params.start_x+Lx/2 );
   double sheet_start_foot = readEntry<double>( pt, "problem", "sheet start foot", sheet_start_head );
   double sheet_thickness = readEntry<double>( pt, "problem", "sheet thickness", Ly/2.0 );

   int border_width = (Ly-sheet_thickness)/2.0/params.dy;
   int uplimit = params.ny-border_width + NYFIRST;
   int dnlimit = border_width + NYFIRST;

   std::vector<int> sheet_start_index(NY);

   double ellip_a  = sheet_start_head - sheet_start_foot;
   double ellip_b  = 0.5*(sheet_thickness-params.dy);
   double ellip_b2 = ellip_b*ellip_b;
   double ellip_x, ellip_y;
   for( int j = dnlimit; j < uplimit; j++ ){
      ellip_y = (j-1.5) * params.dy - Ly/2;
      ellip_x = sheet_start_foot + ellip_a*sqrt( 1.0 - (ellip_y*ellip_y/ellip_b2) );
      sheet_start_index[j] = (-params.start_x+ellip_x)/params.dx + NXFIRST;
   }

   std::string sheet_profile = readEntry<std::string>( pt, "problem", "sheet profile", "constant" );

   if( sheet_profile == "constant" ){
      for( int i = NXFIRST; i < NXLAST; i++ ){
         for( int j = NYFIRST; j < dnlimit; j++ ){
            data.U[0][i][j] = down[0];
            data.u[0][i][j] = down[1];
            data.u[1][i][j] = down[2];
            data.u[2][i][j] = down[3];
            data.U[4][i][j] = down[4];
            data.U[5][i][j] = down[5];
            data.U[6][i][j] = down[6];
            data.p[i][j]    = down[7];
         }
         for( int j = dnlimit; j < uplimit; j++ ){
            if( i < sheet_start_index[j] ){
               data.U[0][i][j] = left[0];
               data.u[0][i][j] = left[1];
               data.u[1][i][j] = left[2];
               data.u[2][i][j] = left[3];
               data.U[4][i][j] = left[4];
               data.U[5][i][j] = left[5];
               data.U[6][i][j] = left[6];
               data.p[i][j]    = left[7];
            } else {
               data.U[0][i][j] = right[0];
               data.u[0][i][j] = right[1];
               data.u[1][i][j] = right[2];
               data.u[2][i][j] = right[3];
               data.U[4][i][j] = right[4];
               data.U[5][i][j] = right[5];
               data.U[6][i][j] = right[6];
               data.p[i][j]    = right[7];
            }
         }
         for( int j = uplimit; j < NYLAST; j++ ){
            data.U[0][i][j] = up[0];
            data.u[0][i][j] = up[1];
            data.u[1][i][j] = up[2];
            data.u[2][i][j] = up[3];
            data.U[4][i][j] = up[4];
            data.U[5][i][j] = up[5];
            data.U[6][i][j] = up[6];
            data.p[i][j]    = up[7];
         }
      }
   } else if( sheet_profile == "tanh" ){
      double left_bb   = left[4]*left[4] + left[5]*left[5] + left[6]*left[6];
      double right_bb   = right[4]*right[4] + right[5]*right[5] + right[6]*right[6];

      for( int i = NXFIRST; i < NXLAST; i++ ){
         for( int j = NYFIRST; j < dnlimit; j++ ){
            data.U[0][i][j] = down[0];
            data.u[0][i][j] = down[1];
            data.u[1][i][j] = down[2];
            data.u[2][i][j] = down[3];
            data.U[4][i][j] = down[4];
            data.U[5][i][j] = down[5];
            data.U[6][i][j] = down[6];
            data.p[i][j]    = down[7];
         }
         for( int j = dnlimit; j < uplimit; j++ ){
            // tanh(y) profile for B, with -1 <= y <= 1
            double coef_b = tanh( 2.0*(j-dnlimit)/((uplimit-1)-dnlimit) - 1.0 );

            if( i < sheet_start_index[j] ){
               data.U[0][i][j] = left[0];
               data.u[0][i][j] = left[1];
               data.u[1][i][j] = left[2];
               data.u[2][i][j] = left[3];
               data.U[4][i][j] = left[4]*coef_b;
               data.U[5][i][j] = left[5]*coef_b;
               data.U[6][i][j] = left[6]*coef_b;
               // calculate pressure so that total pressure is kept constant
               data.p[i][j]    = left[7] + 0.5*left_bb*(1.0-coef_b*coef_b);
            } else {
               data.U[0][i][j] = right[0];
               data.u[0][i][j] = right[1];
               data.u[1][i][j] = right[2];
               data.u[2][i][j] = right[3];
               data.U[4][i][j] = right[4]*coef_b;
               data.U[5][i][j] = right[5]*coef_b;
               data.U[6][i][j] = right[6]*coef_b;
               // calculate pressure so that total pressure is kept constant
               data.p[i][j]    = right[7] + 0.5*right_bb*(1.0-coef_b*coef_b);
            }
         }
         for( int j = uplimit; j < NYLAST; j++ ){
            data.U[0][i][j] = up[0];
            data.u[0][i][j] = up[1];
            data.u[1][i][j] = up[2];
            data.u[2][i][j] = up[3];
            data.U[4][i][j] = up[4];
            data.U[5][i][j] = up[5];
            data.U[6][i][j] = up[6];
            data.p[i][j]    = up[7];
         }
      }
   }

   toConservationData( params, data );

   std::vector<std::string> boundary_name(params.b_count);
   boundary_name[params.b_left  ] = "left";
   boundary_name[params.b_top   ] = "top";
   boundary_name[params.b_right ] = "right";
   boundary_name[params.b_bottom] = "bottom";
   for( int b = 0; b < params.b_count; b++ ){
      std::string boundary_key = "boundary "+boundary_name[b];
      tempstr = readEntry<std::string>( pt, "problem", boundary_key, "open" );

      params.boundary[b] = fromString<BoundaryCondition>( tempstr );
      switch( params.boundary[b] ){
      case BoundaryCondition::Undefined:
         ERROUT << "ERROR: readProblemPlasmaSheet: in section [problem], key \"" << boundary_key << "\":\n"
                << "       Unknown boundary condition: " << tempstr << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case BoundaryCondition::Dirichlet:
         ERROUT << "ERROR: readProblemPlasmaSheet: in section [problem], key \"" << boundary_key << "\":\n"
                << "       Boundary condition not yet implemented: " << tempstr << LF;
         exit( RET_ERR_WRONG_PARAMETER );
   //      params.boundary_left_U = createVector( PRB_DIM );
   //      params.boundary_left_U[0] = readEntry<double>( pt, "problem", "boundary rho", data.U[0][NXFIRST] );
         break;
      case BoundaryCondition::Periodic:
      case BoundaryCondition::Neumann:
      case BoundaryCondition::Open:
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outputBinaryHintFile
   ( const std::string filename
   , const t_params &params
   , const t_output &output_grid
){
   t_output bin_hint;
   bin_hint.binary   = false;
   bin_hint.filename = filename;
   openFile( bin_hint );

   // doubles: x, y, 5 common (rho,Bx,By,Bz,divb), 4 natural (u,v,w,p), 4 conservation (mx,my,mz,e)
   int         field_count      = 2 + 5 + (output_grid.natural?4:0) + (output_grid.conservation?4:0);
   int         field_size       = output_grid.binary_double ? sizeof(double) : sizeof(float);
   const char *binary_data_type = output_grid.binary_double ? "double" : "float";
   int i = 0;
   fprintf( bin_hint.file, "data fields: %d)x  %d)y  %d)rho", i+1, i+2, i+3 ); i+=3;
   if( output_grid.natural ){
      fprintf( bin_hint.file, "  %d)u  %d)v  %d)w  %d)p", i+1, i+2, i+3, i+4 ); i+=4;
   }
   fprintf( bin_hint.file, "  %d)Bx  %d)By  %d)Bz %d)divB", i+1, i+2, i+3, i+4 ); i+=4;
   if( output_grid.conservation ){
      fprintf( bin_hint.file, "  %d)mx  %d)my  %d)mz  %d)e", i+1, i+2, i+3, i+4 ); i+=4;
   }
   fprintf( bin_hint.file, "\n\n" );

   fprintf( bin_hint.file, "nx = %d\nny = %d\n", params.nx, params.ny );
   fprintf( bin_hint.file, "xrange = [%.2f:%.2f]\nyrange = [%.2f:%.2f]\n",
                           params.start_x, params.start_x + params.nx*params.dx,
                           params.start_y, params.start_y + params.ny*params.dy );
   fprintf( bin_hint.file, "block_size = %d*%d\n", field_count, field_size );
   fprintf( bin_hint.file, "record = (ny,nx)\n" );
   fprintf( bin_hint.file, "skip   = index*nx*ny*block_size\n" );
   fprintf( bin_hint.file, "format = \"%%%d%s\"\n", field_count, binary_data_type );

   fprintf( bin_hint.file, "\nexample gnuplot command:\n" );
   fprintf( bin_hint.file, "datafile = \"%s\"\n", output_grid.filename.c_str() );
   fprintf( bin_hint.file, "index = 0\n" );
   fprintf( bin_hint.file, "nx = %d\nny = %d\n", params.nx, params.ny );
   fprintf( bin_hint.file, "set xrange [%.2f:%.2f]\nset yrange [%.2f:%.2f]\n",
                           params.start_x, params.start_x + params.nx*params.dx,
                           params.start_y, params.start_y + params.ny*params.dy );
   fprintf( bin_hint.file, "block_size = %d*%d\n", field_count, field_size );
   fprintf( bin_hint.file, "plot datafile binary record=(ny,nx) skip=index*nx*ny*block_size"
                           " format=\"%%%d%s\" u 1:2:3 w image\n", field_count, binary_data_type );

   closeFile( bin_hint );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outputGridDataBinary
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int    /*step*/ //unused
   , int    /*index*/ //unused
   , int    divb_nxfirst
   , int    divb_nxlast
   , int    divb_nyfirst
   , int    divb_nylast
   , double /*uumax*/ // unused
){
   // local aliases
   t_matrices U = data.U;
   t_matrices u = data.u;
   t_matrix   p = data.p;

   // fields: x, y, 5 common (rho,Bx,By,Bz,divb), 4 natural (u,v,w,p), 4 conservation (mx,my,mz,e)
   int field_count = 2 + 5 + (output.natural?4:0) + (output.conservation?4:0);

   // data records
   int      ri;
   size_t   record_count   = field_count*params.nx*params.ny;
   t_vector records_double = createVector( record_count );
   float   *records_float  = (float *)malloc( record_count*sizeof(*records_float) );

   // used in divb calculation
   const int kbx  = 4;
   const int kby  = 5;
   double    divb;

   if( output.binary_double ){
      ri = 0;
      for( int i = NXFIRST; i < NXLAST; i++ ){
         for( int j = NYFIRST; j < NYLAST; j++ ){
            // calculate divb
            if( ( i >= divb_nxfirst && i < divb_nxlast ) &&
                ( j >= divb_nyfirst && j < divb_nylast ) ){
               divb = fabs( (U[kbx][i+1][j]-U[kbx][i-1][j])/(2.0*params.dx)
                          + (U[kby][i][j+1]-U[kby][i][j-1])/(2.0*params.dy) );
            } else {
               divb = 0.0;
            }

            // fill record
            records_double[ri++] = params.start_x + (i-NXFIRST)*params.dx;
            records_double[ri++] = params.start_y + (j-NYFIRST)*params.dy;
            records_double[ri++] = U[0][i][j];     // rho
            if( output.natural ){
               records_double[ri++] = u[0][i][j];  // u
               records_double[ri++] = u[1][i][j];  // v
               records_double[ri++] = u[2][i][j];  // w
               records_double[ri++] = p[i][j];     // p
            }
            records_double[ri++] = U[4][i][j];     // Bx
            records_double[ri++] = U[5][i][j];     // By
            records_double[ri++] = U[6][i][j];     // Bz
            records_double[ri++] = divb;           // div B
            if( output.conservation ){
               records_double[ri++] = U[1][i][j];  // mx
               records_double[ri++] = U[2][i][j];  // my
               records_double[ri++] = U[3][i][j];  // mz
               records_double[ri++] = U[7][i][j];  // e
            }
         }
      }

      // write all of the records to file at once
      fwrite( records_double, sizeof(*records_double), record_count, output.file );
   } else {
      ri = 0;
      for( int i = NXFIRST; i < NXLAST; i++ ){
         for( int j = NYFIRST; j < NYLAST; j++ ){
            // calculate divb
            if( ( i >= divb_nxfirst && i < divb_nxlast ) &&
                ( j >= divb_nyfirst && j < divb_nylast ) ){
               divb = fabs( (U[kbx][i+1][j]-U[kbx][i-1][j])/(2.0*params.dx)
                          + (U[kby][i][j+1]-U[kby][i][j-1])/(2.0*params.dy) );
            } else {
               divb = 0.0;
            }

            // fill record
            records_float[ri++] = params.start_x + (i-NXFIRST)*params.dx;
            records_float[ri++] = params.start_y + (j-NYFIRST)*params.dy;
            records_float[ri++] = U[0][i][j];     // rho
            if( output.natural ){
               records_float[ri++] = u[0][i][j];  // u
               records_float[ri++] = u[1][i][j];  // v
               records_float[ri++] = u[2][i][j];  // w
               records_float[ri++] = p[i][j];     // p
            }
            records_float[ri++] = U[4][i][j];     // Bx
            records_float[ri++] = U[5][i][j];     // By
            records_float[ri++] = U[6][i][j];     // Bz
            records_float[ri++] = divb;           // div B
            if( output.conservation ){
               records_float[ri++] = U[1][i][j];  // mx
               records_float[ri++] = U[2][i][j];  // my
               records_float[ri++] = U[3][i][j];  // mz
               records_float[ri++] = U[7][i][j];  // e
            }
         }
      }

      // write all of the records to file at once
      fwrite( records_float, sizeof(*records_float), record_count, output.file );
   }

   free( records_float );
   freeVector( records_double );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outputGridDataText
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int    step
   , int    index
   , int    divb_nxfirst
   , int    divb_nxlast
   , int    divb_nyfirst
   , int    divb_nylast
   , double uumax
){
   // local aliases
   t_matrices U = data.U;
   t_matrices u = data.u;
   t_matrix   p = data.p;

   // header
   fprintf( output.file, "# record #%d, step #%d, t = %.5f, u_max = %.2f\n", index, step, data.t_current, sqrt(uumax) );
   int i = 0;
   fprintf( output.file, "# %d)i  %d)j  %d)x  %d)y  %d)rho", i+1, i+2, i+3, i+4, i+5 ); i+=5;
   if( output.natural ){
      fprintf( output.file, "  %d)u  %d)v  %d)w  %d)p", i+1, i+2, i+3, i+4 ); i+=4;
   }
   fprintf( output.file, "  %d)Bx  %d)By  %d)Bz %d)divB", i+1, i+2, i+3, i+4 ); i+=4;
   if( output.conservation ){
      fprintf( output.file, "  %d)mx  %d)my  %d)mz  %d)e", i+1, i+2, i+3, i+4 ); i+=4;
   }
   fprintf( output.file, "\n" );

   // used in divb calculation
   const int kbx = 4, kby = 5;
   double divb;

   // data
   for( int i = NXFIRST; i < NXLAST; i++ ){
      for( int j = NYFIRST; j < NYLAST; j++ ){
         if( ( i >= divb_nxfirst && i < divb_nxlast ) &&
             ( j >= divb_nyfirst && j < divb_nylast ) ){
            divb = fabs( (U[kbx][i+1][j]-U[kbx][i-1][j])/(2.0*params.dx)
                       + (U[kby][i][j+1]-U[kby][i][j-1])/(2.0*params.dy) );
         } else {
            divb = 0.0;
         }
         fprintf( output.file, "%d\t%d\t%+.5e\t%+.5e\t%+.5e",
                  i-NXFIRST, j-NYFIRST,
                  params.start_x + (i-NXFIRST)*params.dx, // x
                  params.start_y + (j-NYFIRST)*params.dy, // y
                  U[0][i][j]      // rho
         );
         if( output.natural )
            fprintf( output.file, "\t%+.5e\t%+.5e\t%+.5e\t%+.5e",
                     u[0][i][j],  // u
                     u[1][i][j],  // v
                     u[2][i][j],  // w
                     p[i][j]      // p
            );

         fprintf( output.file, "\t%+.5e\t%+.5e\t%+.5e\t%+.5e",
                  U[4][i][j],  // Bx
                  U[5][i][j],  // By
                  U[6][i][j],  // Bz
                  divb         // div B
         );
         if( output.conservation )
            fprintf( output.file, "\t%+.5e\t%+.5e\t%+.5e\t%+.5e",
                     U[1][i][j],  // mx
                     U[2][i][j],  // my
                     U[3][i][j],  // mz
                     U[7][i][j]   // e
            );
         fprintf( output.file, "\n" );
      }
      fprintf( output.file, "\n" );
   }
   fprintf( output.file, "\n\n" );
   fflush( output.file );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outputGridData
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int step
   , int index
){
   // local aliases
   t_matrices U = data.U;
   t_matrices u = data.u;

   // apply boundary conditions, required for div B calculation
   const int kbx = 4, kby = 5;
   int divb_nxfirst = NXFIRST;
   int divb_nxlast  = NXLAST;
   int divb_nyfirst = NYFIRST;
   int divb_nylast  = NYLAST;

   // left boundary
   switch( params.boundary[params.b_left] ){
   case BoundaryCondition::Undefined:
      break;
   case BoundaryCondition::Periodic:
      for( int k = kbx; k <= kby; k++ ){
         for( int j = NYFIRST; j < NYLAST; j++ ){
            U[k][NXFIRST-2][j] = U[k][NXLAST-2][j];
            U[k][NXFIRST-1][j] = U[k][NXLAST-1][j];
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      break;
   case BoundaryCondition::Neumann:
      divb_nxfirst = NXFIRST+1;
      break;
   case BoundaryCondition::Open:
      for( int k = kbx; k <= kby; k++ ){
         for( int j = NYFIRST; j < NYLAST; j++ ){
            U[k][NXFIRST-2][j] = U[k][NXFIRST][j];
            U[k][NXFIRST-1][j] = U[k][NXFIRST][j];
         }
      }
      divb_nxfirst = NXFIRST+1;
      break;
   }

   // right boundary
   switch( params.boundary[params.b_right] ){
   case BoundaryCondition::Undefined:
      break;
   case BoundaryCondition::Periodic:
      for( int k = kbx; k <= kby; k++ ){
         for( int j = NYFIRST; j < NYLAST; j++ ){
            U[k][NXLAST][j]    = U[k][NXFIRST][j];
            U[k][NXLAST+1][j]  = U[k][NXFIRST+1][j];
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      break;
   case BoundaryCondition::Neumann:
      divb_nxlast  = NXLAST-1;
      break;
   case BoundaryCondition::Open:
      for( int k = kbx; k <= kby; k++ ){
         for( int j = NYFIRST; j < NYLAST; j++ ){
            U[k][NXLAST][j]    = U[k][NXLAST-1][j];
            U[k][NXLAST+1][j]  = U[k][NXLAST-1][j];
         }
      }
      divb_nxlast  = NXLAST-1;
      break;
   }

   // bottom boundary
   switch( params.boundary[params.b_bottom] ){
   case BoundaryCondition::Undefined:
      break;
   case BoundaryCondition::Periodic:
      for( int k = kbx; k <= kby; k++ ){
         for( int i = NXFIRST; i < NXLAST; i++ ){
            U[k][i][NYFIRST-2] = U[k][i][NYLAST-2];
            U[k][i][NYFIRST-1] = U[k][i][NYLAST-1];
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      break;
   case BoundaryCondition::Neumann:
      divb_nyfirst = NYFIRST+1;
      break;
   case BoundaryCondition::Open:
      for( int k = kbx; k <= kby; k++ ){
         for( int i = NXFIRST; i < NXLAST; i++ ){
            U[k][i][NYFIRST-2] = U[k][i][NYFIRST];
            U[k][i][NYFIRST-1] = U[k][i][NYFIRST];
         }
      }
      divb_nyfirst = NYFIRST+1;
      break;
   }

   // top boundary
   switch( params.boundary[params.b_top] ){
   case BoundaryCondition::Undefined:
      break;
   case BoundaryCondition::Periodic:
      for( int k = kbx; k <= kby; k++ ){
         for( int i = NXFIRST; i < NXLAST; i++ ){
            U[k][i][NYLAST]    = U[k][i][NYFIRST];
            U[k][i][NYLAST+1]  = U[k][i][NYFIRST+1];
         }
      }
      break;
   case BoundaryCondition::Dirichlet:
      break;
   case BoundaryCondition::Neumann:
      divb_nylast  = NYLAST-1;
      break;
   case BoundaryCondition::Open:
      for( int k = kbx; k <= kby; k++ ){
         for( int i = NXFIRST; i < NXLAST; i++ ){
            U[k][i][NYLAST]    = U[k][i][NYLAST-1];
            U[k][i][NYLAST+1]  = U[k][i][NYLAST-1];
         }
      }
      divb_nylast  = NYLAST-1;
      break;
   }

   // find maximum velocity
   double uu, uumax = 0.0;
   for( int i = NXFIRST; i < NXLAST; i++ ){
      for( int j = NYFIRST; j < NYLAST; j++ ){
         uu = u[0][i][j]*u[0][i][j]
            + u[1][i][j]*u[1][i][j]
            + u[2][i][j]*u[2][i][j];
         if( uu > uumax ) uumax = uu;
      }
   }

   // binary or text output
   if( output.binary ){
      outputGridDataBinary( output, params, data, step, index, divb_nxfirst, divb_nxlast, divb_nyfirst, divb_nylast, uumax );
   } else {
      outputGridDataText(   output, params, data, step, index, divb_nxfirst, divb_nxlast, divb_nyfirst, divb_nylast, uumax );
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outputNonGridData
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int step
   , int index
){
   static bool initialized = false;
   t_matrices U = data.U;
   t_matrices u = data.u;
   t_matrix   p = data.p;

   // find maximum velocity and totals for energy, entropy, mass and momentum
   double totenergy = 0.0;
   double entropy, totentropy = 0.0;
   double totmass = 0.0;
   double totmomentum = 0.0;
   double uu, uumax = 0.0;
   for( int i = NXFIRST; i < NXLAST; i++ ){
      for( int j = NYFIRST; j < NYLAST; j++ ){
         uu = u[0][i][j]*u[0][i][j]
            + u[1][i][j]*u[1][i][j]
            + u[2][i][j]*u[2][i][j];
         if( uu > uumax ) uumax = uu;

         entropy = log( p[i][j]*pow( U[0][i][j], -params.gamma ) );
         totentropy  += entropy;
         totenergy   += U[7][i][j];
         totmass     += U[0][i][j];
         totmomentum += U[0][i][j]*sqrt(uu);
      }
   }
   // convert sums into integrals
   totentropy  *= params.dx*params.dy;
   totenergy   *= params.dx*params.dy;
   totmass     *= params.dx*params.dy;
   totmomentum *= params.dx*params.dy;

   // find max |divb|, integral |divb|
   double maxdivb = 0.0, totdivb = 0.0;
   divBCalculation( U, params, maxdivb, totdivb );

   // output header only once
   if( !initialized ){
      fprintf( output.file, "# 1)record_index 2)step 3)t 4)dt 5)u_max 6)|divB|_max 7)sum|divB| 8)total_mass 9)total_momentum 10)total_energy 11)total_entropy\n" );
      initialized = true;
   }

   // data
   fprintf( output.file, "%d\t%d\t%.5e\t%.5e\t%.2e\t%.2e\t%.2e\t%.8e\t%.8e\t%.8e\t%.8e\n",
            index, step, data.t_current, data.dt,
            sqrt(uumax), maxdivb, totdivb,
            totmass, totmomentum, totenergy, totentropy );

   // force write
   fflush( output.file );
}
