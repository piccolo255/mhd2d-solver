/******************************************************************************
 *                                                                            *
 *                Numerical simulation of ideal MHD equations                 *
 *                             in two dimensions                              *
 *                                                                            *
 ******************************************************************************
 * Author: Tretler Rudolf                                                     *
 ******************************************************************************
 * Tatsuno Lab                                                                *
 * University of Electro-Communications, Tokyo                                *
 ******************************************************************************
 * Source code for research.                                                  *
 *                                                                            *
 *                                                                 2014-12-18 *
 ******************************************************************************/

#include "mhd2d.hpp"

#ifdef USE_EXCEPTIONS
/// http://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ThreadException::ThreadException
   (
): Ptr(nullptr) {
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ThreadException::~ThreadException
   (
){
   this->Rethrow();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ThreadException::Rethrow
   (
){
   if(this->Ptr) std::rethrow_exception(this->Ptr);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ThreadException::CaptureException
   (
){
   std::unique_lock<std::mutex> guard(this->Lock);
   this->Ptr = std::current_exception();
}
#endif // USE_EXCEPTIONS

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main
   ( int    argc
   , char **argv
){
   namespace nchrono = std::chrono;

   #ifdef DEBUG
      OUT << "Program arguments, argc = " << argc << LF;
      for( int i = 0; i < argc; i++ ){
         OUT << " - argv[" << i << "] = " << argv[i] << LF;
      }
   #endif // DEBUG

   /* Variable declarations */
   // output file information
   t_output output_grid;
   t_output output_non_grid;
   // parameters
   t_params params;
   // variables
   t_data data;

   /* Clocks */
   nchrono::steady_clock::time_point clstart, clmain;
   nchrono::duration<double> duration;
   clstart = nchrono::steady_clock::now();

   /* Check command line arguments */
   if( argc < 2 ){
      OUT << "Usage: " << argv[0] << " config_file\n";
      exit(0);
   }

   /* Initialize */
   inputData( argv[1], output_grid, output_non_grid, params, data );
   if( output_grid.single_file ){
      openFile( output_grid );
   }
   openFile( output_non_grid );
   outputGridData( output_grid, params, data, 0, 0 );
   outputNonGridData( output_non_grid, params, data, 0, 0 );

   int step_progress = 0;
   if( params.time_mode == TimeStepMode::Constant )
      step_progress = (params.steps>100) ? params.steps/100 : 1;

   // main loop variables
   bool output_to_screen, output_to_file, output_to_non_grid;
   bool done = false;
   int step = 1;
   int record_index = 1;
   int retval;
   double output_time_file   = output_grid.skip_t;
   double output_time_screen = params.t_max/100.0;
   double output_time_non_grid = output_non_grid.skip_t;

   /* main loop clock */
   clmain = nchrono::steady_clock::now();
   duration = clmain - clstart;
   OUT << "@ t = " << duration.count() << " sec : Init finished.\n";

   /* main loop */
   while( !done ){
      // step
      switch( params.time_stepping ){
      case TimeStepMethod::Undefined:
         ERROUT << "ERROR: main: Unknown time stepping method." << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case TimeStepMethod::Euler:
         retval = stepEuler( data.U, data.dt, params );
         break;
      case TimeStepMethod::RungeKutta3_TVD:
         retval = stepRK3TVD( data.U, data.dt, params );
         break;
      }

      // Process the return value
      if( retval == RET_ERR_TIME_UNDERFLOW ){
         ERROUT << "ERROR: main: Time step size smaller than the smallest allowed value.\n"
                << "       Reached at step #" << step << ", at t = " << data.t_current << ".\n"
                << "       The simulation will now terminate." << LF;
         break;
      }
      if( retval < 0 ){
         // other errors
         ERROUT << "ERROR: main: Terminating simulation." << LF;
         break;
      }

      // correct div B
      if( step%params.divb_skip_steps == 0 ){
         switch( params.divb_method ){
         case DivBCorrectionMethod::Undefined:
            break;
         case DivBCorrectionMethod::SOR:
            if( params.divb_sor_steps <= 0 ){
               break;
            }
            double divbsum, divbmax;
            for( int i = 0; i < params.divb_sor_steps; i++ ){
               divBCalculation( data.U, params, divbmax, divbsum );
               if( params.log_params.divb_step ){
                  OUT << "divb corrector, step #" << step << ", loop #" << i << ": divbmax = " << divbmax << "\n";
               }
               if( divbmax < params.divb_sor_max ){
                  if( params.log_params.divb_step ){
                     OUT << "divb corrector, step #" << step << ": divb converged in " << i << " steps.\n";
                  }
                  break;
               }
               retval = divBCorrectionSOR( data.U, params );
               if( retval == RET_ERR_NOT_CONVERGED ){
                  ERROUT << "WARNING: Div B corrector, r failed to converge.\n";
                  break;
               }
            }
            divBCalculation( data.U, params, divbmax, divbsum );
            if( divbmax > params.divb_sor_max ){
               ERROUT << "WARNING: Div B corrector failed to converge: divbmax = " << divbmax << "\n"
                      << "         Simulation resumed, but B may be erroneous." << LF;
               //break;
            }
            break;
         }
      }

      // update time
      data.t_current += data.dt;
#ifdef DEBUG_SHOW_SIMULATED_TIME
      OUT << "*** DEBUG: step = " << step << " / " << params.steps
          << ", t = " << data.t_current << LF;
#endif // DEBUG_SHOW_SIMULATED_TIME

      // time mode-specific behavior
      switch( params.time_mode ){
      case TimeStepMode::Undefined:
         ERROUT << "ERROR: main: Unknown time stepping mode." << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case TimeStepMode::Constant:
         done = step == params.steps;
         output_to_screen = step%step_progress == 0;
         break;
      case TimeStepMode::Variable:
         done = data.t_current >= params.t_max;
         output_to_screen = data.t_current >= output_time_screen;
         while( output_to_screen && ( output_time_screen <= data.t_current ) )
            output_time_screen += params.t_max/100.0;
         break;
      }

      // output mode-specific behavior
      switch( output_grid.skip_mode ){
      case WriteSkipMode::Undefined:
         ERROUT << "ERROR: main: Unknown grid output mode." << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case WriteSkipMode::Step:
         output_to_file = step%output_grid.skip_steps == 0;
         break;
      case WriteSkipMode::Time:
         output_to_file = data.t_current >= output_time_file;
         while( output_to_file && ( output_time_file <= data.t_current ) )
            output_time_file += output_grid.skip_t;
         break;
      }
      switch( output_non_grid.skip_mode ){
      case WriteSkipMode::Undefined:
         ERROUT << "ERROR: main: Unknown non-grid output mode." << LF;
         exit( RET_ERR_WRONG_PARAMETER );
         break;
      case WriteSkipMode::Step:
         output_to_non_grid = step%output_non_grid.skip_steps == 0;
         break;
      case WriteSkipMode::Time:
         output_to_non_grid = data.t_current >= output_time_non_grid;
         while( output_to_non_grid && ( output_time_non_grid <= data.t_current ) )
            output_time_non_grid += output_non_grid.skip_t;
         break;
      }

      // output data to file
      if( output_to_file || output_to_non_grid )
         toNaturalData( params, data );
      if( output_to_file )
         outputGridData( output_grid, params, data, step, record_index++ );
      if( output_to_non_grid )
         outputNonGridData( output_non_grid, params, data, step, record_index );


      // output progress report to screen
      if( output_to_screen ){
         duration = nchrono::steady_clock::now() - clstart;
         switch( params.time_mode ){
         case TimeStepMode::Undefined:
            ERROUT << "ERROR: main: Unknown time stepping mode." << LF;
            exit( RET_ERR_WRONG_PARAMETER );
            break;
         case TimeStepMode::Constant:
            OUT << "@ t = " << duration.count() << " sec : "
                << step << " / " << params.steps << " steps ("
                << (100.0*step/params.steps) << "%) done." << LF;
            break;
         case TimeStepMode::Variable:
            OUT << "@ t = " << duration.count() << " sec : "
                << step << " steps done, "
                << data.t_current << " / " << params.t_max << " simulated time ("
                << (100.0*data.t_current/params.t_max) << "%)." << LF;
            break;
         }
      }

      // stop if NaN encountered
      if( std::isnan( data.U[0][NXFIRST][NYFIRST] ) ){
         ERROUT << "WARNING: NaN detected. Stopping simulation at step " << step
                << ", t = " << data.t_current << ", with "
                << record_index-1 << " grids recorded." << LF;
         done = true;
      }

      step++;
   }

   /* Cleanup */
   closeFile( output_grid );
   closeFile( output_non_grid );

   // output program stats
   duration = nchrono::steady_clock::now() - clstart;
   OUT << "@ t = " << duration.count() << " sec : Simulation finished.\n"
       << " - Grids recorded     : " << record_index-1 << "\n"
       << " - Successful steps   : " << step-1 << "\n"
       << " - Simulated time     : " << data.t_current << "\n";
   duration = nchrono::steady_clock::now() - clmain;
   OUT << " - Main loop duration : " << duration.count() << " seconds\n"
       << " - Time step average  : " << duration.count()/(step-1) << " seconds.\n";

   return 0;
}
