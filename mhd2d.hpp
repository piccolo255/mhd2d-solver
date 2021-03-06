#ifndef MHD2D_HPP_INCLUDED
#define MHD2D_HPP_INCLUDED

// debug groups, uncomment needed ones
#ifdef DEBUG
//#define DEBUG_MIDPOINT_EIGENVALUES
//#define DEBUG_MAX_VELOCITY
//#define DEBUG_DT_UPDATE
//#define DEBUG_SHOW_SIMULATED_TIME
//#define DEBUG_FLUX_LF
#endif // DEBUG

// C headers
#include <cmath>
#include <cstdio>

// C++ headers
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include <chrono>

// Headers required for exception handling
#ifdef USE_THREAD_EXCEPTIONS
#include <thread>
#include <mutex>
#endif // USE_THREAD_EXCEPTIONS

// Headers required for OpenMP parallelization
#ifdef OPENMP
#include <omp.h>
#endif // OPENMP

// Boost headers
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

// Local headers
#include "vector_utilities.hpp"
#include "enums.hpp"

// OUT and IN can be redefined as a filestream
// to enable direct file input/output
// In MinGW already defined in windef.h, so first undefine
#ifdef OUT
#undef OUT
#endif // OUT
#define OUT std::cout
#ifdef IN
#undef IN
#endif
#define IN  std::cin
#define ERROUT std::cerr
#define LF  std::endl
#define ABS std::abs

// Problem size and number of velocity dimensions to keep track of,
// to avoid magic numbers in loops
#define PRB_DIM 8
#define VEL_DIM 3

// Array sizes
//  - data is in [2..n+1] with [0], [1], [n+2] and [n+3] as buffers
#define NX      (params.nx+4)
#define NY      (params.ny+4)
#define NXFIRST (2)
#define NYFIRST (2)
#define NXLAST  (params.nx+2)
#define NYLAST  (params.ny+2)

#define EPS            (1e-20)
#define EPS_EQUAL(a,b) (std::abs((a)-(b))<EPS)
#define EPS_ZERO(x)    (std::abs(x)<EPS)
#define MAX(a,b)       ((a)>(b) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

typedef struct {
   // Binary or text
   bool binary;
   bool binary_double;

   // Single or per-step
   bool single_file;

   // Output categories
   bool natural;
   bool conservation;

   // Skipping parameters
   WriteSkipMode skip_mode;
   int    skip_steps;
   double skip_t;

   // Mask formats
   std::string format_simtime;
   std::string format_index;
   std::string format_step;

   // File
   //TODO
   std::string base_filename;
   std::string filename;
   FILE *file;
} t_output;

typedef struct {
   bool r_step;
   bool r_end;
   bool divb_step;
   bool characteristics;
} t_log_params;

typedef struct {
   bool           isError;
   ReturnStatus   status;
   std::string    message;
} t_status;

typedef struct {
   BoundaryCondition right;
   BoundaryCondition top;
   BoundaryCondition left;
   BoundaryCondition bottom;
} t_boundary;

typedef struct {
   // Grid parameters
   int    nx;
   int    ny;
   double dx;
   double dy;
   double start_x;
   double start_y;

   // Time parameters
   TimeStepMode time_mode;
   int    steps;
   double cfl_number;
   double dt_min;
   double dt_max;
   double t_max;

   // Simulation parameters
   const int b_count  = 4;
   const int b_left   = 0;
   const int b_bottom = 1;
   const int b_right  = 2;
   const int b_top    = 3;
   BoundaryCondition boundary[4];
   t_vectors         boundary_dirichlet_U[4];
   TimeStepMethod    time_stepping;
   IntegrationMethod scheme;

   // Physical parameters
   double gamma;

   // Div B corrector parameters
   DivBCorrectionMethod divb_method;
   int    divb_skip_steps;
   double divb_sor_rmax;
   double divb_sor_rsteps;
   double divb_sor_max;
   int    divb_sor_steps;
   bool   divb_sor_omega_use_static;
   double divb_sor_omega;

   bool break_on_neg_pressure;

   // logging parameters
   t_log_params log_params;
} t_params;

typedef struct {
   t_vectors left;
   t_vectors right;
   t_vectors up;
   t_vectors down;
} borderVectors;

typedef struct {
   // Conserved variables; density, momentum, magnetic field, energy
   t_matrices U;
   // Primitive variables; velocity and pressure
   t_matrices u;
   t_matrix   p;
   // Characteristic velocities
   t_matrices cx;
   t_matrices cy;
   // Characteristic variables
   t_matrices LUx;
   t_matrices LUy;
   // Border fluxes
   borderVectors borderFlux;
   // Time
   double dt;
   double t_current;
} t_data;

// Utility - error handling
void criticalError
   ( ReturnStatus       error
   , const std::string  message );

// File access - open & close
void openFile
   ( t_output &output );

void closeFile
   ( t_output &output );

// File access - input from ini
void inputData
   ( const std::string filename
   , t_output &output_grid
   , t_output &output_non_grid
   , t_params &params
   , t_data   &data );

// File access - output grid data file
void outputGridData
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int step
   , int index );

// File access - output non-grid data to file
void outputNonGridData
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int step
   , int index );

// Output characteristics in binary format
void outputCharacteristicsBinary
   ( const t_output    &output
   , const t_params    &params
   , const t_matrices  &cx
   , const t_matrices  &cy
   , const t_matrices  &LUx
   , const t_matrices  &LUy
   , double time
   , int    step
   , int    index );

// Conversion - conservation to natural variables
void toNaturalPoint
   ( const t_params &params
   , double U[PRB_DIM]
   , double u[VEL_DIM]
   , double &p );

// Conversion - natural to conservation variables
void toConservationPoint
   ( const t_params &params
   , double U[PRB_DIM]
   , double u[VEL_DIM]
   , double &p );

// Conversion - conservation to natural variables
void toNaturalData
   ( const t_params &params
   , t_data &data );

// Conversion - natural to conservation variables
void toConservationData
   ( const t_params &params
   , t_data &data );

// Time stepper - Third order optimal TVD Runge-Kutta time stepping method
t_status stepRK3TVD
   ( t_matrices  U
   , double     &dt
   , const t_params &params );

// Time stepper - Euler time stepping method
t_status stepEuler
   ( t_matrices  U
   , double     &dt
   , const t_params &params );

// Scheme - central FD
t_status methodCentralFD
   ( t_matrices  U
   , t_matrices  UL
   , double     &dt_step
   , const t_params &params );

// Scheme - ENO-Roe with characteristics decomposition
t_status methodENOSystem
   ( t_matrices  U
   , t_matrices  UL
   , double     &dt_step
   , const t_params &params );

// Div B correction - calculated in place, using SOR (Successive Over-Relaxation)
t_status divBCorrectionSOR
   ( t_matrices U
   , const t_params &params );

// Div B calculation - maximum and total
void divBCalculation
   ( t_matrices U
   , const t_params &params
   , double &maxdivb
   , double &totdivb );

#ifdef USE_THREAD_EXCEPTIONS
// Class for exception handling that can be used in OpenMP parallelized code
// Source: http://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp
class ThreadException {
   std::exception_ptr Ptr;
   std::mutex         Lock;

public:
   ThreadException
      ( );
   ~ThreadException
      ( );

   void Rethrow
      ( );

   void CaptureException
      ( );

   template <typename Function, typename... Parameters>
   void ThreadException::Run
      ( Function f
      , Parameters... params
   ){
      try {
         f( params... );
      } catch( ... ) {
         CaptureException();
      }
   }
};
#endif // USE_THREAD_EXCEPTIONS

#endif // MHD2D_HPP_INCLUDED
