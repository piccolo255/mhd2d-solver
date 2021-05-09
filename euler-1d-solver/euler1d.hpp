#ifndef EULER1D_HPP_INCLUDED
#define EULER1D_HPP_INCLUDED

// C includes
#include <cmath>
#include <cstdio>

// C++ includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string>
#include <chrono>

// Boost includes
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

// Local includes
#include "vector_utilities.hpp"

// OUT and IN can be redefined as a filestream
// to enable direct file input/output
// In MinGW already defined in windef.h, so need to undef them first
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
#define PRB_DIM 3
#define VEL_DIM 1

// Array sizes
//  - data is in [2..n+1] with [0], [1], [n+2] and [n+3] as buffers
#define NX      (params.nx+4)
#define NXFIRST (2)
#define NXLAST  (params.nx+2)

// Output modes
#define OUT_MODE_STEP 1
#define OUT_MODE_TIME 2

// Time stepping modes
#define TIME_MODE_CONSTANT 1
#define TIME_MODE_VARIABLE 2

// Time stepping methods
#define STEP_EULER   1
#define STEP_RK3TVD  2
#define STEP_RIEMANN 3

// Spatial integration methods
#define METHOD_CENTRAL_FD 1
#define METHOD_ENO        2

// Boundary conditions
#define BOUNDARY_OPEN        1
#define BOUNDARY_PERIODIC    2

// Special problem types
#define PROBLEM_NORMAL  1
#define PROBLEM_PISTON  2
#define PROBLEM_RIEMANN 3

// Function return values
#define RET_OK              0
#define RET_NO_CHANGE       1
#define RET_UPDATED         2
#define RET_ERR_TIME_UNDERFLOW         -1
#define RET_RIEMANN_FAILED_TO_CONVERGE -2

#define EPS            (1e-20)
#define EPS_EQUAL(a,b) (std::abs((a)-(b))<EPS)
#define EPS_ZERO(x)    (std::abs(x)<EPS)
#define MAX(a,b)       ((a)>(b) ? (a) : (b))

typedef struct {
   // Output categories
   bool natural;
   bool conservation;

   // Skipping parameters
   int    skip_mode;
   int    skip_steps;
   double skip_t;

   // File
   std::string filename;
   FILE *file;
} t_output;

typedef struct {
   // Grid parameters
   int    nx;
   double dx;
   double start_x;

   // Time parameters
   int    time_mode;
   int    steps;
   double cfl_number;
   double dt_min;
   double dt_max;
   double t_max;

   // Simulation parameters
   int boundary;
   int time_stepping;
   int scheme;

   // Physical parameters
   double gamma;

   // Problem-specific parameters
   int problem_type;
   int param_dbl_n;
   int param_int_n;
   t_vector param_dbl;
   int     *param_int;
} t_params;

typedef struct {
   // Conserved variables; density, momentum, magnetic field, energy
   t_vectors U;
   // Primitive variables; velocity and pressure
   t_vectors u;
   t_vector  p;
   double dt;
   double t_current;
} t_data;

// File access - open & close
void openFile(  t_output &output );
void closeFile( t_output  output );

// File access - input from ini
void inputData( const std::string filename,
                t_output &output_grid, t_output &output_non_grid, t_params &params, t_data &data );

// File access - output grid data file
void outputGridData( const t_output &output, const t_params &params, const t_data &data, int step, int index );

// File access - output non-grid data to file
void outputNonGridData( const t_output &output, const t_params &params, const t_data &data, int step, int index );

// Conversion - conservation to natural variables
void toNatural( const t_params &params,
                t_data &data );

// Conversion - natural to conservation variables
void toConservation( const t_params &params,
                     t_data &data );

// Time stepper - Third order optimal TVD Runge-Kutta time stepping method
int rk3tvd( t_vectors U, double &dt, const t_params &params );

// Time stepper - Euler time stepping method
int euler_step( t_vectors U, double &dt, const t_params &params );

// Solver - Riemann solver
int riemann_solver( t_data &data, const t_params &params );

// Scheme - central FD
void central_fd( t_vectors U, t_vectors UL, double &dt_step,
                 const t_params &params );

// Scheme - ENO-Roe with characteristics decomposition
void eno_system_roe( t_vectors U, t_vectors UL, double &dt_step,
                     const t_params &params );

#endif // EULER1D_HPP_INCLUDED
