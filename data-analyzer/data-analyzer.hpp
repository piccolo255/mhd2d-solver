#ifndef DATA_ANALYZER_HPP_INCLUDED
#define DATA_ANALYZER_HPP_INCLUDED

// debug groups, uncomment needed ones
#ifdef DEBUG
#endif // DEBUG

// C++ headers
#include <string>

// Local headers
#include "vector_utilities.hpp"
#include "enums.hpp"

// Problem size and number of velocity dimensions to keep track of,
// to avoid magic numbers in loops
const size_t PRB_DIM = 8;
const size_t VEL_DIM = 3;

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
   bool           isError;
   ReturnStatus   status;
   std::string    message;
} t_status;

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
   double dt_min;
   double dt_max;
   double t_max;
} t_params;

typedef struct {
   // Conserved variables; density, momentum, magnetic field, energy
   t_matrices U;
   // Primitive variables; velocity and pressure
   t_matrices u;
   t_matrix   p;
   double dt;
   double t_current;
} t_data;

// Utility - error handling
void criticalError
   ( ReturnStatus       error
   , const std::string  message );

#endif // DATA_ANALYZER_HPP_INCLUDED
