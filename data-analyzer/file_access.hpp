#ifndef FILE_ACCESS_HPP_INCLUDED
#define FILE_ACCESS_HPP_INCLUDED

// C++ headers
#include <string>

// Local headers
#include "data-analyzer.hpp"

// File access - open & close
void openFile
   ( t_output &output );

void closeFile
   ( t_output &output );

// File access - input from ini
void readConfig
   ( const std::string filename );

void readProblemConfig
   ( const std::string filename
   , t_output &output_grid
   , t_output &output_non_grid
   , t_params &params
   , t_data   &data );

#endif // FILE_ACCESS_HPP_INCLUDED
