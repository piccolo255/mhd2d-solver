#ifndef FILE_ACCESS_HPP_INCLUDED
#define FILE_ACCESS_HPP_INCLUDED

// Read value of a key, show warning and default value if the key is missing
template <class T>
T readEntry
   ( boost::property_tree::ptree pt
   , std::string section
   , std::string name
   , T           defaultValue );

// Read settings from property tree: Div B corrector
void readParamsDivBCorrector
   ( boost::property_tree::ptree &pt
   , t_params &params );

// Read settings from property tree: Shock tube problem
void readProblemShockTube
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data );

// Read settings from property tree: Explosion problem
void readProblemExplosion
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data );

// Read settings from property tree: Orszag-Tang vortex problem
void readProblemOTVortex
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data );

// Read settings from property tree: Plasma sheet problem
void readProblemPlasmaSheet
   ( boost::property_tree::ptree &pt
   , t_params &params
   , t_data   &data );

// Print info about the binary file structure (with hints for plotting with GNUplot)
void outputBinaryHintFile
   ( const std::string filename
   , const t_params &params
   , const t_output &output_grid );

// Output grid data file in binary format
void outputGridDataBinary
   ( const t_output &output
   , const t_params &params
   , const t_data   &data
   , int    step
   , int    index
   , int    divb_nxfirst
   , int    divb_nxlast
   , int    divb_nyfirst
   , int    divb_nylast
   , double uumax );

// Output grid data file in textual format
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
   , double uumax );

#endif // FILE_ACCESS_HPP_INCLUDED
