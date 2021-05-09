#ifndef VECTOR_UTILITIES_HPP_INCLUDED
#define VECTOR_UTILITIES_HPP_INCLUDED

#include <cstdlib>

typedef double * t_vector;
typedef t_vector * t_vectors;

typedef t_vector * t_matrix;
typedef t_matrix * t_matrices;

// create vector & group of vectors
t_vector create_vector(            int size );
t_vectors create_vectors(   int n, int size );
// create matrix & group of matrices
t_matrix create_matrix(            int sizeX, int sizeY );
t_matrices create_matrices( int n, int sizeX, int sizeY );

// free
void free_vector(   t_vector v    );
void free_vectors(  t_vectors m   );
void free_matrix(   t_matrix m    );
void free_matrices( t_matrices mm );

#endif // VECTOR_UTILITIES_HPP_INCLUDED
