#ifndef VECTOR_UTILITIES_HPP_INCLUDED
#define VECTOR_UTILITIES_HPP_INCLUDED

typedef double   * t_vector;
typedef t_vector * t_vectors;

typedef t_vector * t_matrix;
typedef t_matrix * t_matrices;

// create vector & group of vectors
t_vector createVector
   ( int size );

t_vectors createVectors
   ( int n
   , int size );

// create matrix & group of matrices
t_matrix createMatrix
   ( int sizeX
   , int sizeY );

t_matrices createMatrices
   ( int n
   , int sizeX
   , int sizeY );

// free
void freeVector
   ( t_vector v );

void freeVectors
   ( t_vectors m );

void freeMatrix
   ( t_matrix m );

void freeMatrices
   ( t_matrices mm );

#endif // VECTOR_UTILITIES_HPP_INCLUDED
