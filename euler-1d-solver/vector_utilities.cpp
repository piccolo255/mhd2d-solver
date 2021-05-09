/******************************************************************************
 *                                                                            *
 *                 Helper functions for C vectors & matrices                  *
 *                                                                            *
 ******************************************************************************
 * Author: Tretler Rudolf                                                     *
 ******************************************************************************
 * Tatsuno Lab                                                                *
 * University of Electro-Communications, Tokyo                                *
 ******************************************************************************
 *                                                                 2015-01-13 *
 ******************************************************************************/

#include "vector_utilities.hpp"

t_vector create_vector( int size ){
   t_vector v;
   v = (t_vector)malloc( size*sizeof(*v) );
   return v;
}

void free_vector( t_vector v ){
   if( v ) free( v );
}

t_vectors create_vectors( int n, int size ){
   t_vector v;
   v = (t_vector)malloc( n*size*sizeof(*v) );
   t_vectors m;
   m = (t_vectors)malloc( n*sizeof(*m) );
   for( int i = 0; i < n; i++ )
      m[i] = v + i*size;
   return m;
}

void free_vectors( t_vectors m ){
   if( m ){
      if( m[0] ) free( m[0] );
      free( m );
   }
}

t_matrix create_matrix( int sizeX, int sizeY ){
   t_matrix m;
   m = (t_matrix)malloc( sizeX*sizeof(*m) );
   t_vector v;
   v = (t_vector)malloc( sizeX*sizeY*sizeof(*v) );
   for( int i = 0; i < sizeX; i++ )
      m[i] = v + i*sizeY;
   return m;
}

void free_matrix( t_matrix m ){
   if( m ){
      if( m[0] ) free( m[0] );
      free( m );
   }
}

t_matrices create_matrices( int n, int sizeX, int sizeY ){
   t_matrices mm;
   mm = (t_matrices)malloc( n*sizeof(*mm) );
   t_matrix m;
   m = (t_matrix)malloc( n*sizeX*sizeof(*m) );
   t_vector v;
   v = (t_vector)malloc( n*sizeX*sizeY*sizeof(*v) );
   for( int k = 0; k < n; k++ ){
      mm[k] = m + k*sizeX;
      for( int i = 0; i < sizeX; i++ )
         mm[k][i] = v + (k*sizeX+i)*sizeY;
   }
   return mm;
}

void free_matrices( t_matrices mm ){
   if( mm ){
      if( mm[0] ){
         if( mm[0][0] ) free( mm[0][0] );
         free( mm[0] );
      }
      free( mm );
   }
}
