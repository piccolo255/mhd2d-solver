#ifndef SCHEME_ENO_H_INCLUDED
#define SCHEME_ENO_H_INCLUDED

// Eigenvalues lambda at point U
t_status getEigenvalues
   ( double U[PRB_DIM]
   , double gamma
   , double lambda[PRB_DIM]
   , bool   break_on_neg_pressure );

// Eigenvalues lambda at midpoint between U1 and U2
t_status getEigenvalues
   ( double U1[PRB_DIM]
   , double U2[PRB_DIM]
   , double gamma
   , double lambda[PRB_DIM]
   , bool   break_on_neg_pressure );

// Eigenvalues / eigenvectors between points 1 (left) and 2 (right) for x-direction flux F
t_status getEigens_F
   ( double U1[PRB_DIM]
   , double U2[PRB_DIM]
   , double gamma
   , double lambda[PRB_DIM]
   , double lv[PRB_DIM][PRB_DIM]
   , double rv[PRB_DIM][PRB_DIM]
   , bool   break_on_neg_pressure );

// Eigenvalues / eigenvectors between points 1 (up) and 2 (down) for y-direction flux G
t_status getEigens_G
   ( double U1[PRB_DIM]
   , double U2[PRB_DIM]
   , double gamma
   , double lambda[PRB_DIM]
   , double lv[PRB_DIM][PRB_DIM]
   , double rv[PRB_DIM][PRB_DIM]
   , bool   break_on_neg_pressure );

// Boundary conditions
void applyBoundaryConditions
   ( t_matrices U
   , const t_params &params );

// Calculate horizontal and vertical physical fluxes, F and G, from physical values U
void getFluxes
   ( t_matrices U
   , t_matrices F
   , t_matrices G
   , const t_params &params );

// Use ENO-Roe method to find numerical flux in the x direction (flux F)
t_status findNumericalFluxRoe_F
   ( double  U1[PRB_DIM]
   , double  U2[PRB_DIM]
   , double  F[PRB_DIM][8]
   , double  VUF[PRB_DIM][7]
   , double  F_[PRB_DIM]
   , double &a_max
   , int     first_valid
   , int     last_valid
   , const t_params &params );

// Use ENO-Roe method to find numerical flux in the y direction (flux G)
t_status findNumericalFluxRoe_G
   ( double  U1[PRB_DIM]
   , double  U2[PRB_DIM]
   , double  G[PRB_DIM][8]
   , double  VUG[PRB_DIM][7]
   , double  G_[PRB_DIM]
   , double &a_max
   , int     first_valid
   , int     last_valid
   , const t_params &params );

// Use ENO-LF method to find numerical flux in the x direction (flux F)
t_status findNumericalFluxLF_F
   ( double  U1[PRB_DIM]
   , double  U2[PRB_DIM]
   , double  U[PRB_DIM][8]
   , double  F[PRB_DIM][8]
   , double  F_[PRB_DIM]
   , double &a_max
   , double  alpha[PRB_DIM]
   , int     first_valid
   , int     last_valid
   , const t_params &params );

// Use ENO-LF method to find numerical flux in the y direction (flux G)
t_status findNumericalFluxLF_G
   ( double  U1[PRB_DIM]
   , double  U2[PRB_DIM]
   , double  U[PRB_DIM][8]
   , double  G[PRB_DIM][8]
   , double  G_[PRB_DIM]
   , double &a_max
   , double  alpha[PRB_DIM]
   , int     first_valid
   , int     last_valid
   , const t_params &params );

#endif // SCHEME_ENO_H_INCLUDED
