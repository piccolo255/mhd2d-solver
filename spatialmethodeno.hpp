#ifndef SPATIALMETHODENO_HPP
#define SPATIALMETHODENO_HPP

#include "spatialintegrationmethod.hpp"

class SpatialMethodEno : public SpatialIntegrationMethod
{
   public:
      explicit SpatialMethodEno( size_t      nx
                               , size_t      ny
                               , size_t      bufferWidth
                               , double      dx
                               , double      dy
                               , t_boundary  boundary
                               , double      gamma
                               );
      virtual ~SpatialMethodEno();

      static size_t requiredBufferWidth();

      bool getCharacteristicsX( t_matrices cx
                              , t_matrices LUx
                              ) override;
      bool getCharacteristicsY( t_matrices cy
                              , t_matrices LUy
                              ) override;

   protected:
      bool breakOnNegativePressure  = false;
      bool shownPressureWarning     = false;

      t_matrices F;
      t_matrices G;

      // x and y component of characteristic velocities;
      // see eqs. (2.92)&(2.94) in C.-W. Shu, "Essentially Non-Oscillatory and Weighted
      // Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws",
      // Lecture Notes in Mathematics 1697, pp.325-432, 1998
      // and Section 2.2 in R. Tretler, "Two-dimensional MHD simulation
      // of the piston-induced rarefaction wave in the Earth's plasma sheet",
      // Master's Thesis, University of Electro-Communications, 2015
      t_matrices _cx;
      t_matrices _cy;
      t_matrices _LUx;
      t_matrices _LUy;

      // Calculate horizontal and vertical physical fluxes, F and G, from physical values U
      void updateFluxes( const t_matrices U );

      // Eigenvalues lambda at point U
      t_status getEigenvaluesAt( const double   U[PRB_DIM]
                               , double         lambda[PRB_DIM] );
      // Eigenvalues lambda at midpoint between U1 and U2
      t_status getEigenvaluesBetween( const double U1[PRB_DIM]
                                    , const double U2[PRB_DIM]
                                    , double       lambda[PRB_DIM] );

      // Eigenvalues / eigenvectors between points 1 (left) and 2 (right) for x-direction flux F
      t_status getEigensF( double U1[PRB_DIM]
                         , double U2[PRB_DIM]
                         , double lambda[PRB_DIM]
                         , double lv[PRB_DIM][PRB_DIM]
                         , double rv[PRB_DIM][PRB_DIM] );
      // Eigenvalues / eigenvectors between points 1 (up) and 2 (down) for y-direction flux G
      t_status getEigensG( double U1[PRB_DIM]
                         , double U2[PRB_DIM]
                         , double lambda[PRB_DIM]
                         , double lv[PRB_DIM][PRB_DIM]
                         , double rv[PRB_DIM][PRB_DIM] );

   private:
      static const size_t minimumBufferWidth = 2;
};

#endif // SPATIALMETHODCENTRALFD_HPP
