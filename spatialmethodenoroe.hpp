#ifndef SPATIALMETHODENOROE_HPP
#define SPATIALMETHODENOROE_HPP

#include "spatialmethodeno.hpp"

class SpatialMethodEnoRoe : public SpatialMethodEno
{
   public:
      explicit SpatialMethodEnoRoe( size_t      nx
                                  , size_t      ny
                                  , size_t      bufferWidth
                                  , double      dx
                                  , double      dy
                                  , t_boundary  boundary
                                  , double      gamma
                                  );
      virtual ~SpatialMethodEnoRoe();

      t_status integrate( t_matrices    U
                        , t_matrices    UL
                        , double       &dtIdeal
                        ) override;

   protected:

   private:
      t_matrices VUF; // undivided differences, flux F
      t_matrices VUG; // undivided differences, flux G
      t_matrices F_;  // numerical flux F
      t_matrices G_;  // numerical flux G

      t_status getNumericalFluxF( double  U1[PRB_DIM]
                                , double  U2[PRB_DIM]
                                , double  F[PRB_DIM][8]
                                , double  VUF[PRB_DIM][7]
                                , double  F_[PRB_DIM]
                                , double &maxWaveSpeed );
      t_status getNumericalFluxG( double  U1[PRB_DIM]
                                , double  U2[PRB_DIM]
                                , double  G[PRB_DIM][8]
                                , double  VUG[PRB_DIM][7]
                                , double  G_[PRB_DIM]
                                , double &maxWaveSpeed );
};

#endif // SPATIALMETHODENOROE_HPP
