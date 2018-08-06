#ifndef SPATIALMETHODENOLF_HPP
#define SPATIALMETHODENOLF_HPP

#include "spatialmethodeno.hpp"

class SpatialMethodEnoLF : public SpatialMethodEno
{
   public:
      explicit SpatialMethodEnoLF( size_t       nx
                                 , size_t       ny
                                 , size_t       bufferWidth
                                 , double       dx
                                 , double       dy
                                 , t_boundary   boundary
                                 , double       gamma
                                 );
      virtual ~SpatialMethodEnoLF();

      t_status integrate( t_matrices    U
                        , t_matrices    UL
                        , double       &dtIdeal
                        ) override;

   protected:

   private:
      t_matrices F_;  // numerical flux F
      t_matrices G_;  // numerical flux G

      t_status getNumericalFluxF( double  U1[PRB_DIM]
                                , double  U2[PRB_DIM]
                                , double  F[PRB_DIM][8]
                                , double  alpha[PRB_DIM]
                                , double  F_[PRB_DIM]
                                , double &maxWaveSpeed );
      t_status getNumericalFluxG( double  U1[PRB_DIM]
                                , double  U2[PRB_DIM]
                                , double  G[PRB_DIM][8]
                                , double  alpha[PRB_DIM]
                                , double  G_[PRB_DIM]
                                , double &maxWaveSpeed );
};

#endif // SPATIALMETHODENOLF_HPP
