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
};

#endif // SPATIALMETHODENOROE_HPP
