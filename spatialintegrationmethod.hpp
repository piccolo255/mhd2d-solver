#ifndef SPATIALINTEGRATIONMETHOD_HPP
#define SPATIALINTEGRATIONMETHOD_HPP

#include "mhd2d.hpp"

class SpatialIntegrationMethod
{
   public:
      explicit SpatialIntegrationMethod( size_t       nx
                                       , size_t       ny
                                       , size_t       bufferWidth
                                       , double       dx
                                       , double       dy
                                       , t_boundary   boundary
                                       , double       gamma
                                       );
      virtual ~SpatialIntegrationMethod();

      void initializeDirichletBoundaries( t_matrices U
                                        );

      virtual t_status integrate( t_matrices    U
                                , t_matrices    UL
                                , t_vectors     borderFluxLRUD
                                , double       &dtIdeal
                                ) = 0;

      static size_t requiredBufferWidth();

   protected:
      static const size_t minimumBufferWidth = 0;

      const size_t      nxProblem;
      const size_t      nyProblem;
      const size_t      bufferWidth;
      const double      dx;
      const double      dy;
      const t_boundary  boundary;
      const double      gamma;

      const size_t nxFirst;
      const size_t nyFirst;
      const size_t nxLast;
      const size_t nyLast;
      const size_t nxTotal;
      const size_t nyTotal;

      bool requireBoundaryInitialization;
      t_vectors dirichletBoundaryLeft   = nullptr;
      t_vectors dirichletBoundaryRight  = nullptr;
      t_vectors dirichletBoundaryTop    = nullptr;
      t_vectors dirichletBoundaryBottom = nullptr;

      virtual t_status applyBoundaryConditions( t_matrices U
                                              );

   private:
};

#endif // SPATIALINTEGRATIONMETHOD_HPP
