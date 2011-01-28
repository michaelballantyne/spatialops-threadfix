#ifndef SpatialOps_structured_StencilBuilder_h
#define SpatialOps_structured_StencilBuilder_h

/**
 *  \file StencilBuilder.h
 */

namespace SpatialOps{

  class OperatorDatabase;

namespace structured{

  /**
   *  \fn void build_stencils( const unsigned int nx,
   *                           const unsigned int ny,
   *                           const unsigned int nz,
   *                           const double Lx,
   *                           const double Ly,
   *                           const double Lz,
   *                           OperatorDatabase& opdb );
   *
   *  \param nx number of points in the x-direction
   *  \param ny number of points in the y-direction
   *  \param nz number of points in the z-direction
   *  \param Lx length in x-direction
   *  \param Ly length in y-direction
   *  \param Lz length in z-direction
   *  \param opdb the OperatorDatabase to register the operators on
   */
  void build_stencils( const unsigned int nx,
                       const unsigned int ny,
                       const unsigned int nz,
                       const double Lx,
                       const double Ly,
                       const double Lz,
                       OperatorDatabase& opdb );

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_StencilBuilder_h
