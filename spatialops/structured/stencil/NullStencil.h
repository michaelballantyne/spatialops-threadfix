#ifndef SpatialOps_structured_NullStencil.h
#define SpatialOps_structured_NullStencil

namespace SpatialOps{
namespace structured{

  /**
   *  \class NullStencil
   *  \author James C. Sutherland
   *
   *  \brief Direct copy operator.
   *
   *  For some operations, we are simply converting field types for
   *  uniform structured meshes.  Examples include:
   *   - interpolation from the scalar volume to staggered surfaces
   *   - interpolation from staggered volumes to scalar surfaces
   */
  template< typename OpT, typename SrcFieldT, typename DestFieldT >
  struct NullStencil
  {
    NullStencil();
    void apply_to_field( const SrcFieldT& src, DestFieldT& dest ) const;
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_NullStencil
