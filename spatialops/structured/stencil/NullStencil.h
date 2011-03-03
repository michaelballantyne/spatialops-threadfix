#ifndef SpatialOps_structured_NullStencil_h
#define SpatialOps_structured_NullStencil_h

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
  template< typename OperatorT, typename SrcFieldT, typename DestFieldT >
  struct NullStencil
  {
    typedef OperatorT  OpT;
    typedef SrcFieldT  SrcFieldType;
    typedef DestFieldT DestFieldType;

    NullStencil();
    void apply_to_field( const SrcFieldT& src, DestFieldT& dest ) const;
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_NullStencil_h
