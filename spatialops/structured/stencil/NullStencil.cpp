#include "NullStencil.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{

  template< typename OpT, typename SrcT, typename DestT >
  NullStencil<OpT,SrcT,DestT>::NullStencil()
  {}

  template< typename OpT, typename SrcT, typename DestT >
  void
  NullStencil<OpT,SrcT,DestT>::
  apply_to_field( const SrcT& src, DestT& dest ) const
  {
#   ifndef NDEBUG
    assert( src.window_with_ghost() == dest.window_with_ghost() );
#   endif
    typename SrcT::const_iterator isrc = src.begin();
    typename DestT::iterator idest = dest.begin();
    const typename DestT::iterator ideste = dest.end();
    for( ; idest!=ideste; ++isrc, ++idest ){
      *idest == *isrc;
    }
  }

  //==================================================================
  // Explicit template instantiation
  //
  template struct NullStencil< Interpolant, XVolField, SSurfXField >;  // x-advecting velocity on scalar x-surface
  template struct NullStencil< Interpolant, YVolField, SSurfYField >;  // y-advecting velocity on scalar y-surface
  template struct NullStencil< Interpolant, ZVolField, SSurfZField >;  // z-advecting velocity on scalar z-surface

  template struct NullStencil< Interpolant, SVolField, XSurfXField >;  // viscosity on x-x-face
  template struct NullStencil< Interpolant, SVolField, YSurfYField >;  // viscosity on y-y-face
  template struct NullStencil< Interpolant, SVolField, ZSurfZField >;  // viscosity on z-z-face
  //
  //==================================================================


} // namespace structured
} // namespace SpatialOps
