#include "Stencil4.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/FieldExpressionsStencil4.h>

namespace SpatialOps{
namespace structured{

  //==================================================================

  template< typename OpT, typename SrcT, typename DestT >
  Stencil4<OpT,SrcT,DestT>::
  Stencil4( const double coef1,
            const double coef2,
            const double coef3,
            const double coef4 )
    : coef1_( coef1 ),
      coef2_( coef2 ),
      coef3_( coef3 ),
      coef4_( coef4 )
  {}

  //------------------------------------------------------------------

  template< typename OpT, typename SrcT, typename DestT >
  void
  Stencil4<OpT,SrcT,DestT>::
  apply_to_field( const SrcT& src, DestT& dest ) const
  {
      stencil_4_apply_to_field_general_execute<OpT,SrcT,DestT>(src,
                                                               dest,
                                                               coef1_,
                                                               coef2_,
                                                               coef3_,
                                                               coef4_);
  }

  //==================================================================
  // Explicit template instantiation
  //
#define DECLARE_STENCIL( OP, SRC, DEST )                                \
  template struct Stencil4< OP, SRC, DEST >;

  // viscosity from scalar cells to staggered surfaces for stress
  DECLARE_STENCIL( Interpolant, SVolField, XSurfYField ) 
  DECLARE_STENCIL( Interpolant, SVolField, XSurfZField )

  DECLARE_STENCIL( Interpolant, SVolField, YSurfXField )
  DECLARE_STENCIL( Interpolant, SVolField, YSurfZField )

  DECLARE_STENCIL( Interpolant, SVolField, ZSurfXField )
  DECLARE_STENCIL( Interpolant, SVolField, ZSurfYField )
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
