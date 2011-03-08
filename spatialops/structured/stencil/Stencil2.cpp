#include "Stencil2.h"

namespace SpatialOps{
namespace structured{


  template< typename OperatorT, typename SrcT, typename DestT >
  Stencil2<OperatorT,SrcT,DestT>::
  Stencil2( const double coefLo, const double coefHi )
    : coefLo_( coefLo ),
      coefHi_( coefHi )
  {}

  //------------------------------------------------------------------

  template< typename OperatorT, typename SrcT, typename DestT >
  Stencil2<OperatorT,SrcT,DestT>::
  ~Stencil2()
  {}

  //------------------------------------------------------------------

  template< typename OperatorT, typename SrcT, typename DestT >
  void
  Stencil2<OperatorT,SrcT,DestT>::
  apply_to_field( const SrcT& src, DestT& dest ) const
  {
    const Stencil2Helper<SrcT,DestT> helper( src.window_with_ghost(),
                                             dest.window_with_ghost() );

    const IntVec sinc = helper.src_increment();
    const IntVec dinc = helper.dest_increment();

    typename DestT::iterator idest = dest.begin() + helper.dest_offset();
    typename SrcT::const_iterator isrcm = src.begin() + helper.src_offset_lo();
    typename SrcT::const_iterator isrcp = src.begin() + helper.src_offset_hi();

    const IntVec lo = helper.low ();
    const IntVec hi = helper.high();

    for( int k=lo[2]; k<hi[2]; ++k ){
      for( int j=lo[1]; j<hi[1]; ++j ){
        for( int i=lo[0]; i<hi[0]; ++i ){
          *idest = coefLo_ * *isrcm +  coefHi_ * *isrcp;
          idest += dinc[0];
          isrcm += sinc[0];
          isrcp += sinc[0];
        }
        idest += dinc[1];
        isrcm += sinc[1];
        isrcp += sinc[1];
      }
      idest += dinc[2];
      isrcm += sinc[2];
      isrcp += sinc[2];
    }
  }

  //==================================================================
  // Explicit template instantiation
  //
#define DECLARE_STENCIL( OP, SRC, DEST )        \
  template class Stencil2< OP, SRC, DEST >;

# define DECLARE_BASIC_VARIANTS( VOL )					\
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::XFace )          \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::YFace )          \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::ZFace )          \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::XFace )          \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::YFace )          \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::ZFace )          \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::XFace, VOL )          \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::YFace, VOL )          \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::ZFace, VOL )
  
  DECLARE_BASIC_VARIANTS( SVolField );
  DECLARE_BASIC_VARIANTS( XVolField );
  DECLARE_BASIC_VARIANTS( YVolField );
  DECLARE_BASIC_VARIANTS( ZVolField );

  DECLARE_STENCIL( Interpolant, XVolField, YSurfXField )  // advecting velocity
  DECLARE_STENCIL( Gradient   , XVolField, YSurfXField )  // stress
  DECLARE_STENCIL( Interpolant, XVolField, ZSurfXField )  // advecting velocity
  DECLARE_STENCIL( Gradient   , XVolField, ZSurfXField )  // stress

  DECLARE_STENCIL( Interpolant, YVolField, XSurfYField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    YVolField, XSurfYField )  // stress
  DECLARE_STENCIL( Interpolant, YVolField, ZSurfYField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    YVolField, ZSurfYField )  // stress

  DECLARE_STENCIL( Interpolant, ZVolField, XSurfZField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    ZVolField, XSurfZField )  // stress
  DECLARE_STENCIL( Interpolant, ZVolField, YSurfZField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    ZVolField, YSurfZField )  // stress

  DECLARE_STENCIL( Interpolant, SVolField, XVolField )  // density, dp/dx
  DECLARE_STENCIL( Interpolant, SVolField, YVolField )  // density, dp/dy
  DECLARE_STENCIL( Interpolant, SVolField, ZVolField )  // density, dp/dz

  DECLARE_STENCIL( Interpolant, XVolField, SVolField )  // pressure projection RHS
  DECLARE_STENCIL( Interpolant, YVolField, SVolField )  // pressure projection RHS
  DECLARE_STENCIL( Interpolant, ZVolField, SVolField )  // pressure projection RHS

  DECLARE_STENCIL( Gradient, XVolField, SVolField )  // dilatation
  DECLARE_STENCIL( Gradient, YVolField, SVolField )  // dilatation
  DECLARE_STENCIL( Gradient, ZVolField, SVolField )  // dilatation

  DECLARE_STENCIL( Gradient, SVolField, XVolField )  // pressure
  DECLARE_STENCIL( Gradient, SVolField, YVolField )  // pressure
  DECLARE_STENCIL( Gradient, SVolField, ZVolField )  // pressure
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
