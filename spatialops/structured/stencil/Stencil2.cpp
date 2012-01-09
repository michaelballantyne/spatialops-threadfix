#include "Stencil2.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/FieldExpressionsStencil2.h>

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
      stencil_2_apply_to_field_general_execute<OperatorT,SrcT,DestT>(src,
                                                                     dest,
                                                                     coefLo_,
                                                                     coefHi_);
  }

  //==================================================================
  // Explicit template instantiation
  //
#define DECLARE_STENCIL( OP, SRC, DEST )        \
  template class Stencil2< OP, SRC, DEST >;

#define DECLARE_BASIC_VARIANTS( VOL )                          \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::XFace )   \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::YFace )   \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::ZFace )   \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::XFace )   \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::YFace )   \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::ZFace )   \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::XFace, VOL )   \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::YFace, VOL )   \
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

  DECLARE_STENCIL( Interpolant, SSurfXField, SVolField )  // ODT colocated mesh

  DECLARE_STENCIL( Interpolant, XSurfXField, XVolField )  // BC operator for tau_xx
  DECLARE_STENCIL( Interpolant, YSurfYField, YVolField )  // BC operator for tau_yy
  DECLARE_STENCIL( Interpolant, ZSurfZField, ZVolField )  // BC operator for tau_zz
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
