#include "Stencil2.h"
#include <spatialops/FieldExpressions.h>

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
    typedef s2detail::ExtentsAndOffsets<SrcT,DestT> Extents;

    const MemoryWindow& ws = src.window_with_ghost();

    const MemoryWindow ws1( ws.glob_dim(),
                            ws.offset() + Extents::Src1Offset::int_vec(),
                            ws.extent() + Extents::Src1Extent::int_vec() + ws.has_bc()*Extents::Src1ExtentBC::int_vec(),
                            ws.has_bc(0), ws.has_bc(1), ws.has_bc(2) );

    const MemoryWindow ws2( ws.glob_dim(),
                            ws.offset() + Extents::Src2Offset::int_vec(),
                            ws.extent() + Extents::Src2Extent::int_vec() + ws.has_bc()*Extents::Src2ExtentBC::int_vec(),
                            ws.has_bc(0), ws.has_bc(1), ws.has_bc(2) );

    const MemoryWindow& wdest = dest.window_with_ghost();

    const MemoryWindow wd( wdest.glob_dim(),
                           wdest.offset() + Extents::DestOffset::int_vec(),
                           wdest.extent() + Extents::DestExtent::int_vec() + wdest.has_bc()*Extents::DestExtentBC::int_vec(),
                           wdest.has_bc(0), wdest.has_bc(1), wdest.has_bc(2) );

//    std::cout << "apply_to_field2 info: " << std::endl
//        << wdest << std::endl
//        << wdest.has_bc() *Extents::DestExtentBC::int_vec() << std::endl
//        << "s1 : " << ws1 << std::endl
//        << "s2 : " << ws2 << std::endl
//        << "d  : " << wd  << std::endl
//        << "s1o: " << Extents::Src1Offset::print() << std::endl
//        << "s2o: " << Extents::Src2Offset::print() << std::endl
//        << "do : " << Extents::DestOffset::print() << std::endl
//        << "s2e: " << Extents::Src2Extent::print() << std::endl
//        << "de : " << Extents::DestExtent::print() << std::endl
//        << "s2bcaug: " << Extents::Src2ExtentBC::print() << std::endl
//        << "d bcaug: " << Extents::DestExtentBC::print() << std::endl;

    assert( ws1.extent() == ws2.extent() && ws1.extent() == wd.extent() );

    // build fields using these newly created windows to do the stencil operation.
    DestT  d( wd, &dest[0], ExternalStorage );
    SrcT  s1( ws1, &src[0], ExternalStorage );
    SrcT  s2( ws2, &src[0], ExternalStorage );

    typename DestT::iterator      id  = d .begin();
    typename DestT::iterator      ide = d .end();
    typename SrcT::const_iterator is1 = s1.begin();
    typename SrcT::const_iterator is2 = s2.begin();
    for( ; id!=ide; ++id, ++is1, ++is2 ){
      *id = *is1 * coefLo_ + *is2 * coefHi_;
    }
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
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
