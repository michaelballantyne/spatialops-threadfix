#include "Stencil4.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>

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
    typedef s4detail::ExtentsAndOffsets<SrcT,DestT> Extents;

    const MemoryWindow& wsrc  =  src.window_with_ghost();
    const MemoryWindow& wdest = dest.window_with_ghost();

    const MemoryWindow ws1( wsrc.glob_dim(),
                            wsrc.offset() + Extents::Src1Offset::int_vec(),
                            wsrc.extent() + Extents::Src1Extent::int_vec() + wsrc.has_bc()*Extents::Src1ExtentBC::int_vec(),
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );

    const MemoryWindow ws2( wsrc.glob_dim(),
                            wsrc.offset() + Extents::Src2Offset::int_vec(),
                            wsrc.extent() + Extents::Src2Extent::int_vec() + wsrc.has_bc()*Extents::Src2ExtentBC::int_vec(),
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );

    const MemoryWindow ws3( wsrc.glob_dim(),
                            wsrc.offset() + Extents::Src3Offset::int_vec(),
                            wsrc.extent() + Extents::Src3Extent::int_vec() + wsrc.has_bc()*Extents::Src3ExtentBC::int_vec(),
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );

    const MemoryWindow ws4( wsrc.glob_dim(),
                            wsrc.offset() + Extents::Src4Offset::int_vec(),
                            wsrc.extent() + Extents::Src4Extent::int_vec() + wsrc.has_bc()*Extents::Src4ExtentBC::int_vec(),
                            wsrc.has_bc(0), wsrc.has_bc(1), wsrc.has_bc(2) );

    const MemoryWindow  wd( wdest.glob_dim(),
                            wdest.offset() + Extents::DestOffset::int_vec(),
                            wdest.extent() + Extents::DestExtent::int_vec() + wsrc.has_bc()*Extents::DestExtentBC::int_vec(),
                            wdest.has_bc(0), wdest.has_bc(1), wdest.has_bc(2) );

    // ensure that all field extents are equal - a minimum requirement for sanity here.
    assert( ws1.extent() == ws2.extent() &&
            ws3.extent() == ws4.extent() &&
            ws1.extent() == wd.extent()  );

    // build new windowed fields
    DestT  d( wd, &dest[0], ExternalStorage );
    SrcT  s1( ws1, &src[0], ExternalStorage );
    SrcT  s2( ws2, &src[0], ExternalStorage );
    SrcT  s3( ws3, &src[0], ExternalStorage );
    SrcT  s4( ws4, &src[0], ExternalStorage );

    // implement the stencil
    typename DestT::iterator      id =d .begin();
    typename DestT::iterator      ide=d .end();
    typename SrcT::const_iterator is1=s1.begin(), is2=s2.begin(), is3=s3.begin(), is4=s4.begin();
    for( ; id!=ide; ++id, ++is1, ++is2, ++is3, ++is4 ){
      *id = *is1 *coef1_ + *is2*coef2_ + *is3*coef3_ + *is4*coef4_;
    }

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
