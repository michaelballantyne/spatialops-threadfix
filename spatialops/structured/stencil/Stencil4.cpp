#include "Stencil4.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \brief specialization of Stencil4Helper for SVol->XSurfY (viscosity)
   */
  template<> struct Stencil4Helper< SVolField, XSurfYField >
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    unsigned int src_offset_1() const{ return 0; }
    unsigned int src_offset_2() const{ return 1; }
    unsigned int src_offset_3() const{ return wsrc_.extent(0); }
    unsigned int src_offset_4() const{ return wsrc_.extent(0)+1; }

    unsigned int dest_offset() const{ return wdest_.extent(0)+1; }

    IntVec src_increment()  const{ return IntVec( 1, 1, wsrc_.extent(0) ); }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low()  const{ return IntVec(1,1,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, destInc_;
  };

  /**
   *  \brief specialization of Stencil4Helper for SVol->XSurfZ (viscosity)
   */
  template<> struct Stencil4Helper< SVolField, XSurfZField >
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        --hiBounds_[2];
      }
    }

    unsigned int src_offset_1() const{ return 0; }
    unsigned int src_offset_2() const{ return 1; }
    unsigned int src_offset_3() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    unsigned int src_offset_4() const{ return src_offset_3()+1; }

    unsigned int dest_offset() const{ return wdest_.extent(0)*wdest_.extent(1)+1; }

    IntVec src_increment()  const{ return IntVec( 1, 1, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 1, 0 ); }

    IntVec low()  const{ return IntVec(1,0,1); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, destInc_;
  };


  /**
   *  \brief specialization of Stencil4Helper for SVol->YSurfX (viscosity)
   */
  template<> struct Stencil4Helper< SVolField, YSurfXField >
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    unsigned int src_offset_1() const{ return 0; }
    unsigned int src_offset_2() const{ return 1; }
    unsigned int src_offset_3() const{ return wsrc_.extent(0); }
    unsigned int src_offset_4() const{ return src_offset_3()+1; }

    unsigned int dest_offset() const{ return wdest_.extent(0)+1; }

    IntVec src_increment()  const{ return IntVec( 1, 1, wsrc_.extent(0) ); }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low()  const{ return IntVec(1,1,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, destInc_;
  };

  /**
   *  \brief specialization of Stencil4Helper for SVol->YSurfZ (viscosity)
   */
  template<> struct Stencil4Helper< SVolField, YSurfZField >
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        --hiBounds_[0];
      }
    }

    unsigned int src_offset_1() const{ return 0; }
    unsigned int src_offset_2() const{ return wsrc_.extent(0); }
    unsigned int src_offset_3() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    unsigned int src_offset_4() const{ return src_offset_3()+wsrc_.extent(0); }

    unsigned int dest_offset() const{ return wdest_.extent(0)*wdest_.extent(1)+wdest_.extent(0); }

    IntVec src_increment()  const{ return IntVec( 1, 0, wsrc_.extent(0) ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, wdest_.extent(0) ); }

    IntVec low()  const{ return IntVec(0,1,1); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_;
  };

  /**
   *  \brief specialization of Stencil4Helper for SVol->ZSurfX (viscosity)
   */
  template<> struct Stencil4Helper< SVolField, ZSurfXField >
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    unsigned int src_offset_1() const{ return 0; }
    unsigned int src_offset_2() const{ return 1; }
    unsigned int src_offset_3() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    unsigned int src_offset_4() const{ return src_offset_3()+1; }

    unsigned int dest_offset() const{ return 1 + wdest_.extent(0)*wdest_.extent(1); }

    IntVec src_increment()  const{ return IntVec( 1, 1, 0 ); }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low()  const{ return IntVec(1,0,1); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, destInc_;
  };

  /**
   *  \brief specialization of Stencil4Helper for SVol->ZSurfY (viscosity)
   */
  template<> struct Stencil4Helper< SVolField, ZSurfYField >
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        destInc_( 1, 0, wdest.extent(0) )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    unsigned int src_offset_1() const{ return 0; }
    unsigned int src_offset_2() const{ return wsrc_.extent(0); }
    unsigned int src_offset_3() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    unsigned int src_offset_4() const{ return src_offset_3()+wsrc_.extent(0); }

    unsigned int dest_offset() const{ return wdest_.extent(0)*wdest_.extent(1) + wdest_.extent(0); }

    IntVec src_increment()  const{ return IntVec( 1, 0, wsrc_.extent(0) ); }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low()  const{ return IntVec(0,1,1); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, destInc_;
  };


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
    const Stencil4Helper<SrcT,DestT> helper( src.window_with_ghost(),
                                             dest.window_with_ghost() );

    const IntVec sinc = helper.src_increment();
    const IntVec dinc = helper.dest_increment();

    typename  DestT::iterator idest = dest.begin() + helper.dest_offset();

    typedef typename SrcT::const_iterator SrcIter;
    SrcIter isrc1 = src.begin() + helper.src_offset_1();
    SrcIter isrc2 = src.begin() + helper.src_offset_2();
    SrcIter isrc3 = src.begin() + helper.src_offset_3();
    SrcIter isrc4 = src.begin() + helper.src_offset_4();

    const IntVec lo = helper.low ();
    const IntVec hi = helper.high();

    for( int k=lo[2]; k<hi[2]; ++k ){
      for( int j=lo[1]; j<hi[1]; ++j ){
        for( int i=lo[0]; i<hi[0]; ++i ){
          *idest = coef1_ * *isrc1
                 + coef2_ * *isrc2
                 + coef3_ * *isrc3
                 + coef4_ * *isrc4;
          idest += dinc[0];
          isrc1 += sinc[0];
          isrc2 += sinc[0];
          isrc3 += sinc[0];
          isrc4 += sinc[0];
        }
        idest += dinc[1];
        isrc1 += sinc[1];
        isrc2 += sinc[1];
        isrc3 += sinc[1];
        isrc4 += sinc[1];
      }
      idest += dinc[2];
      isrc1 += sinc[2];
      isrc2 += sinc[2];
      isrc3 += sinc[2];
      isrc4 += sinc[2];
    }
  }


  //==================================================================
  // Explicit template instantiation
  //
#define DECLARE_STENCIL( OP, SRC, DEST )                                \
  template struct Stencil4< OP, SRC, DEST >;                            \
  template struct Stencil4Helper< SRC, DEST >;

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
