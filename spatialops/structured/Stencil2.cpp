#include "Stencil2.h"

#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{


  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from volume to x-surfaces
   *  \author James C. Sutherland
   *
   *  SVol -> SSurfX  (Grad, Interp)
   *  XVol -> XSurfX  (Grad, Interp)
   *  YVol -> YSurfX  (Grad, Interp)
   *  ZVol -> ZSurfX  (Grad, Interp)
   */
  template<typename VolT>
  struct Stencil2Helper< VolT, typename FaceTypes<VolT>::XFace >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    unsigned int dest_offset  () const{ return 1; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(1,0,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    IntVec srcInc_, destInc_, hiBounds_;
  };


  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from volume to y-surfaces
   *  \author James C. Sutherland
   *
   *  SVol -> SSurfY  (Grad, Interp)
   *  XVol -> XSurfY  (Grad, Interp)
   *  YVol -> YSurfY  (Grad, Interp)
   *  ZVol -> ZSurfY  (Grad, Interp)
   */
  template<typename VolT>
  struct Stencil2Helper< VolT, typename FaceTypes<VolT>::YFace >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, 0 ),
        hiBounds_( wdest.extent() )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }
    unsigned int dest_offset  () const{ return wdest_.extent(0); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,1,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec srcInc_, destInc_, hiBounds_;
  };


  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from volume to z-surfaces
   *  \author James C. Sutherland
   *
   *  SVol -> SSurfZ  (Grad, Interp)
   *  XVol -> XSurfZ  (Grad, Interp)
   *  YVol -> YSurfZ  (Grad, Interp)
   *  ZVol -> ZSurfZ  (Grad, Interp)
   */
  template<typename VolT>
  struct Stencil2Helper< VolT,typename FaceTypes<VolT>::ZFace >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        // physical boundary present
        --hiBounds_[2];
      }
    }
    unsigned int dest_offset  () const{ return wdest_.extent(0)*wdest_.extent(1); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return IntVec(0,0,1); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_;
  };



  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from x-surfaces to volumes
   *  \author James C. Sutherland
   *
   *  SSurfX -> SVol  (Div)
   *  XSurfX -> XVol  (Div)
   *  YSurfX -> YVol  (Div)
   *  ZSurfX -> ZVol  (Div)
   */
  template<typename VolT>
  struct Stencil2Helper< typename FaceTypes<VolT>::XFace, VolT >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : hiBounds_( wdest.extent() ),
        destInc_( 1, 0, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,0); }
    IntVec high() const{ return hiBounds_; }

  private:
    IntVec hiBounds_, destInc_;
  };

  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from y-surfaces to volumes
   *  \author James C. Sutherland
   *
   *  SSurfY -> SVol  (Div)
   *  XSurfY -> XVol  (Div)
   *  YSurfY -> YVol  (Div)
   *  ZSurfY -> ZVol  (Div)
   */
  template<typename VolT>
  struct Stencil2Helper< typename FaceTypes<VolT>::YFace, VolT >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ),  wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, wdest.extent(0) )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    unsigned int dest_offset  () const{ return wdest_.extent(0); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,1,0); }
    IntVec high() const{ return hiBounds_; }

  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, srcInc_, destInc_;
  };

  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from z-surfaces to volumes
   *  \author James C. Sutherland
   *
   *  SSurfZ -> SVol  (Div)
   *  XSurfZ -> XVol  (Div)
   *  YSurfZ -> YVol  (Div)
   *  ZSurfZ -> ZVol  (Div)
   */
  template<typename VolT>
  struct Stencil2Helper< typename FaceTypes<VolT>::ZFace, VolT >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        ++hiBounds_[2];
      }
    }

    unsigned int dest_offset  () const{ return wdest_.extent(0)*wdest_.extent(1); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return IntVec(0,0,1); }
    IntVec high() const{ return hiBounds_; }

  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_;
  };

  /**
   *  \brief Specialization for XVol->YSurfX (advecting velocity)
   */
  template<>
  struct Stencil2Helper< XVolField, YSurfXField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    unsigned int dest_offset  () const{ return wdest_.extent(0); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,1,0); }
    IntVec high() const{ return hiBounds_; }

  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, srcInc_, destInc_;
  };

  /**
   *  \brief Specialization for XVol->ZSurfX (advecting velocity)
   */
  template<>
  struct Stencil2Helper< XVolField, ZSurfXField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        destInc_( 1, 0, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    unsigned int dest_offset  () const{ return wdest_.extent(0) * wdest_.extent(1); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0) * wdest_.extent(1); }
    
    IntVec src_increment () const{ return IntVec(1,0,wsrc_.extent(0)); }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,1); }
    IntVec high() const{ return hiBounds_; }

  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, destInc_;
  };

  /**
   *  \brief Specialization for YVol->XSurfY (advecting velocity)
   */
  template<>
  struct Stencil2Helper< YVolField, XSurfYField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        --hiBounds_[0];
        destInc_[2] += wdest.extent(0);
      }
    }

    unsigned int dest_offset  () const{ return 1; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(1,0,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    IntVec hiBounds_, srcInc_, destInc_;
  };

  /**
   *  \brief Specialization for YVol->ZSurfY (advecting velocity)
   */
  template<>
  struct Stencil2Helper< YVolField, ZSurfYField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 0, 0 ),
        destInc_( 1, 0, 0 )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present - skip the extra face
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    unsigned int dest_offset  () const{ return wdest_.extent(0)*wdest_.extent(1); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,1); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, srcInc_, destInc_;
  };

  /**
   *  \brief Specialization for ZVol->XSurfZ (advecting velocity)
   */
  template<>
  struct Stencil2Helper< ZVolField, XSurfZField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        // physical boundary present
        --hiBounds_[2];
      }
    }

    unsigned int dest_offset  () const{ return 1; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(1,0,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, srcInc_, destInc_;
  };

  /**
   *  \brief Specialization for ZVol->YSurfZ (advecting velocity)
   */
  template<>
  struct Stencil2Helper< ZVolField, YSurfZField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, wdest.extent(0) )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        // physical boundary present
        --hiBounds_[2];
      }
    }

    unsigned int dest_offset  () const{ return wdest_.extent(0); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,1,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, srcInc_, destInc_;
  };

  //==================================================================

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
# define DECLARE_BASIC_VARIANTS( VOL )					\
  template class Stencil2< Interpolant, VOL, FaceTypes<VOL>::XFace >;	\
  template class Stencil2< Interpolant, VOL, FaceTypes<VOL>::YFace >;	\
  template class Stencil2< Interpolant, VOL, FaceTypes<VOL>::ZFace >;	\
  template class Stencil2< Gradient,    VOL, FaceTypes<VOL>::XFace >;	\
  template class Stencil2< Gradient,    VOL, FaceTypes<VOL>::YFace >;	\
  template class Stencil2< Gradient,    VOL, FaceTypes<VOL>::ZFace >;	\
  template class Stencil2< Divergence, FaceTypes<VOL>::XFace, VOL >;	\
  template class Stencil2< Divergence, FaceTypes<VOL>::YFace, VOL >;	\
  template class Stencil2< Divergence, FaceTypes<VOL>::ZFace, VOL >;
  
  DECLARE_BASIC_VARIANTS( SVolField );
  DECLARE_BASIC_VARIANTS( XVolField );
  DECLARE_BASIC_VARIANTS( YVolField );
  DECLARE_BASIC_VARIANTS( ZVolField );

  template class Stencil2< Interpolant, XVolField, YSurfXField >;  // advecting velocity
  template class Stencil2< Interpolant, XVolField, ZSurfXField >;  // advecting velocity

  template class Stencil2< Interpolant, YVolField, XSurfYField >;  // advecting velocity
  template class Stencil2< Interpolant, YVolField, ZSurfYField >;  // advecting velocity

  template class Stencil2< Interpolant, ZVolField, XSurfZField >;  // advecting velocity
  template class Stencil2< Interpolant, ZVolField, YSurfZField >;  // advecting velocity
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
