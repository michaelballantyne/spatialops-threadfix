#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

namespace SpatialOps{
namespace structured{


  template< typename SrcT, typename DestT > IntVec high( const MemoryWindow& );  ///< upper bounds for mesh loop
  template< typename SrcT, typename DestT > IntVec  low( const MemoryWindow& );  ///< lower bounds for mesh loop

  template< typename SrcT, typename DestT > IntVec src_increment( const MemoryWindow& );  ///< increment count after each loop for source field
  template< typename SrcT, typename DestT > IntVec dest_increment( const MemoryWindow& );  ///< increment count after each loop for dest field

  template< typename SrcT, typename DestT > size_t src_offset_lo( const MemoryWindow& );
  template< typename SrcT, typename DestT > size_t src_offset_hi( const MemoryWindow& );

  template< typename SrcT, typename DestT > size_t dest_offset( const MemoryWindow& );


  /**
   *  \class Stencil
   *  \tparam OpT - the type of operator
   *  \tparam SrcT - the type of field the operator is applied to
   *  \tparam DestT - the type of field the operator produces
   */
  template< typename OperatorT, typename SrcFieldT, typename DestFieldT >
  class Stencil
  {
    const double coefLo_, coefHi_;
  public:

    typedef OperatorT  OpT;
    typedef SrcFieldT  SrcT;
    typedef DestFieldT DestT;

    Stencil( const double coefLo, const double coefHi )
      : coefLo_( coefLo ),
        coefHi_( coefHi )
    {}

    ~Stencil(){}

    void
    apply_to_field( const SrcT& src, DestT& dest ) const
    {
      const MemoryWindow& wsrc  = src.window_with_ghost();
      const MemoryWindow& wdest = dest.window_with_ghost();

      const IntVec sinc =  src_increment<SrcT,DestT>( wsrc  );
      const IntVec dinc = dest_increment<SrcT,DestT>( wdest );

      typename DestT::iterator idest = dest.begin() + dest_offset<SrcT,DestT>(wdest);
      typename SrcT::const_iterator isrcm = src.begin() + src_offset_lo<SrcT,DestT>(wsrc);
      typename SrcT::const_iterator isrcp = src.begin() + src_offset_hi<SrcT,DestT>(wsrc);

      const IntVec lo = low<SrcT,DestT>(wdest);
      const IntVec hi = high<SrcT,DestT>(wdest);

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
  };



  // --- SVol -> SSurfX
  template<> size_t dest_offset  <SVolField,SSurfXField>( const MemoryWindow& w ){ return 0; }
  template<> size_t src_offset_lo<SVolField,SSurfXField>( const MemoryWindow& w ){ return 0; }
  template<> size_t src_offset_hi<SVolField,SSurfXField>( const MemoryWindow& w ){ return 1; }

  template<> IntVec src_increment <SVolField,SSurfXField>( const MemoryWindow& w ){ return IntVec( w.extent(0)>1 ? 1 : 0, w.extent(1)>1 ? 1 : 0, w.extent(2)>1 ? 1 : 0 ); }
  template<> IntVec dest_increment<SVolField,SSurfXField>( const MemoryWindow& w ){ return IntVec( 1, 0, 0 ); }

  template<> IntVec low <SVolField,SSurfXField>( const MemoryWindow& w ){ return w.offset(); }
  template<> IntVec high<SVolField,SSurfXField>( const MemoryWindow& w ){ return w.extent(); }

  // --- SVol -> SSurfY
  template<> size_t dest_offset  <SVolField,SSurfYField>( const MemoryWindow& w ){ return w.extent(0); }
  template<> size_t src_offset_lo<SVolField,SSurfYField>( const MemoryWindow& w ){ return 0; }
  template<> size_t src_offset_hi<SVolField,SSurfYField>( const MemoryWindow& w ){ return w.extent(0); }

  template<> IntVec src_increment <SVolField,SSurfYField>( const MemoryWindow& w ){ return IntVec( w.extent(0)>1 ? 1 : 0, w.extent(1)>1 ? 1 : 0, w.extent(2)>1 ? 1 : 0 ); }
  template<> IntVec dest_increment<SVolField,SSurfYField>( const MemoryWindow& w ){ return IntVec( 1, 0, 0 ); }

  template<> IntVec low <SVolField,SSurfYField>( const MemoryWindow& w ){ return w.offset(); }
  template<> IntVec high<SVolField,SSurfYField>( const MemoryWindow& w ){ return w.extent() - IntVec(0,2,0); }

  // --- SVol -> SSurfZ
  template<> size_t dest_offset  <SVolField,SSurfZField>( const MemoryWindow& w ){ return w.extent(0); }
  template<> size_t src_offset_lo<SVolField,SSurfZField>( const MemoryWindow& w ){ return 0; }
  template<> size_t src_offset_hi<SVolField,SSurfZField>( const MemoryWindow& w ){ return w.extent(0)*w.extent(1); }

  template<> IntVec src_increment <SVolField,SSurfZField>( const MemoryWindow& w ){ return IntVec( 1, 0, 1 ); }
  template<> IntVec dest_increment<SVolField,SSurfZField>( const MemoryWindow& w ){ return IntVec( 1, 0, 0 ); }

  template<> IntVec low <SVolField,SSurfZField>( const MemoryWindow& w ){ return w.offset(); }
  template<> IntVec high<SVolField,SSurfZField>( const MemoryWindow& w ){ return w.extent() - IntVec(0,0,2); }



  // --- XVol -> XSurfX
  template<> size_t dest_offset  <XVolField,XSurfXField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfXField>(w); }
  template<> size_t src_offset_lo<XVolField,XSurfXField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfXField>(w); }
  template<> size_t src_offset_hi<XVolField,XSurfXField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfXField>(w); }

  template<> IntVec src_increment <XVolField,XSurfXField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfXField>(w); }
  template<> IntVec dest_increment<XVolField,XSurfXField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfXField>(w); }

  template<> IntVec low <XVolField,XSurfXField>( const MemoryWindow& w ){ return low<SVolField,SSurfXField>(w); }
  template<> IntVec high<XVolField,XSurfXField>( const MemoryWindow& w ){ return high<SVolField,SSurfXField>(w); }

  // --- XVol -> XSurfY
  template<> size_t dest_offset  <XVolField,XSurfYField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfYField>(w); }
  template<> size_t src_offset_lo<XVolField,XSurfYField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfYField>(w); }
  template<> size_t src_offset_hi<XVolField,XSurfYField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfYField>(w); }

  template<> IntVec src_increment <XVolField,XSurfYField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfYField>(w); }
  template<> IntVec dest_increment<XVolField,XSurfYField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfYField>(w); }

  template<> IntVec low <XVolField,XSurfYField>( const MemoryWindow& w ){ return low<SVolField,SSurfYField>(w); }
  template<> IntVec high<XVolField,XSurfYField>( const MemoryWindow& w ){ return high<SVolField,SSurfYField>(w); }

  // --- XVol -> XSurfZ
  template<> size_t dest_offset  <XVolField,XSurfZField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfZField>(w); }
  template<> size_t src_offset_lo<XVolField,XSurfZField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfZField>(w); }
  template<> size_t src_offset_hi<XVolField,XSurfZField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfZField>(w); }

  template<> IntVec src_increment <XVolField,XSurfZField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfZField>(w); }
  template<> IntVec dest_increment<XVolField,XSurfZField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfZField>(w); }

  template<> IntVec low <XVolField,XSurfZField>( const MemoryWindow& w ){ return low<SVolField,SSurfZField>(w); }
  template<> IntVec high<XVolField,XSurfZField>( const MemoryWindow& w ){ return high<SVolField,SSurfZField>(w); }



  // --- YVol -> YSurfX
  template<> size_t dest_offset  <YVolField,YSurfXField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfXField>(w); }
  template<> size_t src_offset_lo<YVolField,YSurfXField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfXField>(w); }
  template<> size_t src_offset_hi<YVolField,YSurfXField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfXField>(w); }

  template<> IntVec src_increment <YVolField,YSurfXField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfXField>(w); }
  template<> IntVec dest_increment<YVolField,YSurfXField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfXField>(w); }

  template<> IntVec low <YVolField,YSurfXField>( const MemoryWindow& w ){ return low<SVolField,SSurfXField>(w); }
  template<> IntVec high<YVolField,YSurfXField>( const MemoryWindow& w ){ return high<SVolField,SSurfXField>(w); }

  // --- YVol -> YSurfY
  template<> size_t dest_offset  <YVolField,YSurfYField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfYField>(w); }
  template<> size_t src_offset_lo<YVolField,YSurfYField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfYField>(w); }
  template<> size_t src_offset_hi<YVolField,YSurfYField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfYField>(w); }

  template<> IntVec src_increment <YVolField,YSurfYField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfYField>(w); }
  template<> IntVec dest_increment<YVolField,YSurfYField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfYField>(w); }

  template<> IntVec low <YVolField,YSurfYField>( const MemoryWindow& w ){ return low<SVolField,SSurfYField>(w); }
  template<> IntVec high<YVolField,YSurfYField>( const MemoryWindow& w ){ return high<SVolField,SSurfYField>(w); }

  // --- YVol -> YSurfZ
  template<> size_t dest_offset  <YVolField,YSurfZField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfZField>(w); }
  template<> size_t src_offset_lo<YVolField,YSurfZField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfZField>(w); }
  template<> size_t src_offset_hi<YVolField,YSurfZField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfZField>(w); }

  template<> IntVec src_increment <YVolField,YSurfZField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfZField>(w); }
  template<> IntVec dest_increment<YVolField,YSurfZField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfZField>(w); }

  template<> IntVec low <YVolField,YSurfZField>( const MemoryWindow& w ){ return low<SVolField,SSurfZField>(w); }
  template<> IntVec high<YVolField,YSurfZField>( const MemoryWindow& w ){ return high<SVolField,SSurfZField>(w); }



  // --- ZVol -> ZSurfX
  template<> size_t dest_offset  <ZVolField,ZSurfXField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfXField>(w); }
  template<> size_t src_offset_lo<ZVolField,ZSurfXField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfXField>(w); }
  template<> size_t src_offset_hi<ZVolField,ZSurfXField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfXField>(w); }

  template<> IntVec src_increment <ZVolField,ZSurfXField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfXField>(w); }
  template<> IntVec dest_increment<ZVolField,ZSurfXField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfXField>(w); }

  template<> IntVec low <ZVolField,ZSurfXField>( const MemoryWindow& w ){ return low<SVolField,SSurfXField>(w); }
  template<> IntVec high<ZVolField,ZSurfXField>( const MemoryWindow& w ){ return high<SVolField,SSurfXField>(w); }

  // --- ZVol -> ZSurfY
  template<> size_t dest_offset  <ZVolField,ZSurfYField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfYField>(w); }
  template<> size_t src_offset_lo<ZVolField,ZSurfYField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfYField>(w); }
  template<> size_t src_offset_hi<ZVolField,ZSurfYField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfYField>(w); }

  template<> IntVec src_increment <ZVolField,ZSurfYField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfYField>(w); }
  template<> IntVec dest_increment<ZVolField,ZSurfYField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfYField>(w); }

  template<> IntVec low <ZVolField,ZSurfYField>( const MemoryWindow& w ){ return low<SVolField,SSurfYField>(w); }
  template<> IntVec high<ZVolField,ZSurfYField>( const MemoryWindow& w ){ return high<SVolField,SSurfYField>(w); }

  // --- ZVol -> ZSurfZ
  template<> size_t dest_offset  <ZVolField,ZSurfZField>( const MemoryWindow& w ){ return dest_offset<SVolField,SSurfZField>(w); }
  template<> size_t src_offset_lo<ZVolField,ZSurfZField>( const MemoryWindow& w ){ return src_offset_lo<SVolField,SSurfZField>(w); }
  template<> size_t src_offset_hi<ZVolField,ZSurfZField>( const MemoryWindow& w ){ return src_offset_hi<SVolField,SSurfZField>(w); }

  template<> IntVec src_increment <ZVolField,ZSurfZField>( const MemoryWindow& w ){ return src_increment<SVolField,SSurfZField>(w); }
  template<> IntVec dest_increment<ZVolField,ZSurfZField>( const MemoryWindow& w ){ return dest_increment<SVolField,SSurfZField>(w); }

  template<> IntVec low <ZVolField,ZSurfZField>( const MemoryWindow& w ){ return low<SVolField,SSurfZField>(w); }
  template<> IntVec high<ZVolField,ZSurfZField>( const MemoryWindow& w ){ return high<SVolField,SSurfZField>(w); }

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
