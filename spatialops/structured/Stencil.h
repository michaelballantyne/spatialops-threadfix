#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

namespace SpatialOps{
namespace structured{

  /**
   *  \struct Stencil2Helper
   *
   *  \brief Provides methods to adjust iterator positions when
   *         applying 2-point stencils to fields.  See also Stencil2
   */
  template< typename SrcT, typename DestT >
  struct Stencil2Helper
  {
    Stencil2Helper( const MemoryWindow& wsrc, ///< MemoryWindow for the source field
                    const MemoryWindow& wdest ///< MemoryWindow for the destination field
                    );
    IntVec high() const;           ///< upper bounds for mesh loop over dest field
    IntVec low() const;            ///< lower bounds for mesh loop over dest field
    IntVec src_increment() const;  ///< source field increment count after each loop (x,y,z)
    IntVec dest_increment() const; ///< destination field increment count after each loop (x,y,z)
    size_t src_offset_lo() const;  ///< offset for the "low" (minus-side) source field iterator
    size_t src_offset_hi() const;  ///< offset for the "high" (plus-side) source field iterator
    size_t dest_offset() const;    ///< offset for the destination field iterator (nominally 0)
  };

  /**
   *  \class Stencil2
   *
   *  \brief Support for implementing simple two-point stencils in
   *         one-dimension on structured meshes.
   *
   *  \tparam OpT - the type of operator
   *  \tparam SrcT - the type of field the operator is applied to
   *  \tparam DestT - the type of field the operator produces
   *
   *  See also Stencil2Helper
   */
  template< typename OperatorT, typename SrcFieldT, typename DestFieldT >
  class Stencil2
  {
    const double coefLo_, coefHi_;
  public:

    typedef OperatorT  OpT;
    typedef SrcFieldT  SrcT;
    typedef DestFieldT DestT;

    /**
     *  \brief construct a stencil with the specified coefficients
     *  \param coefLo the coefficient to multiply the (-) side field by
     *  \param coefHi the coefficient to multiply the (+) side field by
     */
    Stencil2( const double coefLo, const double coefHi )
      : coefLo_( coefLo ),
        coefHi_( coefHi )
    {}

    ~Stencil2(){}

    void
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
  };


  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from volume to x-surfaces
   *
   *  SVol -> SSurfX  (Grad, Interp)
   *  XVol -> XSurfX  (Grad, Interp)
   *  YVol -> YSurfX  (Grad, Interp)
   *  ZVol -> ZSurfX  (Grad, Interp)
   *
   *  \todo this breaks on dest_increment if we have physical boundaries...
   */
  template<typename VolT>
  struct Stencil2Helper< VolT,typename FaceTypes<VolT>::XFace >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest ) : wsrc_( wsrc ), wdest_( wdest ){}
    size_t dest_offset  () const{ return 1; }
    size_t src_offset_lo() const{ return 0; }
    size_t src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return IntVec( wsrc_.extent(0)>1 ? 1 : 0, wsrc_.extent(1)>1 ? 1 : 0, wsrc_.extent(2)>1 ? 1 : 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return wdest_.offset(); }
    IntVec high() const{ return wdest_.extent(); }
  private:
    const MemoryWindow &wsrc_, &wdest_;
  };


  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from volume to y-surfaces
   *
   *  SVol -> SSurfY  (Grad, Interp)
   *  XVol -> XSurfY  (Grad, Interp)
   *  YVol -> YSurfY  (Grad, Interp)
   *  ZVol -> ZSurfY  (Grad, Interp)
   *
   *  \todo this breaks on dest_increment if we have physical boundaries...
   */
  template<typename VolT>
  struct Stencil2Helper< VolT,typename FaceTypes<VolT>::YFace >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest ) : wsrc_( wsrc ), wdest_( wdest ){}
    size_t dest_offset  () const{ return 0; }
    size_t src_offset_lo() const{ return 0; }
    size_t src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return wdest_.offset(); }
    IntVec high() const{ return wdest_.extent(); }
  private:
    const MemoryWindow &wsrc_, &wdest_;
  };


  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from volume to z-surfaces
   *
   *  SVol -> SSurfZ  (Grad, Interp)
   *  XVol -> XSurfZ  (Grad, Interp)
   *  YVol -> YSurfZ  (Grad, Interp)
   *  ZVol -> ZSurfZ  (Grad, Interp)
   *
   *  \todo this breaks on dest_increment if we have physical boundaries...
   */
  template<typename VolT>
  struct Stencil2Helper< VolT,typename FaceTypes<VolT>::ZFace >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest ) : wsrc_( wsrc ), wdest_( wdest ){}
    size_t dest_offset  () const{ return 0; }
    size_t src_offset_lo() const{ return 0; }
    size_t src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return wdest_.offset(); }
    IntVec high() const{ return wdest_.extent(); }
  private:
    const MemoryWindow &wsrc_, &wdest_;
  };



  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from x-surfaces to volumes
   *
   *  SSurfX -> SVol  (Div)
   *  XSurfX -> XVol  (Div)
   *  YSurfX -> YVol  (Div)
   *  ZSurfX -> ZVol  (Div)
   *
   *  \todo this breaks on src_increment if we have physical boundaries.
   */
  template<typename VolT>
  struct Stencil2Helper< typename FaceTypes<VolT>::XFace, VolT >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest ) : wsrc_( wsrc ), wdest_( wdest ){}

    size_t dest_offset  () const{ return 0; }
    size_t src_offset_lo() const{ return 0; }
    size_t src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return wdest_.offset(); }
    IntVec high() const{ return wdest_.extent(); }

  private:
    const MemoryWindow &wsrc_, &wdest_;
  };

  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from y-surfaces to volumes
   *
   *  SSurfY -> SVol  (Div)
   *  XSurfY -> XVol  (Div)
   *  YSurfY -> YVol  (Div)
   *  ZSurfY -> ZVol  (Div)
   *
   *  \todo this breaks on src_increment if we have physical boundaries.
   */
  template<typename VolT>
  struct Stencil2Helper< typename FaceTypes<VolT>::YFace, VolT >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest ) : wsrc_( wsrc ), wdest_( wdest ){}

    size_t dest_offset  () const{ return 0; }
    size_t src_offset_lo() const{ return 0; }
    size_t src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return wdest_.offset(); }
    IntVec high() const{ return wdest_.extent(); }

  private:
    const MemoryWindow &wsrc_, &wdest_;
  };

  /**
   *  \brief Stencil2Helper for Stencil2 operators moving from z-surfaces to volumes
   *
   *  SSurfZ -> SVol  (Div)
   *  XSurfZ -> XVol  (Div)
   *  YSurfZ -> YVol  (Div)
   *  ZSurfZ -> ZVol  (Div)
   *
   *  \todo this breaks on src_increment if we have physical boundaries.
   */
  template<typename VolT>
  struct Stencil2Helper< typename FaceTypes<VolT>::ZFace, VolT >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest ) : wsrc_( wsrc ), wdest_( wdest ){}

    size_t dest_offset  () const{ return 0; }
    size_t src_offset_lo() const{ return 0; }
    size_t src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return wdest_.offset(); }
    IntVec high() const{ return wdest_.extent(); }

  private:
    const MemoryWindow &wsrc_, &wdest_;
  };


} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
