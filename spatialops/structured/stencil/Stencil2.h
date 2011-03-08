#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \class Stencil2
   *  \author James C. Sutherland
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
    typedef SrcFieldT  SrcFieldType;
    typedef DestFieldT DestFieldType;

    /**
     *  \brief construct a stencil with the specified coefficients
     *  \param coefLo the coefficient to multiply the (-) side field by
     *  \param coefHi the coefficient to multiply the (+) side field by
     */
    Stencil2( const double coefLo, const double coefHi );

    ~Stencil2();

    void apply_to_field( const SrcFieldType& src, DestFieldType& dest ) const;

    double get_minus_coef() const{ return coefLo_; }
    double get_plus_coef () const{ return coefHi_; }
  };





  /**
   *  \struct Stencil2Helper
   *  \author James C. Sutherland
   *
   *  \brief Provides methods to adjust iterator positions when
   *         applying 2-point stencils to fields.  See also Stencil2.
   *
   *  Note that the Stencil2Helper declaration below is not
   *  implemented.  Only template specializations are implemented.
   *  This just serves as a way to show the required methods that a
   *  specialization should implement.
   */
  template< typename SrcT, typename DestT >
  struct Stencil2Helper;



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

    ~Stencil2Helper(){}

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
        destInc_( 1, 0, wdest.extent(0) ),
        hiBounds_( wdest.extent() )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    ~Stencil2Helper(){}

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

    ~Stencil2Helper(){}

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
      : hiBounds_( wdest.extent() - IntVec(1,0,0) ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present - include the "ghost" cell on the right
        ++hiBounds_[0];
        --destInc_[1];
      }
    }

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,0); }
    IntVec high() const{ return hiBounds_; }

  private:
    IntVec hiBounds_, srcInc_, destInc_;
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
        hiBounds_( wdest.extent() - IntVec(0,1,0) ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, wdest.extent(0) )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        ++hiBounds_[1];
        destInc_[2] -= wdest.extent(0);
      }
    }

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,0); }
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
        hiBounds_( wdest.extent() - IntVec(0,0,1) )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        // physical boundary present
        ++hiBounds_[2];
      }
    }

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return IntVec(0,0,0); }
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
        destInc_( 1, 0, wdest.extent(0) )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    ~Stencil2Helper(){}

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

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return wdest_.extent(0) * wdest_.extent(1); }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0) * wdest_.extent(1); }
    
    IntVec src_increment () const{ return IntVec(1,0,0); }
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
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    ~Stencil2Helper(){}

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

    ~Stencil2Helper(){}

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

    ~Stencil2Helper(){}

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

    ~Stencil2Helper(){}

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
   *  \brief Specialization for SVol->XVol (density)
   */
  template<>
  struct Stencil2Helper< SVolField, XVolField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        --hiBounds_[0];
        ++destInc_[1];
      }
    }

    ~Stencil2Helper(){}

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
   *  \brief Specialization for SVol->YVol (density)
   */
  template<>
  struct Stencil2Helper< SVolField, YVolField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, wdest.extent(0) ),
        hiBounds_( wdest.extent() )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        --hiBounds_[1];
        destInc_[2] += wdest.extent(0);
      }
    }

    ~Stencil2Helper(){}

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
   *  \brief Specialization for SVol->ZVol (density)
   */
  template<>
  struct Stencil2Helper< SVolField, ZVolField >
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

    ~Stencil2Helper(){}

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

  //==================================================================


  /**
   *  \brief Specialization for XVol->SVol (pressure solve)
   *
   *  Note that this is the same as the SSurfX->SVol helper which has
   *  already been defined.
   */
  template<>
  struct Stencil2Helper< XVolField, SVolField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : hiBounds_( wdest.extent() - IntVec(1,0,0) ),
        srcInc_ ( 1, 1, 0 ),
        destInc_( 1, 1, 0 )
    {
      if( wsrc.extent(0) != wdest.extent(0) ){
        // physical boundary present
        ++hiBounds_[0];
        --destInc_[1];
      }
    }

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return 1; }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,0); }
    IntVec high() const{ return hiBounds_; }
  private:
    IntVec hiBounds_, srcInc_, destInc_;
  };

  /**
   *  \brief Specialization for YVol->SVol
   *
   *  Note that this is the same as the SSurfY->SVol helper which has
   *  already been defined.
   */
  template<>
  struct Stencil2Helper< YVolField, SVolField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ),  wdest_( wdest ),
        hiBounds_( wdest.extent() - IntVec(0,1,0) ),
        srcInc_ ( 1, 0, wsrc.extent(0) ),
        destInc_( 1, 0, wdest.extent(0) )
    {
      if( wsrc.extent(1) != wdest.extent(1) ){
        // physical boundary present
        ++hiBounds_[1];
        destInc_[2] -= wdest.extent(0);
      }
    }

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0); }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return IntVec(0,0,0); }
    IntVec high() const{ return hiBounds_; }

  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_, srcInc_, destInc_;
  };


  /**
   *  \brief Specialization for ZVol->SVol
   *
   *  Note that this is the same as the SSurfZ->SVol helper which has
   *  already been defined.
   */
  template<>
  struct Stencil2Helper< ZVolField, SVolField >
  {
    Stencil2Helper( const MemoryWindow& wsrc, const MemoryWindow& wdest )
      : wsrc_( wsrc ), wdest_( wdest ),
        hiBounds_( wdest.extent() - IntVec(0,0,1) )
    {
      if( wsrc.extent(2) != wdest.extent(2) ){
        // physical boundary present
        ++hiBounds_[2];
      }
    }

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return 0; }
    unsigned int src_offset_lo() const{ return 0; }
    unsigned int src_offset_hi() const{ return wsrc_.extent(0)*wsrc_.extent(1); }
    
    IntVec src_increment () const{ return IntVec( 1, 0, 0 ); }
    IntVec dest_increment() const{ return IntVec( 1, 0, 0 ); }

    IntVec low () const{ return IntVec(0,0,0); }
    IntVec high() const{ return hiBounds_; }

  private:
    const MemoryWindow &wsrc_, &wdest_;
    IntVec hiBounds_;
  };

  //==================================================================

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
