#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

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


  /*******************************************************************
   *
   * NOTE: all information in the s2detail namespace is meant only for
   *       use within the Stencil2Helper class and should not be used
   *       elsewhere!
   *
   ******************************************************************/
  namespace s2detail{

    /**
     *  \struct ActiveDir
     *
     *  \brief determines the direction that the operator is acting along
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     *
     *  Note that only specialized versions of this struct exist.
     */
    template< typename SrcT, typename DestT > struct ActiveDir;

    template< typename T > struct ActiveDir<T,T>; // invalid - all Stencil2 must do something.

    template< typename VolT > struct ActiveDir< VolT, typename FaceTypes<VolT>::XFace >{ typedef XDIR type; };
    template< typename VolT > struct ActiveDir< VolT, typename FaceTypes<VolT>::YFace >{ typedef YDIR type; };
    template< typename VolT > struct ActiveDir< VolT, typename FaceTypes<VolT>::ZFace >{ typedef ZDIR type; };

    template< typename SurfT > struct ActiveDir<SurfT,typename VolType<SurfT>::VolField>{ typedef typename SurfT::Location::FaceDir type; };

    template<> struct ActiveDir<XVolField,YSurfXField>{ typedef YDIR type; };
    template<> struct ActiveDir<XVolField,ZSurfXField>{ typedef ZDIR type; };

    template<> struct ActiveDir<YVolField,XSurfYField>{ typedef XDIR type; };
    template<> struct ActiveDir<YVolField,ZSurfYField>{ typedef XDIR type; };

    template<> struct ActiveDir<ZVolField,XSurfYField>{ typedef XDIR type; };
    template<> struct ActiveDir<ZVolField,YSurfYField>{ typedef YDIR type; };


    /**
     *  \struct HasPlusBoundary
     *
     *  \brief determine if there is a boundary present in the (+)
     *         direction of the application of this operator.
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     *
     *  The default implementation functions properly for most
     *  operators.  There are specializations for several operators,
     *  however.
     */
    template< typename SrcT, typename DestT >
    struct HasPlusBoundary{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        const size_t dir = ActiveDir<SrcT,DestT>::type::value;
        return ( src.extent(dir) != dest.extent(dir) );
      }
    };

    /**
     *  \struct SrcIncrement
     *
     *  \brief Obtain the increment counters when an iterator bump occurs on the source field.
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     */
    template< typename SrcT, typename DestT >
    struct SrcIncrement
    {
      typedef typename s2detail::ActiveDir<SrcT,DestT>::type ActiveDir;
      static IntVec value( const MemoryWindow& mw ){
        return IntVec( stride<XDIR>(mw),
                       CompareDirection<ActiveDir,XDIR>::same() ? stride<XDIR>(mw) : 0,
                       CompareDirection<ActiveDir,YDIR>::same() ? stride<YDIR>(mw) : 0 );
      }
    };

    /**
     *  \struct DestIncrement
     *
     *  \brief Obtain the increment counters when an iterator bump occurs on the destination field.
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     */
    template< typename SrcT, typename DestT >
    struct DestIncrement
    {
      typedef typename s2detail::ActiveDir<SrcT,DestT>::type ActiveDir;
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec inc( SrcIncrement<SrcT,DestT>::value(mw) );
        if( hasPlusBoundary ){
          const int sign = IndexStagger<SrcT,ActiveDir>::value-IndexStagger<DestT,ActiveDir>::value;
          switch( ActiveDir::value ){
          case 0: inc[1] += sign; break;
          case 1: inc[2] += sign*mw.extent(0); break;
          case 2: break;
          }
        }
        return inc;
      }
    };

    /**
     *  \struct SetLowBounds
     *
     *  \brief Obtain the starting (lower) index for looping through the fields.
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     */
    template< typename SrcT, typename DestT >
    struct SetLowBounds{
      static IntVec value(){
        // either 0 or 1 in each entry, corresponding to the active
        // direction and the staggering of the src/dest fields.
        return IntVec( std::max( IndexStagger<SrcT,XDIR>::value-IndexStagger<DestT,XDIR>::value, 0 ),
                       std::max( IndexStagger<SrcT,YDIR>::value-IndexStagger<DestT,YDIR>::value, 0 ),
                       std::max( IndexStagger<SrcT,ZDIR>::value-IndexStagger<DestT,ZDIR>::value, 0 ) );
      }
    };

    /**
     *  \struct SetLowBounds
     *
     *  \brief Obtain the ending (upper) index for looping through the fields.
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     */
    template< typename SrcT, typename DestT >
    struct SetHiBounds
    {
      typedef typename s2detail::ActiveDir<SrcT,DestT>::type ActiveDir;
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec vec( mw.extent() );
        vec[ ActiveDir::value ] -= std::max( IndexStagger<DestT,ActiveDir>::value-IndexStagger<SrcT,ActiveDir>::value, 0 );
        if( hasPlusBoundary ){
          const int sign = IndexStagger<DestT,ActiveDir>::value-IndexStagger<SrcT,ActiveDir>::value;
          vec[ ActiveDir::value ] += sign;
        }
        return vec;
      }
    };

    /**
     *  \struct InitialOffsets
     *
     *  \brief Determines the number of bumps to apply to source and
     *         destination iterators to set up for the stencil application.
     *
     *  \tparam SrcT the source field type for this 2-point stencil operator
     *  \tparam DestT the destination field type for this 2-point stencil operator
     */
    template< typename SrcT, typename DestT >
    struct InitialOffsets{
      typedef typename s2detail::ActiveDir<SrcT,DestT>::type ActiveDir;
      static unsigned int dest( const MemoryWindow& mw ){
        const int tmp = std::max( IndexStagger<SrcT,ActiveDir>::value-IndexStagger<DestT,ActiveDir>::value, 0 );
        switch( ActiveDir::value ){
        case XDIR::value: 
          return tmp;
          break;
        case YDIR::value:
          return tmp * mw.extent(0);
          break;
        case ZDIR::value:
          return tmp * mw.extent(0)*mw.extent(1);
          break;
        }
        assert( 0 );
      }
      static unsigned int src_lo( const MemoryWindow& mw ){
        return 0;
      }
      static unsigned int src_hi( const MemoryWindow& mw ){
        switch( ActiveDir::value ){
        case XDIR::value: return 1;
        case YDIR::value: return mw.extent(0);
        case ZDIR::value: return mw.extent(0)*mw.extent(1);
        } 
        assert( 0 );
      }
    };


    // specializations for XVolField -> YSurfXField (advecting velocity)
    template<> struct SrcIncrement<XVolField,YSurfXField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, mw.extent(0) ); }
    };
    template<> struct InitialOffsets<XVolField,YSurfXField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return mw.extent(0); }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0); }
    };
    template<> struct DestIncrement<XVolField,YSurfXField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, hasPlusBoundary ? 1 : 0, mw.extent(0) );
      }
    };
    template<> struct SetHiBounds<XVolField,YSurfXField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[0];
        return bounds;
      }
    };
    template<> struct SetLowBounds<XVolField,YSurfXField>{
      static IntVec value(){ return IntVec(0,1,0); }
    };
    template<> struct HasPlusBoundary<XVolField,YSurfXField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(0) != dest.extent(0) );
      }
    };


    // specializations for XVolField -> ZSurfXField (advecting velocity)
    template<> struct SrcIncrement<XVolField,ZSurfXField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, 0 ); }
    };
    template<> struct InitialOffsets<XVolField,ZSurfXField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
    };
    template<> struct DestIncrement<XVolField,ZSurfXField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, hasPlusBoundary ? 1 : 0, 0 );
      }
    };
    template<> struct SetHiBounds<XVolField,ZSurfXField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[0];
        return bounds;
      }
    };
    template<> struct SetLowBounds<XVolField,ZSurfXField>{
      static IntVec value(){ return IntVec(0,0,1); }
    };
    template<> struct HasPlusBoundary<XVolField,ZSurfXField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(0) != dest.extent(0) );
      }
    };


    // specializations for YVolField -> XSurfYField (advecting velocity)
    template<> struct SrcIncrement<YVolField,XSurfYField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 1, 0 ); }
    };
    template<> struct DestIncrement<YVolField,XSurfYField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 1, hasPlusBoundary ? mw.extent(0) : 0 );
      }
    };
    template<> struct InitialOffsets<YVolField,XSurfYField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return 1; }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return 1; }
    };
    template<> struct SetHiBounds<YVolField,XSurfYField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[1];
        return bounds;
      }
    };
    template<> struct SetLowBounds<YVolField,XSurfYField>{
      static IntVec value(){ return IntVec(1,0,0); }
    };
    template<> struct HasPlusBoundary<YVolField,XSurfYField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(1) != dest.extent(1) );
      }
    };


    // specializations for YVolField -> ZSurfYField (advecting velocity)
    template<> struct SrcIncrement<YVolField,ZSurfYField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, 0 ); }
    };
    template<> struct DestIncrement<YVolField,ZSurfYField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 0, hasPlusBoundary ? mw.extent(0) : 0 );
      }
    };
    template<> struct InitialOffsets<YVolField,ZSurfYField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
    };
    template<> struct SetHiBounds<YVolField,ZSurfYField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[1];
        return bounds;
      }
    };
    template<> struct SetLowBounds<YVolField,ZSurfYField>{
      static IntVec value(){ return IntVec(0,0,1); }
    };
    template<> struct HasPlusBoundary<YVolField,ZSurfYField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(1) != dest.extent(1) );
      }
    };


    // specializations for ZVolField -> XSurfZField (advecting velocity)
    template<> struct SrcIncrement<ZVolField,XSurfZField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 1, 0 ); }
    };
    template<> struct DestIncrement<ZVolField,XSurfZField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 1, 0 );
      }
    };
    template<> struct InitialOffsets<ZVolField,XSurfZField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return 1; }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return 1; }
    };
    template<> struct SetHiBounds<ZVolField,XSurfZField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[2];
        return bounds;
      }
    };
    template<> struct SetLowBounds<ZVolField,XSurfZField>{
      static IntVec value(){ return IntVec(1,0,0); }
    };
    template<> struct HasPlusBoundary<ZVolField,XSurfZField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(2) != dest.extent(2) );
      }
    };


    // specializations for ZVolField -> YSurfZField (advecting velocity)
    template<> struct SrcIncrement<ZVolField,YSurfZField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, mw.extent(0) ); }
    };
    template<> struct DestIncrement<ZVolField,YSurfZField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 0, mw.extent(0) );
      }
    };
    template<> struct InitialOffsets<ZVolField,YSurfZField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return mw.extent(0); }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0); }
    };
    template<> struct SetHiBounds<ZVolField,YSurfZField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[2];
        return bounds;
      }
    };
    template<> struct SetLowBounds<ZVolField,YSurfZField>{
      static IntVec value(){ return IntVec(0,1,0); }
    };
    template<> struct HasPlusBoundary<ZVolField,YSurfZField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(2) != dest.extent(2) );
      }
    };


    // specializations for SVolField -> XVolField (density)
    template<> struct SrcIncrement<SVolField,XVolField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 1, 0 ); }
    };
    template<> struct DestIncrement<SVolField,XVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 1, 0 );
      }
    };
    template<> struct InitialOffsets<SVolField,XVolField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return 1; }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return 1; }
    };
    template<> struct SetHiBounds<SVolField,XVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[0];
        return bounds;
      }
    };
    template<> struct SetLowBounds<SVolField,XVolField>{
      static IntVec value(){ return IntVec(1,0,0); }
    };
    template<> struct HasPlusBoundary<SVolField,XVolField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(0) != dest.extent(0) );
      }
    };


    // specializations for SVolField -> YVolField (density)
    template<> struct SrcIncrement<SVolField,YVolField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, mw.extent(0) ); }
    };
    template<> struct DestIncrement<SVolField,YVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 0, mw.extent(0) );
      }
    };
    template<> struct InitialOffsets<SVolField,YVolField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return mw.extent(0); }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0); }
    };
    template<> struct SetHiBounds<SVolField,YVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[1];
        return bounds;
      }
    };
    template<> struct SetLowBounds<SVolField,YVolField>{
      static IntVec value(){ return IntVec(0,1,0); }
    };
    template<> struct HasPlusBoundary<SVolField,YVolField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(1) != dest.extent(1) );
      }
    };


    // specializations for SVolField -> ZVolField (density)
    template<> struct SrcIncrement<SVolField,ZVolField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, 0 ); }
    };
    template<> struct DestIncrement<SVolField,ZVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 0, 0 );
      }
    };
    template<> struct InitialOffsets<SVolField,ZVolField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
    };
    template<> struct SetHiBounds<SVolField,ZVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() );
        if( hasPlusBoundary ) --bounds[2];
        return bounds;
      }
    };
    template<> struct SetLowBounds<SVolField,ZVolField>{
      static IntVec value(){ return IntVec(0,0,1); }
    };
    template<> struct HasPlusBoundary<SVolField,ZVolField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(2) != dest.extent(2) );
      }
    };


    // specializations for XVolField -> SVolField (pressure projection)
    template<> struct SrcIncrement<XVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 1, 0 ); }
    };
    template<> struct DestIncrement<XVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( hasPlusBoundary ? 0 : 1, 1, 0 );
      }
    };
    template<> struct InitialOffsets<XVolField,SVolField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return 1; }
    };
    template<> struct SetHiBounds<XVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() - IntVec(1,0,0) );
        if( hasPlusBoundary ) ++bounds[0];
        return bounds;
      }
    };
    template<> struct SetLowBounds<XVolField,SVolField>{
      static IntVec value(){ return IntVec(0,0,0); }
    };
    template<> struct HasPlusBoundary<XVolField,SVolField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(0) != dest.extent(0) );
      }
    };


    // specializations for YVolField -> SVolField (pressure projection)
    template<> struct SrcIncrement<YVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, mw.extent(0) ); }
    };
    template<> struct DestIncrement<YVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 0, hasPlusBoundary ? 0 : mw.extent(0) );
      }
    };
    template<> struct InitialOffsets<YVolField,SVolField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0); }
    };
    template<> struct SetHiBounds<YVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() - IntVec(0,1,0) );
        if( hasPlusBoundary ) ++bounds[1];
        return bounds;
      }
    };
    template<> struct SetLowBounds<YVolField,SVolField>{
      static IntVec value(){ return IntVec(0,0,0); }
    };
    template<> struct HasPlusBoundary<YVolField,SVolField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(1) != dest.extent(1) );
      }
    };


    // specializations for ZVolField -> SVolField (pressure projection)
    template<> struct SrcIncrement<ZVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw ){ return IntVec( 1, 0, 0 ); }
    };
    template<> struct DestIncrement<ZVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        return IntVec( 1, 0, 0 );
      }
    };
    template<> struct InitialOffsets<ZVolField,SVolField>{
      static unsigned int dest  ( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_lo( const MemoryWindow& mw ){ return 0; }
      static unsigned int src_hi( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }
    };
    template<> struct SetHiBounds<ZVolField,SVolField>{
      static IntVec value( const MemoryWindow& mw, const bool hasPlusBoundary ){
        IntVec bounds( mw.extent() - IntVec(0,0,1) );
        if( hasPlusBoundary ) ++bounds[2];
        return bounds;
      }
    };
    template<> struct SetLowBounds<ZVolField,SVolField>{
      static IntVec value(){ return IntVec(0,0,0); }
    };
    template<> struct HasPlusBoundary<ZVolField,SVolField>{
      static bool value( const MemoryWindow& src, const MemoryWindow& dest ){
        return ( src.extent(2) != dest.extent(2) );
      }
    };

  } // namespace s2detail


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
  template< typename SrcT, typename DestT >  struct Stencil2Helper;

  template< typename SrcT, typename DestT >
  struct Stencil2Helper
  {
    Stencil2Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest )
      : bcPlus_  ( s2detail::HasPlusBoundary<SrcT,DestT>::value(wsrc,wdest) ),
        doffset_ ( s2detail::InitialOffsets <SrcT,DestT>::dest(wdest)  ),
        soffsetL_( s2detail::InitialOffsets <SrcT,DestT>::src_lo(wsrc) ),
        soffsetH_( s2detail::InitialOffsets <SrcT,DestT>::src_hi(wsrc) ),
        srcInc_  ( s2detail::SrcIncrement   <SrcT,DestT>::value(wsrc ) ),
        destInc_ ( s2detail::DestIncrement  <SrcT,DestT>::value(wdest,bcPlus_ ) ),
        loBounds_( s2detail::SetLowBounds   <SrcT,DestT>::value() ),
        hiBounds_( s2detail::SetHiBounds    <SrcT,DestT>::value(wdest,bcPlus_) )
    {}

    ~Stencil2Helper(){}

    unsigned int dest_offset  () const{ return doffset_ ; }
    unsigned int src_offset_lo() const{ return soffsetL_; }
    unsigned int src_offset_hi() const{ return soffsetH_; }
    
    IntVec src_increment () const{ return srcInc_; }
    IntVec dest_increment() const{ return destInc_; }

    IntVec low () const{ return loBounds_; }
    IntVec high() const{ return hiBounds_; }

  private:
    const int bcPlus_;
    const unsigned int doffset_, soffsetL_, soffsetH_;
    const IntVec srcInc_, destInc_, loBounds_, hiBounds_;
  };


  //==================================================================

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
