#ifndef SpatialOps_Structured_Stencil4_h
#define SpatialOps_Structured_Stencil4_h

#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \struct Stencil4
   *  \author James C. Sutherland
   *
   *  \brief Intended for use on 4-point stencils on structured, uniform meshes
   *
   *  \tparam OpT the type of operator (Interpolant, Gradient)
   *  \tparam SrcFieldT the type of field we apply the operator to
   *  \tparam DestFieldT the type of field produced by the operator
   *
   *  Examples of 4-point stencils:
   *   X-volume to y-surface or z-surface (advecting velocity)
   *   Y-volume to x-surface or z-surface (advecting velocity)
   *   Z-volume to x-surface or y-surface (advecting velocity)
   */
  template< typename OpT, typename SrcFieldT, typename DestFieldT >
  struct Stencil4
  {
    typedef OpT         Type;           ///< The operator type (Interpolant, Gradient, Divergence)
    typedef SrcFieldT   SrcFieldType;   ///< The source field type
    typedef DestFieldT  DestFieldType;  ///< The destination field type


    Stencil4( const double coef1,
              const double coef2,
              const double coef3,
              const double coef4 );

    void apply_to_field( const SrcFieldT& src, DestFieldT& dest ) const;

  private:
    const double coef1_, coef2_, coef3_, coef4_;
  };



  /*******************************************************************
   *
   * NOTE: all information in the s2detail namespace is meant only for
   *       use within the Stencil2Helper class and should not be used
   *       elsewhere!
   *
   ******************************************************************/
  namespace s4detail{

    typedef IndexTriplet<0,0,0>  ZeroTriplet;

    template< typename IT1, typename IT2, typename IT3> struct SelectNonzeroTriplet;

    template< typename IT1, typename IT2 > struct SelectNonzeroTriplet< IT1, IT2, ZeroTriplet >{ typedef IT1 type1;  typedef IT2 type2; };
    template< typename IT1, typename IT3 > struct SelectNonzeroTriplet< IT1, ZeroTriplet, IT3 >{ typedef IT1 type1;  typedef IT3 type2; };
    template< typename IT2, typename IT3 > struct SelectNonzeroTriplet< ZeroTriplet, IT2, IT3 >{ typedef IT2 type1;  typedef IT3 type2; };

    template< typename SrcT, typename DestT >
    struct ActiveDirs
    {
      typedef typename Subtract<
          typename  SrcT::Location::Offset,
          typename DestT::Location::Offset >::result                    Difference;
    private:
      typedef typename Multiply< typename UnitTriplet<XDIR>::type, Difference >::result  ActiveX;     // {1,0,0} if active
      typedef typename Multiply< typename UnitTriplet<YDIR>::type, Difference >::result  ActiveY;     // {0,1,0} if active
      typedef typename Multiply< typename UnitTriplet<ZDIR>::type, Difference >::result  ActiveZ;     // {0,0,1} if active

    public:
      typedef typename SelectNonzeroTriplet<ActiveX,ActiveY,ActiveZ>::type1     Dir1Vec;  ///< the first active direction UnitTriplet
      typedef typename SelectNonzeroTriplet<ActiveX,ActiveY,ActiveZ>::type2     Dir2Vec;  ///< the second active direction UnitTriplet

      typedef typename GetNonzeroDir< Dir1Vec >::DirT                           Dir1;     ///< the first active direction
      typedef typename GetNonzeroDir< Dir2Vec >::DirT                           Dir2;     ///< the second active direction
    };

    template< typename T > struct ActiveDirs<T,T>;  // invalid - all Stencil4 must do something.

    /**
     * \struct ExtentsAndOffsets
     * \author James C. Sutherland
     * \brief  Information about the extents and offsets for Stencil4
     */
    template< typename SrcT, typename DestT >
    struct ExtentsAndOffsets
    {
    private:
      typedef typename  SrcT::Location::Offset              SFO;            ///< Offset information for Source field
      typedef typename  SrcT::Location::BCExtra             SFBCExtra;      ///< Extra cell information for Source field
      typedef typename DestT::Location::Offset              DFO;            ///< Offset information for Destination field
      typedef typename DestT::Location::BCExtra             DFBCExtra;      ///< Extra cell information for Destination field

    public:

      typedef typename ActiveDirs<SrcT,DestT>::Dir1         Dir1;            ///< The first active direction
      typedef typename ActiveDirs<SrcT,DestT>::Dir2         Dir2;            ///< The second active direction

      typedef typename ActiveDirs<SrcT,DestT>::Dir1Vec      Dir1Vec;         ///< The first active direction
      typedef typename ActiveDirs<SrcT,DestT>::Dir2Vec      Dir2Vec;         ///< The second active direction

      typedef typename Add<Dir1Vec,Dir2Vec>::result::Negate UpperLoopShift;  ///< shift for uppper bounds

      typedef typename Multiply<SFBCExtra,Dir1Vec>::result  UpperLoopBCAug1; ///< shift for dir1 upper bounds when BC is present
      typedef typename Multiply<SFBCExtra,Dir2Vec>::result  UpperLoopBCAug2; ///< shift for dir2 upper bounds when BC is present

      typedef IndexTriplet<0,0,0>                           Src1Offset;      ///< offset for the first source field
      typedef UpperLoopShift                                Src1Extent;      ///< extent modification for the first source field
      typedef typename Subtract<SFBCExtra,
          DFBCExtra>::result::PositiveOrZero::Negate        Src1ExtentBC;    ///< amount to augment source1 extent by if a BC is present

      typedef typename Add<Dir1Vec,Src1Offset>::result      Src2Offset;      ///< offset for the second source field
      typedef Src1Extent                                    Src2Extent;      ///< extent modification for the second source field
      typedef Src1ExtentBC                                  Src2ExtentBC;    ///< amount to augment source2 extent by if a BC is present

      typedef typename Add<Dir2Vec,Src1Offset>::result      Src3Offset;      ///< offset for the third source field
      typedef UpperLoopShift                                Src3Extent;      ///< extent modification for the third source field
      typedef typename Subtract<SFBCExtra,
          DFBCExtra>::result::PositiveOrZero::Negate        Src3ExtentBC;    ///< amount to augment source3 extent by if a BC is present

      typedef typename Add<Dir1Vec,Src3Offset>::result      Src4Offset;      ///< offset for the fourth source field
      typedef Src3Extent                                    Src4Extent;      ///< extent modification for the fourth source field
      typedef Src3ExtentBC                                  Src4ExtentBC;    ///< amount to augment source4 extent by if a BC is present

      typedef typename Multiply<
          typename ActiveDirs<SrcT,DestT>::Difference,
          typename Subtract<SFO,DFO>::result
          >::result::PositiveOrZero                         DestOffset;      ///< the offset for the destination field
      typedef UpperLoopShift                                DestExtent;      ///< the extent for the destination field
      typedef typename Subtract<DFBCExtra,
          SFBCExtra>::result::PositiveOrZero::Negate        DestExtentBC;    ///< amount to augment destination extents by if a BC is present

    };

  } // namespace s4detail


} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil4_h
