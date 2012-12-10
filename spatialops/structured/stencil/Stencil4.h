/*
 * Copyright (c) 2011 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef SpatialOps_Structured_Stencil4_h
#define SpatialOps_Structured_Stencil4_h

#include <spatialops/structured/IndexTriplet.h>
#include <spatialops/FieldExpressions.h>

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
    typedef OpT         type;           ///< The operator type (Interpolant, Gradient, Divergence)
    typedef SrcFieldT   SrcFieldType;   ///< The source field type
    typedef DestFieldT  DestFieldType;  ///< The destination field type
    typedef typename SrcFieldType::Location::Offset
                        SrcOffset;      ///< The source field offset
    typedef typename DestFieldType::Location::Offset
                        DestOffset;     ///< The destination field offset
    typedef structured::IndexTriplet<1, 0, 0>
                        XUnit;          ///< X-axis unit vector
    typedef structured::IndexTriplet<0, 1, 0>
                        YUnit;          ///< Y-axis unit vector
    typedef structured::IndexTriplet<0, 0, 1>
                        ZUnit;          ///< Z-axis unit vector
    typedef typename TemplateIf<((int)(SrcOffset::X) != (int)(DestOffset::X)),
                                XUnit,
                                YUnit>::result
                        FirstDir;       ///< The first direction (unit vector)
    typedef typename TemplateIf<((int)(SrcOffset::Z) != (int)(DestOffset::Z)),
                                ZUnit,
                                YUnit>::result
                        SecondDir;      ///< The second direction (unit vector)
    typedef typename structured::Multiply<SrcOffset, FirstDir>::result
                        SrcInFirstDir;  ///< The source offset in the first direction
    typedef typename structured::Multiply<SrcOffset, SecondDir>::result
                        SrcInSecondDir; ///< The source offset in the second direction
    typedef typename structured::Multiply<DestOffset, FirstDir>::result
                        DestInFirstDir; ///< The destination offset in the first direction
    typedef typename structured::Multiply<DestOffset, SecondDir>::result
                        DestInSecondDir; ///< The destination offset in the second direction
    typedef typename structured::GreaterThan<SrcInFirstDir, DestInFirstDir>::result::Negate
                        LoValInFirstDir; ///< The low value in the first direction
    typedef typename structured::LessThan<SrcInFirstDir, DestInFirstDir>::result
                        HiValInFirstDir; ///< The high value in the first direction
    typedef typename structured::GreaterThan<SrcInSecondDir, DestInSecondDir>::result::Negate
                        LoValInSecondDir; ///< The low value in the second direction
    typedef typename structured::LessThan<SrcInSecondDir, DestInSecondDir>::result
                        HiValInSecondDir; ///< The high value in the second direction
    typedef typename structured::Add<LoValInFirstDir, LoValInSecondDir>::result
                        StPt1;          ///< First stencil point location
                                        ///  (relative to the destination point)
    typedef typename structured::Add<HiValInFirstDir, LoValInSecondDir>::result
                        StPt2;          ///< Second stencil point location
                                        ///  (relative to the destination point)
    typedef typename structured::Add<LoValInFirstDir, HiValInSecondDir>::result
                        StPt3;          ///< Third stencil point location
                                        ///  (relative to the destination point)
    typedef typename structured::Add<HiValInFirstDir, HiValInSecondDir>::result
                        StPt4;          ///< Fourth stencil point location
                                        ///  (relative to the destination point)
    typedef typename BuildFourPointList<StPt1, StPt2, StPt3, StPt4>::Result
                        StPtList;       ///< The list of all stencil points in this stencil

    typedef typename DestFieldType::value_type AtomicType;  // scalar type

    // Nebo-related internal struct
    // argument is a Nebo expression
    template<typename Arg>
    struct ResultConstructor {
        typedef NeboStencil<Initial, StPtList, Arg, DestFieldType> Stencil;
        typedef NeboExpression<Stencil, DestFieldType> Result;
    };

    // Nebo-related internal struct
    // argument is a field
    typedef NeboConstField<Initial, SrcFieldType> FieldArg;
    typedef NeboStencil<Initial, StPtList, FieldArg, DestFieldType> FieldStencil;
    typedef NeboExpression<FieldStencil, DestFieldType> FieldResult;

    Stencil4( const double coef1,
              const double coef2,
              const double coef3,
              const double coef4 );

    void apply_to_field( const SrcFieldT& src, DestFieldT& dest ) const;

    /**
     * \brief Nebo's inline operator for scalar values
     * \param src the scalar to which the operator is applied
     */
    inline AtomicType operator ()( const AtomicType src ) const
    {
        return get_coef1() * src
            + get_coef2() * src
            + get_coef3() * src
            + get_coef4() * src;
    }

    /**
     * \brief Nebo's inline operator for field values
     * \param src the field to which the operator is applied
     */
    inline FieldResult operator ()( const SrcFieldType & src ) const
    {
        return FieldResult(FieldStencil(FieldArg(src), coefList_));
    }

    /**
     * \brief Nebo's inline operator for Nebo expressions
     * \param src the Nebo expression to which the operator is applied
     */
    template<typename Arg>
    inline typename ResultConstructor<Arg>::Result operator ()( const NeboExpression<Arg, SrcFieldType> & src ) const
    {
        typedef typename ResultConstructor<Arg>::Stencil Stencil;
        typedef typename ResultConstructor<Arg>::Result Result;
        return Result(Stencil(src.expr(), coefList_));
    }

    inline double get_coef1() const{ return coef1_; } ///< get the first coefficient
    inline double get_coef2() const{ return coef2_; } ///< get the second coefficient
    inline double get_coef3() const{ return coef3_; } ///< get the third coefficient
    inline double get_coef4() const{ return coef4_; } ///< get the fourth coefficient

  private:
    const double coef1_, coef2_, coef3_, coef4_;
    const NeboStencilCoefList<4> coefList_;
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
      typedef typename Multiply< typename UnitTriplet<XDIR>::type,
                                 Difference >::result::Abs              ActiveX;     // {1,0,0} if active
      typedef typename Multiply< typename UnitTriplet<YDIR>::type,
                                 Difference >::result::Abs              ActiveY;     // {0,1,0} if active
      typedef typename Multiply< typename UnitTriplet<ZDIR>::type,
                                 Difference >::result::Abs              ActiveZ;     // {0,0,1} if active

    public:
      typedef typename SelectNonzeroTriplet<ActiveX,ActiveY,ActiveZ>::type1  Dir1Vec;  ///< the first active direction UnitTriplet
      typedef typename SelectNonzeroTriplet<ActiveX,ActiveY,ActiveZ>::type2  Dir2Vec;  ///< the second active direction UnitTriplet

      typedef typename GetNonzeroDir< Dir1Vec >::DirT                        Dir1;     ///< the first active direction
      typedef typename GetNonzeroDir< Dir2Vec >::DirT                        Dir2;     ///< the second active direction
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

      typedef typename Multiply<SFBCExtra,Dir1Vec>::result  UpperLoopBCAug1; ///< shift for dir1 upper bounds when BC is present
      typedef typename Multiply<SFBCExtra,Dir2Vec>::result  UpperLoopBCAug2; ///< shift for dir2 upper bounds when BC is present

      typedef IndexTriplet<0,0,0>                           Src1Offset;      ///< offset for the first source field
      typedef typename Add<Dir1Vec,Dir2Vec>::result::Negate Src1Extent;      ///< extent modification for the first source field
      typedef typename Subtract<
          SFBCExtra,
          typename Add<Dir1Vec,Dir2Vec>::result
          >::result::PositiveOrZero::Negate                 Src1ExtentBC;    ///< amount to augment source1 extent by if a BC is present

      typedef typename Add<Dir1Vec,Src1Offset>::result      Src2Offset;      ///< offset for the second source field
      typedef Src1Extent                                    Src2Extent;      ///< extent modification for the second source field
      typedef Src1ExtentBC                                  Src2ExtentBC;    ///< amount to augment source2 extent by if a BC is present

      typedef typename Add<Dir2Vec,Src1Offset>::result      Src3Offset;      ///< offset for the third source field
      typedef Src1Extent                                    Src3Extent;      ///< extent modification for the third source field
      typedef Src1ExtentBC                                  Src3ExtentBC;    ///< amount to augment source3 extent by if a BC is present

      typedef typename Add<Dir1Vec,Src3Offset>::result      Src4Offset;      ///< offset for the fourth source field
      typedef Src1Extent                                    Src4Extent;      ///< extent modification for the fourth source field
      typedef Src3ExtentBC                                  Src4ExtentBC;    ///< amount to augment source4 extent by if a BC is present

      typedef typename Add<
          typename Multiply< Dir1Vec,
            typename Subtract<
              typename Multiply<Dir1Vec,SFO>::result,
              typename Multiply<Dir1Vec,DFO>::result
              >::result
            >::result::PositiveOrZero,
            typename Multiply< Dir2Vec,
              typename Subtract<
                typename Multiply<Dir2Vec,SFO>::result,
                typename Multiply<Dir2Vec,DFO>::result
                >::result
              >::result::PositiveOrZero
          >::result                                         DestOffset;      ///< the offset for the destination field
      typedef Src1Extent                                    DestExtent;      ///< the extent for the destination field
      typedef typename Subtract<
          typename Multiply<DFO,DFBCExtra>::result,
          typename Multiply< typename Add<Dir1Vec,Dir2Vec>::result,
                             typename Multiply<SFO,SFBCExtra>::result
                             >::result
          >::result                                         DestExtentBC;

    };

  } // namespace s4detail


} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil4_h
