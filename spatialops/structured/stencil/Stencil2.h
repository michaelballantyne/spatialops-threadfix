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

#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

#include <spatialops/structured/IndexTriplet.h>
#include <spatialops/Nebo.h>

namespace SpatialOps {
  namespace structured {

    /**
     *  \class Stencil2
     *  \author James C. Sutherland
     *
     *  \brief Support for implementing simple two-point stencils in
     *         one-dimension on structured meshes.
     *
     *  \tparam OpT - the type of operator (Interpolant, Gradient, Divergence)
     *  \tparam SrcT - the type of field the operator is applied to
     *  \tparam DestT - the type of field the operator produces
     */
    template<typename OperatorT, typename SrcFieldT, typename DestFieldT>
    class Stencil2
    {
      const double coefLo_, coefHi_;
      const NeboStencilCoefCollection<2> coefCollection_;

    public:

      typedef OperatorT   type;           ///< The operator type (Interpolant, Gradient, Divergence)
      typedef SrcFieldT   SrcFieldType;   ///< The source field type
      typedef DestFieldT  DestFieldType;  ///< The destination field type
      typedef typename SrcFieldType::Location::Offset
                          SrcOffset;      ///< The source field offset
      typedef typename DestFieldType::Location::Offset
                          DestOffset;     ///< The destination field offset
      typedef typename structured::GreaterThan<SrcOffset, DestOffset>::result::Negate
                          LowStPt;        ///< The low (first) stencil point location
                                          ///  (relative to the destination point)
      typedef typename structured::LessThan<SrcOffset, DestOffset>::result
                          HighStPt;       ///< The high (second) stencil point location
                                          ///  (relative to the destination point)
      typedef typename BuildTwoPointCollection<LowStPt, HighStPt>::Result
                          StPtCollection;       ///< The collection of all stencil points in this stencil

      typedef typename DestFieldType::value_type  AtomicType;  // scalar type

      // Nebo-related internal struct
      // argument is a Nebo expression
      template<typename Arg>
      struct ResultConstructor {
          typedef NeboStencil<Initial, StPtCollection, Arg, DestFieldType> Stencil;
          typedef NeboExpression<Stencil, DestFieldType> Result;
      };

      // Nebo-related internal struct
      // argument is a field
      typedef NeboConstField<Initial, SrcFieldType> FieldArg;
      typedef NeboStencil<Initial, StPtCollection, FieldArg, DestFieldType> FieldStencil;
      typedef NeboExpression<FieldStencil, DestFieldType> FieldResult;

      /**
       *  \brief construct a stencil with the specified coefficients
       *  \param coefLo the coefficient to multiply the (-) side field by
       *  \param coefHi the coefficient to multiply the (+) side field by
       */
      Stencil2( const double coefLo, const double coefHi );

      ~Stencil2();

      /**
       * \brief Apply this operator to the supplied source field to produce the supplied destination field
       * \param src the field that the operator is applied to
       * \param dest the resulting field.
       */
      void apply_to_field( const SrcFieldType& src, DestFieldType& dest ) const;

      /**
       * \brief Nebo's inline operator for scalar values
       * \param src the scalar to which the operator is applied
       */
      inline AtomicType operator ()( const AtomicType src ) const
      {
          return get_minus_coef() * src + get_plus_coef() * src;
      }

      /**
       * \brief Nebo's inline operator for field values
       * \param src the field to which the operator is applied
       */
      inline FieldResult operator ()( const SrcFieldType & src ) const
      {
          return FieldResult(FieldStencil(FieldArg(src), coefCollection_));
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
          return Result(Stencil(src.expr(), coefCollection_));
      }

      inline double get_minus_coef() const{ return coefLo_; } ///< get the (-) coefficient
      inline double  get_plus_coef() const{ return coefHi_; } ///< get the (+) coefficient
    };

    /*******************************************************************
     *
     * NOTE: all information in the s2detail namespace is meant only for
     *       use within the Stencil2 class and should not be used elsewhere!
     *
     ******************************************************************/
    namespace s2detail {

      /**
       *  \struct ActiveDir
       *
       *  \brief determines the direction that the operator is acting along
       *
       *  \tparam SrcT the source field type for this 2-point stencil operator
       *  \tparam DestT the destination field type for this 2-point stencil operator
       *
       *  Example:
       *  \code typedef ActiveDir<SrcT,DestT>::type OpDir; \endcode
       */
      template<typename SrcT, typename DestT>
      struct ActiveDir
      {
      private:
        typedef typename Subtract< typename  SrcT::Location::Offset,
                                   typename DestT::Location::Offset >::result  Difference;
      public:
        typedef typename GetNonzeroDir< Difference >::DirT type; ///< the direction that the operator acts in.
      };

      template<typename T> struct ActiveDir<T,T>;  ///< invalid - all Stencil2 must do something.


      /**
       * \struct ExtentsAndOffsets
       * \author James C. Sutherland
       * \brief Provides typedefs for dealing with extents and offsets for Stencil2 operators.
       */
      template<typename SrcT, typename DestT>
      struct ExtentsAndOffsets
      {
      private:
        typedef typename  SrcT::Location::Offset                SFO;            ///< Offset information for Source field
        typedef typename  SrcT::Location::BCExtra               SFBCExtra;      ///< Extra cell information for Source field
        typedef typename DestT::Location::Offset                DFO;            ///< Offset information for Destination field
        typedef typename DestT::Location::BCExtra               DFBCExtra;      ///< Extra cell information for Destination field

      public:

        typedef typename ActiveDir<SrcT,DestT>::type            Dir;            ///< The direction that this operator acts in

        typedef typename UnitTriplet<Dir>::type                 DirUnitVec;     ///< unit vector for the direction that this operator acts in.

        typedef typename Subtract<
            typename Multiply<SFBCExtra,DFBCExtra>::result,
            DirUnitVec>::result::PositiveOrZero                 BCMatchNotOnDir;///< amount to augment for BCExtra along axis that is not Dir

        typedef IndexTriplet<0,0,0>                             Src1Offset;     ///< offset for the first source field
        typedef typename DirUnitVec::Negate                     Src1Extent;     ///< extent modification for the first source field
        typedef typename Add<BCMatchNotOnDir,
            typename Subtract<SFBCExtra,
                DirUnitVec>::result::PositiveOrZero::Negate
            >::result                                           Src1ExtentBC;   ///< amount to augment source extents by if a BC is present

        typedef typename Add<DirUnitVec,Src1Offset>::result     Src2Offset;     ///< offset for the second source field
        typedef Src1Extent                                      Src2Extent;     ///< extent modification for the second source field
        typedef Src1ExtentBC                                    Src2ExtentBC;   ///< additional extent modification if a BC is present

        typedef typename Multiply< DirUnitVec,
            typename Subtract<
              typename Multiply<DirUnitVec,SFO>::result,
              typename Multiply<DirUnitVec,DFO>::result
              >::result
            >::result::PositiveOrZero                           DestOffset;     ///< the offset for the destination field
        typedef Src1Extent                                      DestExtent;     ///< the extent for the destination field
        typedef typename Add<BCMatchNotOnDir,
            typename Subtract<
                typename Multiply<DFO,DFBCExtra>::result,
                typename Multiply< DirUnitVec,
                    typename Multiply<SFO,SFBCExtra>::result
                    >::result
                >::result
            >::result                                           DestExtentBC;   ///< amount to augment destination extents by if a BC is present
      };

    } // namespace s2detail

  }// namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
