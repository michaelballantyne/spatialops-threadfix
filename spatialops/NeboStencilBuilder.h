/*
 * Copyright (c) 2013 The University of Utah
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

#ifndef NEBO_STENCIL_BUILDER_H
#define NEBO_STENCIL_BUILDER_H

#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps {

  /**
   * \struct NeboStencilBuilder
   * \brief Supports definition of new Nebo stencils.
   *
   * \tparam OperatorType the type of operator, e.g. SpatialOps::Gradient, SpatialOps::Interpolant, etc.
   * \tparam PntCltnT     defines the stencil points
   * \tparam SrcFieldT    the type of field that this operator acts on
   * \tparam DestFieldT   the type of field that this operator produces
   */
    template<typename OperatorType, typename PntCltnT, typename SrcFieldT, typename DestFieldT>
    struct NeboStencilBuilder {
    public:
        typedef OperatorType  type;                 ///< operator type (Interpolant, Gradient, Divergence)
        typedef PntCltnT      PointCollectionType;  ///< collection of stencil points
        typedef SrcFieldT     SrcFieldType;         ///< source field type
        typedef DestFieldT    DestFieldType;        ///< destination field type
        typedef NeboStencilCoefCollection<PointCollectionType::length>
                              CoefCollection;       ///< collection of coefficients

        // typedefs for when argument is a Nebo expression
        template<typename Arg>
        struct WithArg {
            typedef NeboStencil<Initial, PointCollectionType, Arg, DestFieldType> Stencil;
            typedef NeboExpression<Stencil, DestFieldType> Result;
        };

        // typedefs for when argument is a field
        typedef NeboConstField<Initial, SrcFieldType> FieldArg;
        typedef NeboStencil<Initial, PointCollectionType, FieldArg, DestFieldType> FieldStencil;
        typedef NeboExpression<FieldStencil, DestFieldType> FieldResult;

        /**
         *  \brief construct a stencil with the specified coefficients
         */
        NeboStencilBuilder(const CoefCollection & coefs)
         : coefCollection_(coefs)
        {}

        ~NeboStencilBuilder() {}

        /**
         * \brief Return coefficient collection
         */
        CoefCollection const & coefs(void) const { return coefCollection_; };

        /**
         * \brief Apply this operator to the supplied source field to produce the supplied destination field
         * \param src the field that the operator is applied to
         * \param dest the resulting field.
         */
        void apply_to_field( const SrcFieldType & src, DestFieldType & dest ) const {
            dest <<= operator()(src);
        }

        /**
         * \brief Nebo's inline operator for field values
         * \param src the field to which the operator is applied
         */
        inline FieldResult operator ()( const SrcFieldType & src ) const {
            return FieldResult(FieldStencil(FieldArg(src), coefs()));
        }

        /**
         * \brief Nebo's inline operator for Nebo expressions
         * \param src the Nebo expression to which the operator is applied
         */
        template<typename Arg>
        inline typename WithArg<Arg>::Result
        operator ()( const NeboExpression<Arg, SrcFieldType> & src ) const {
            typedef typename WithArg<Arg>::Stencil Stencil;
            typedef typename WithArg<Arg>::Result Result;
            return Result(Stencil(src.expr(), coefs()));
        }

    private:
        const CoefCollection coefCollection_;
    };

    template<typename OperatorType, typename SrcFieldType, typename DestFieldType>
    struct Stencil2Collection {
        // source field offset
        typedef typename SrcFieldType::Location::Offset                      SrcOffset;
        // destination field offset
        typedef typename DestFieldType::Location::Offset                     DestOffset;
        // low (first) stencil point location (relative to the destination point)
        typedef typename structured::GreaterThan<SrcOffset,
                                                 DestOffset>::result::Negate LowStPt;
        // high (second) stencil point location (relative to the destination point)
        typedef typename structured::LessThan<SrcOffset, DestOffset>::result HighStPt;
        // collection of all stencil points in this stencil
        typedef typename BuildTwoPointCollection<LowStPt, HighStPt>::Result  StPtCollection;
    };

    template<typename OperatorType, typename SrcFieldType, typename DestFieldType>
    struct NullStencilCollection {
        typedef structured::IndexTriplet<0, 0, 0> StencilPoint;
        typedef BuildStencilPointCollection<StencilPoint>::Result StPtCollection;
    };

    template<typename OperatorType, typename SrcFieldType, typename DestFieldType>
    struct Stencil4Collection {
        template<bool Boolean, typename True, typename False>
        struct TemplateIf;

        template<typename True, typename False>
        struct TemplateIf<true, True, False> { True typedef result; };

        template<typename True, typename False>
        struct TemplateIf<false, True, False> { False typedef result; };

        // source field offset
        typedef typename SrcFieldType::Location::Offset                             SrcOffset;
        // destination field offset
        typedef typename DestFieldType::Location::Offset                            DestOffset;
        // unit vectors
        typedef structured::IndexTriplet<1, 0, 0>                                   XUnit;
        typedef structured::IndexTriplet<0, 1, 0>                                   YUnit;
        typedef structured::IndexTriplet<0, 0, 1>                                   ZUnit;
        // first direction (unit vector)
        typedef typename TemplateIf<((int)(SrcOffset::X) != (int)(DestOffset::X)),
                                    XUnit,
                                    YUnit>::result                                  FirstDir;
        // second direction (unit vector)
        typedef typename TemplateIf<((int)(SrcOffset::Z) != (int)(DestOffset::Z)),
                                    ZUnit,
                                    YUnit>::result                                  SecondDir;
        // source offset in the first direction
        typedef typename structured::Multiply<SrcOffset, FirstDir>::result          SrcInFirstDir;
        // source offset in the second direction
        typedef typename structured::Multiply<SrcOffset, SecondDir>::result         SrcInSecondDir;
        // destination offset in the first direction
        typedef typename structured::Multiply<DestOffset, FirstDir>::result         DestInFirstDir;
        // destination offset in the second direction
        typedef typename structured::Multiply<DestOffset, SecondDir>::result        DestInSecondDir;
        // low value in the first direction
        typedef typename structured::GreaterThan<SrcInFirstDir,
                                                 DestInFirstDir>::result::Negate    LoValInFirstDir;
        // high value in the first direction
        typedef typename structured::LessThan<SrcInFirstDir,
                                              DestInFirstDir>::result               HiValInFirstDir;
        // low value in the second direction
        typedef typename structured::GreaterThan<SrcInSecondDir,
                                            DestInSecondDir>::result::Negate        LoValInSecondDir;
        // high value in the second direction
        typedef typename structured::LessThan<SrcInSecondDir,
                                              DestInSecondDir>::result              HiValInSecondDir;
        // stencil point locations (relative to the destination point)
        typedef typename structured::Add<LoValInFirstDir, LoValInSecondDir>::result StPt1;
        typedef typename structured::Add<HiValInFirstDir, LoValInSecondDir>::result StPt2;
        typedef typename structured::Add<LoValInFirstDir, HiValInSecondDir>::result StPt3;
        typedef typename structured::Add<HiValInFirstDir, HiValInSecondDir>::result StPt4;
        // collection of all stencil points in this stencil
        typedef typename BuildFourPointCollection<StPt1,
                                                  StPt2,
                                                  StPt3,
                                                  StPt4>::Result                    StPtCollection;
    };

    template<typename OperatorType, typename SrcFieldType, typename DestFieldType>
    struct FDStencilCollection {
        typedef typename OperatorType::DirT                                 DirT;
        typedef typename structured::UnitTriplet<DirT>::type                DirVec;
        typedef typename DirVec::Negate                                     LowStPt;
        typedef DirVec                                                      HighStPt;
        typedef typename BuildTwoPointCollection<LowStPt, HighStPt>::Result StPtCollection;
    };

} // namespace SpatialOps
#endif // NEBO_STENCIL_BUILDER_H
