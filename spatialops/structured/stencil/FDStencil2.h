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

#ifndef FDSTENCIL2_H_
#define FDSTENCIL2_H_

#include <spatialops/structured/IndexTriplet.h>
#include <spatialops/FieldExpressions.h>

namespace SpatialOps{
  namespace structured{

    /**
     *  \author James C. Sutherland
     *  \class FDStencil2
     *  \brief provides finite different stencil for interpolants and
     *         derivatives at second order on a uniform mesh.
     */
    template< typename OpT, typename FieldT, typename DirT >
    class FDStencil2{
      const double coefLo_, coefHi_;
      const NeboStencilCoefList<2> coefList_;

    public:
      typedef OpT       type;           ///< Operation achieved by this operator (e.g. grad, interp, etc.)
      typedef FieldT    SrcFieldType;   ///< Type of the source field
      typedef FieldT    DestFieldType;  ///< Type of the destination field
      typedef DirT      Dir;            ///< Direction that this operator acts in

      typedef typename UnitTriplet<DirT>::type  DirUnitVec;

      typedef typename DirUnitVec::Negate LowStPt;  ///< The offset for the low point
      typedef DirUnitVec                  HighStPt; ///< The offset for the high point

      typedef typename BuildTwoPointList<LowStPt, HighStPt>::Result StPtList;

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

      FDStencil2( const double coefLo, const double coefHi );
      ~FDStencil2();

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

      inline double get_minus_coef() const{ return coefLo_; } ///< get the (-) coefficient
      inline double  get_plus_coef() const{ return coefHi_; } ///< get the (+) coefficient
    };

  }
}


#endif /* FDSTENCIL2_H_ */
