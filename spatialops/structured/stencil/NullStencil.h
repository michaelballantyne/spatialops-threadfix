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

#ifndef SpatialOps_structured_NullStencil_h
#define SpatialOps_structured_NullStencil_h

#include <spatialops/FieldExpressions.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \class NullStencil
   *  \author James C. Sutherland
   *
   *  \brief Direct copy operator.
   *
   *  For some operations, we are simply converting field types for
   *  uniform structured meshes.  Examples include:
   *   - interpolation from the scalar volume to staggered surfaces
   *   - interpolation from staggered volumes to scalar surfaces
   */
  template< typename OperatorT, typename SrcFieldT, typename DestFieldT >
  struct NullStencil
  {
    typedef OperatorT  OpT;
    typedef SrcFieldT  SrcFieldType;
    typedef DestFieldT DestFieldType;

    NullStencil();
    void apply_to_field( const SrcFieldT& src, DestFieldT& dest ) const;

    typedef typename DestFieldType::value_type AtomicType;  // scalar type

    /**
     * \brief Nebo's inline operator for scalar values
     * \param src the scalar to which the operator is applied
     */
    inline AtomicType operator ()( const AtomicType src ) const
    {
        return src;
    }

    // Nebo-related typedefs
    // argument is a field
    typedef NeboConstField<Initial, SrcFieldType> FieldArg;
    typedef NeboExpression<FieldArg, DestFieldType> FieldResult;

    /**
     * \brief Nebo's inline operator for field values
     * \param src the field to which the operator is applied
     */
    inline FieldResult operator ()( const SrcFieldType & src ) const
    {
        return FieldResult(FieldArg(src));
    }

    /**
     * \brief Nebo's inline operator for Nebo expressions
     * \param src the Nebo expression to which the operator is applied
     */
    template<typename Arg>
    inline NeboExpression<Arg, DestFieldType> operator ()( const NeboExpression<Arg, SrcFieldType> & src ) const
    {
        return NeboExpression<Arg, DestFieldType>(src.expr());
    }
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_NullStencil_h
