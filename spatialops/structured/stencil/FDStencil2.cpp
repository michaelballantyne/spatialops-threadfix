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

#include "FDStencil2.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/WriteMatlab.h>
#include <spatialops/FieldExpressionsFDStencil2.h>

namespace SpatialOps{
namespace structured{

template< typename OpT, typename FieldT, typename DirT >
FDStencil2<OpT,FieldT,DirT>::FDStencil2( const double coefLo, const double coefHi )
: coefLo_( coefLo ),
  coefHi_( coefHi )
{}

template< typename OpT, typename FieldT, typename DirT >
FDStencil2<OpT,FieldT,DirT>::~FDStencil2()
{}

template< typename OpT, typename FieldT, typename DirT >
void
FDStencil2<OpT,FieldT,DirT>::apply_to_field( const FieldT& src, FieldT& dest ) const
{
    NeboFDStencilConstructor<FDStencil2<OpT,FieldT,DirT> > typedef Stencil;
    const Stencil stencil(this);
    dest <<= stencil(src);
    //fd_stencil_2_apply_to_field_general_execute<OpT,FieldT,DirVec>(src, dest, coefLo_, coefHi_);
}

// Explicit template instantiation
#define DECLARE_CLASS( OP, FIELDT )     \
  template class FDStencil2<OP,FIELDT,OP::DirT>;
#define DECLARE_FIELD_VARIANTS(OP)      \
    DECLARE_CLASS(OP,SVolField)         \
    DECLARE_CLASS(OP,XVolField)         \
    DECLARE_CLASS(OP,YVolField)         \
    DECLARE_CLASS(OP,ZVolField)

DECLARE_FIELD_VARIANTS( InterpolantX )
DECLARE_FIELD_VARIANTS( InterpolantY )
DECLARE_FIELD_VARIANTS( InterpolantZ )

DECLARE_FIELD_VARIANTS( GradientX )
DECLARE_FIELD_VARIANTS( GradientY )
DECLARE_FIELD_VARIANTS( GradientZ )

} // namespace structured
} // namespace SpatialOps
