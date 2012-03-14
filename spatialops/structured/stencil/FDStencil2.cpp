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
  const MemoryWindow& w = src.window_with_ghost();

  const IntVec shift = DirVec::int_vec() + DirVec::int_vec();

  const MemoryWindow ws1( w.glob_dim(),
                          w.offset(),
                          w.extent() - shift,
                          w.has_bc(0), w.has_bc(1), w.has_bc(2) );
  const MemoryWindow ws2( w.glob_dim(),
                          w.offset() + shift,
                          w.extent() - shift,
                          w.has_bc(0), w.has_bc(1), w.has_bc(2) );
  const MemoryWindow wd(  w.glob_dim(),
                          w.offset() + DirVec::int_vec(),
                          w.extent() - shift,
                          w.has_bc(0), w.has_bc(1), w.has_bc(2) );

  FieldT  d( wd, dest.field_values(), ExternalStorage );
  FieldT s1( ws1, src.field_values(), ExternalStorage );
  FieldT s2( ws2, src.field_values(), ExternalStorage );

  typename FieldT::const_iterator is1=s1.begin(), is2=s2.begin();
  typename FieldT::iterator id=d.begin();
  const typename FieldT::iterator ide=d.end();
  for( ; id!=ide; ++id, ++is1, ++is2 ){
    *id = *is1*coefLo_ + *is2*coefHi_;
  }
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
