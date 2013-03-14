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

#include "NullStencil.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{

  template< typename OpT, typename SrcT, typename DestT >
  NullStencil<OpT,SrcT,DestT>::NullStencil()
  {}

  template< typename OpT, typename SrcT, typename DestT >
  void
  NullStencil<OpT,SrcT,DestT>::
  apply_to_field( const SrcT& src, DestT& dest ) const
  {
#   ifndef NDEBUG
    assert( src.window_with_ghost() == dest.window_with_ghost() );
#   endif

    dest <<= operator()(src);
  }

  //==================================================================
  // Explicit template instantiation
  //
  template struct NullStencil< Interpolant, SVolField, SVolField >;  // no-op interpolant
  template struct NullStencil< Interpolant, XVolField, XVolField >;  // no-op interpolant
  template struct NullStencil< Interpolant, YVolField, YVolField >;  // no-op interpolant
  template struct NullStencil< Interpolant, ZVolField, ZVolField >;  // no-op interpolant

  template struct NullStencil< Interpolant, XVolField, SSurfXField >;  // x-advecting velocity on scalar x-surface
  template struct NullStencil< Interpolant, YVolField, SSurfYField >;  // y-advecting velocity on scalar y-surface
  template struct NullStencil< Interpolant, ZVolField, SSurfZField >;  // z-advecting velocity on scalar z-surface

  template struct NullStencil< Interpolant, SVolField, XSurfXField >;  // viscosity on x-x-face
  template struct NullStencil< Interpolant, SVolField, YSurfYField >;  // viscosity on y-y-face
  template struct NullStencil< Interpolant, SVolField, ZSurfZField >;  // viscosity on z-z-face
  //
  //==================================================================


} // namespace structured
} // namespace SpatialOps
