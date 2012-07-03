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

#include "Stencil2_def.h"


namespace SpatialOps{
namespace structured{

  DECLARE_BASIC_VARIANTS( SVolField );
  DECLARE_BASIC_VARIANTS( XVolField );
  DECLARE_BASIC_VARIANTS( YVolField );
  DECLARE_BASIC_VARIANTS( ZVolField );

  DECLARE_STENCIL( Interpolant, XVolField, YSurfXField )  // advecting velocity
  DECLARE_STENCIL( Gradient   , XVolField, YSurfXField )  // stress
  DECLARE_STENCIL( Interpolant, XVolField, ZSurfXField )  // advecting velocity
  DECLARE_STENCIL( Gradient   , XVolField, ZSurfXField )  // stress

  DECLARE_STENCIL( Interpolant, YVolField, XSurfYField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    YVolField, XSurfYField )  // stress
  DECLARE_STENCIL( Interpolant, YVolField, ZSurfYField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    YVolField, ZSurfYField )  // stress

  DECLARE_STENCIL( Interpolant, ZVolField, XSurfZField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    ZVolField, XSurfZField )  // stress
  DECLARE_STENCIL( Interpolant, ZVolField, YSurfZField )  // advecting velocity
  DECLARE_STENCIL( Gradient,    ZVolField, YSurfZField )  // stress

  DECLARE_STENCIL( Interpolant, SVolField, XVolField )  // density, dp/dx
  DECLARE_STENCIL( Interpolant, SVolField, YVolField )  // density, dp/dy
  DECLARE_STENCIL( Interpolant, SVolField, ZVolField )  // density, dp/dz

  DECLARE_STENCIL( Interpolant, XVolField, SVolField )  // pressure projection RHS
  DECLARE_STENCIL( Interpolant, YVolField, SVolField )  // pressure projection RHS
  DECLARE_STENCIL( Interpolant, ZVolField, SVolField )  // pressure projection RHS

  DECLARE_STENCIL( Gradient, XVolField, SVolField )  // dilatation
  DECLARE_STENCIL( Gradient, YVolField, SVolField )  // dilatation
  DECLARE_STENCIL( Gradient, ZVolField, SVolField )  // dilatation

  DECLARE_STENCIL( Gradient, SVolField, XVolField )  // pressure
  DECLARE_STENCIL( Gradient, SVolField, YVolField )  // pressure
  DECLARE_STENCIL( Gradient, SVolField, ZVolField )  // pressure

  DECLARE_STENCIL( Interpolant, SSurfXField, SVolField )  // ODT colocated mesh

  DECLARE_STENCIL( Interpolant, XSurfXField, XVolField )  // BC operator for tau_xx
  DECLARE_STENCIL( Interpolant, YSurfYField, YVolField )  // BC operator for tau_yy
  DECLARE_STENCIL( Interpolant, ZSurfZField, ZVolField )  // BC operator for tau_zz
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
