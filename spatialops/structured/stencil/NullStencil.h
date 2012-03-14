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
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_NullStencil_h
