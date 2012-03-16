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

namespace SpatialOps{
  namespace structured{

    template< typename OpT, typename FieldT, typename DirT >
    class FDStencil2{
      const double coefLo_, coefHi_;

      typedef typename UnitTriplet<DirT>::type  DirVec;

    public:
      typedef OpT       Type;
      typedef FieldT    SrcFieldType;
      typedef FieldT    DestFieldType;

      FDStencil2( const double coefLo, const double coefHi );
      ~FDStencil2();

      void apply_to_field( const SrcFieldType& src, DestFieldType& dest ) const;

      inline double get_minus_coef() const{ return coefLo_; } ///< get the (-) coefficient
      inline double  get_plus_coef() const{ return coefHi_; } ///< get the (+) coefficient
    };

  }
}


#endif /* FDSTENCIL2_H_ */
