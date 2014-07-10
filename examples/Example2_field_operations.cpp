/**
 *  \file   Example2_field_operations.cpp
 *  \date   Jul 6, 2014
 *  \author "James C. Sutherland"
 *
 *
 * The MIT License
 *
 * Copyright (c) 2014 The University of Utah
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
 *
 */

#include <spatialops/structured/FVStaggered.h>
#include <spatialops/structured/Grid.h> // convenient way to define coordinates
using namespace SpatialOps;

// If we are compiling with GPU CUDA support, create fields on the device.
// Otherwise, create them on the host.
#ifdef ENABLE_CUDA
# define LOCATION GPU_INDEX
#else
# define LOCATION CPU_INDEX
#endif

int main()
{
  // Define the size of the field (nx,ny,nz)
  const IntVec fieldDim( 10, 9, 8 );


  //----------------------------------------------------------------------------
  // Create fields
  typedef SpatialOps::SVolField FieldT;

  // set some objects that are requried to construct a field.
  // Don't worry about these too much for now. Just use these values as defaults.
  const bool bcx=true, bcy=true, bcz=true;
  const GhostData nghost(0);
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( bcx, bcy, bcz );
  const MemoryWindow window( get_window_with_ghost( fieldDim, nghost, bcInfo ) );

  FieldT x( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );
  FieldT y( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );
  FieldT z( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );

  FieldT f( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );
  FieldT g( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );

  //----------------------------------------------------------------------------
  // Build a grid. This is a convenient way to set coordinate values that will
  // be used below.
  std::vector<double> domainLength(3,1.0);  // a cube of unit length
  const Grid grid( fieldDim, domainLength );
  grid.set_coord<XDIR>(x);
  grid.set_coord<YDIR>(y);
  grid.set_coord<ZDIR>(z);

  //----------------------------------------------------------------------------
  // Perform operations on fields

  // Assign field values.  Note that this is a vectorized statement that will
  // work on GPU, serial CPU and multicore CPU.
  f <<= sin(x) + cos(y) + tanh(z);

  // Field reduction operations - currently only supported on CPU
  const double fmax = field_max( f );                   // maximum of a field
  const double max2 = field_max( sin(x)*cos(x) + 2.0 ); // maximum of an expression
  const double fnorm = field_norm( f );                 // L2 norm of f

  g <<= field_max(f) + field_min(f) + exp(-x-y-z);      // combine several field operations

  // conditional (if/then/else...)
  g <<= cond( f >  0, x+z )          // if     ( f[i] >  0 ) g[i] = x[i]+z[i];
            ( f < -2, y   )          // else if( f[i] < -2 ) g[i] = y[i];
            ( f < -1, 3.4 )          // else if( f[i] < -1 ) g[i] = 3.45
            ( f );                   // else                 g[i] = f[i];

  return 0;
}
