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
#include <spatialops/structured/FieldHelper.h> // convenient way to view small fields

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
  const IntVec fieldDim( 5, 5, 1 );

  //----------------------------------------------------------------------------
  // Create fields:
  typedef SpatialOps::SVolField FieldT;

  // Use default values to create objects that are required to construct a field.
  const GhostData nghost(0);
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( true, true, true );
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
  // Print fields to standard output:
  std::cout << "x:" << std::endl;
  print_field(x, std::cout);
  std::cout << "y:" << std::endl;
  print_field(y, std::cout);
  std::cout << "z:" << std::endl;
  print_field(z, std::cout);

  //----------------------------------------------------------------------------
  // Perform operations on fields

  // Assign field values.  Note that this is a vectorized statement that will
  // work on GPU, serial CPU and multicore CPU.
  f <<= sin(x) + cos(y) + tanh(z);

#ifdef ENABLE_CUDA
  //If f uses GPU memory, to print f, f needs to be copied to CPU memory.
  f.add_device_sync(CPU_INDEX);
#endif
  std::cout << "f:" << std::endl;
  print_field(f, std::cout);

  //----------------------------------------------------------------------------
  // Conditional (if/then/else...) evaluation.
  // cond in Nebo creates conditional expressions. Each cond clause, except the last,
  // must have two arguments.  The first argument is the condition (if this), and the
  // second is the result (then that).  The final clause takes only one arguement
  // (else other), which is returned only if all previous conditions fail.
  //
  // The conditions are evaluated first to last, and evaluation stops when a
  // condition returns true. Thus, the order of conditions can and will effect
  // the results.
  g <<= cond( f > 2.0, x+z )          // if     ( f[i] > 2.0 ) g[i] = x[i]+z[i];
            ( f > 1.5, y   )          // else if( f[i] > 1.5 ) g[i] = y[i];
            ( f > 1.0, -f  )          // else if( f[i] > 1.0 ) g[i] = -f[i];
            ( f );                    // else                  g[i] = f[i];

#ifdef ENABLE_CUDA
  //If g uses GPU memory, to print g, g needs to be copied to CPU memory.
  g.add_device_sync(CPU_INDEX);
#endif
  std::cout << "g:" << std::endl;
  print_field(g, std::cout);

  return 0;
}
