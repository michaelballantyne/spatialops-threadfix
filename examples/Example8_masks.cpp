/**
 *  \file   Example8_masks.cpp
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
#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldHelper.h>

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
  //----------------------------------------------------------------------------
  // Define the domain size and number of points
  std::vector<double> domainLength(3,7.0);  // a cube of length 7.0
  // a 5 x 5 x 1 problem with a ghost cell on each edge
  const IntVec fieldDim( 7, 7, 1 );

  //----------------------------------------------------------------------------
  // Create fields
  const GhostData nghost(0);
  const BoundaryCellInfo sVolBCInfo = BoundaryCellInfo::build<SVolField>( true, true, true );
  const MemoryWindow sVolWindow( get_window_with_ghost( fieldDim, nghost, sVolBCInfo) );

  SVolField x( sVolWindow, sVolBCInfo, nghost, NULL, InternalStorage, LOCATION );
  SVolField y( sVolWindow, sVolBCInfo, nghost, NULL, InternalStorage, LOCATION );
  SVolField f( sVolWindow, sVolBCInfo, nghost, NULL, InternalStorage, LOCATION );

  //----------------------------------------------------------------------------
  // Build a coordinates.
  const Grid grid( fieldDim, domainLength );
  grid.set_coord<XDIR>(x);
  grid.set_coord<YDIR>(y);

  //----------------------------------------------------------------------------
  // Assign f:
  f <<= x + y;

  //----------------------------------------------------------------------------
  // Build a mask.
  // A mask is a set of indicies where something different should happen, such
  // as a boundary condition: Edge of a flame, fuel injection or exhaust pipe.
  // Indexing for masks has the origin at the first interior point. (Ghosts are
  // negative indices.)
  std::vector<IntVec> maskSet;
  for(int i = 0; i < 7; i++)
    maskSet.push_back(IntVec(0, i, 0));

  // Creating a mask needs a prototype field and a std::vector of indices
  SpatialMask<SVolField> mask(f, maskSet);

  //----------------------------------------------------------------------------
  // Use mask:

  std::cout << "f before applying mask:" << std::endl;
  print_field(f, std::cout);

  f <<= cond(mask, 10.0 - f)
            (f);

  std::cout << "f after applying mask:" << std::endl;
  print_field(f, std::cout);

  return 0;
}
