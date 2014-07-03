/**
 *  \file   Example1_field_creation.cpp
 *  \date   Jul 2, 2014
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
 */

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

using namespace SpatialOps;
using namespace structured;

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

  // Determine if we have physical boundaries present on each (+) face.
  const bool bcx=true, bcy=true, bcz=true;

  typedef SpatialOps::structured::SVolField FieldT;

  const GhostData nghost(1);
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( bcx, bcy, bcz );
  const MemoryWindow window( get_window_with_ghost( fieldDim, nghost, bcInfo ) );


  // Create a field with the specified size and have memory managed internally
  // (hence the NULL argument and the InternalStorage flag)
  // Also create this field in the appropriate location (GPU/CPU) depending on
  // how this was configured.
  FieldT f( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );

  return 0;
}


