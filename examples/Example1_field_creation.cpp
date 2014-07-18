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

#include <spatialops/structured/FVStaggered.h>  // everything required to build fields on structured meshes

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
  // SVolField = Scalar Volume Field (non-staggered, cell-centered field)
  typedef SVolField FieldT;

  //----------------------------------------------------------------------------
  // Use default values to create objects required to build a field:

  // Ghost cells are needed by some applications and operations
  // In general, set ghost cells to zero, unless they are needed:
  const GhostData nghost(0);

  // Determine if we have physical boundaries present on each right/positive (+) face.
  const bool bcx=true, bcy=true, bcz=true;
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( bcx, bcy, bcz );

  // Define the size of the field (nx,ny,nz)
  const IntVec fieldDim( 10, 9, 8 );

  // Construct a memory window (logical extents of a field) from dimensions,
  //  ghost cell data, and boundary cell information
  const MemoryWindow window( get_window_with_ghost( fieldDim, nghost, bcInfo ) );

  //----------------------------------------------------------------------------
  // Create a field from scratch
  //  parameters:
  //   window          : logical extents for the field
  //   bcInfo          : boundary cell information for the field
  //   nghost          : ghost cell data for the field
  //   NULL            : externally allocated memory (not needed here; hence, NULL)
  //   InternalStorage : internally manage memory (could also be ExternalStorage)
  //   LOCATION        : CPU or GPU memory
  FieldT f( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );

  //----------------------------------------------------------------------------
  // Create a field from a prototype using the "SpatialFieldStore." SpatFldPtr has
  // regular pointer semantics but is a reference-counted pointer.
  SpatFldPtr<FieldT> f2 = SpatialFieldStore::get<FieldT>(f); // field with same layout as "f"

  return 0;
}


