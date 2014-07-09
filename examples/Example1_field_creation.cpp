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
  // Define the size of the field (nx,ny,nz)
  const IntVec fieldDim( 10, 9, 8 );

  // Determine if we have physical boundaries present on each (+) face.
  const bool bcx=true, bcy=true, bcz=true;


  //----------------------------------------------------------------------------
  // Create a field with the specified size and have memory managed internally
  // (hence the NULL argument and the InternalStorage flag)
  // Also create this field in the appropriate location (GPU/CPU) depending on
  // how this was configured.
  typedef SVolField FieldT;  // SVolField = Scalar Volume Field (non-staggered, cell-centered field)
  const GhostData nghost(1);
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( bcx, bcy, bcz );
  const MemoryWindow window( get_window_with_ghost( fieldDim, nghost, bcInfo ) );

  FieldT f( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );


  //----------------------------------------------------------------------------
  // Create some fields from the "SpatialFieldStore" using the previously
  // created field as a prototype.  SpatFldPtr has regular pointer semantics
  // but is a reference-counted pointer.

  SpatFldPtr<FieldT> f2 = SpatialFieldStore::get<FieldT>(f); // field with same layout as "f"


  //----------------------------------------------------------------------------
  // f and f2 are volume fields on the scalar (non-staggered) mesh.  Now let's
  // create a few surface fields.  We will use type inference to get the face
  // field associated with the cell field type that we are using

  typedef FaceTypes<FieldT>::XFace XFaceT;  // X-face field on the scalar mesh
  typedef FaceTypes<FieldT>::YFace YFaceT;  // Y-face field on the scalar mesh
  typedef FaceTypes<FieldT>::ZFace ZFaceT;  // Z-face field on the scalar mesh

  SpatFldPtr<XFaceT> fx = SpatialFieldStore::get<XFaceT>(f); // field on x-face of the same grid as "f"
  SpatFldPtr<YFaceT> fy = SpatialFieldStore::get<YFaceT>(f); // field on y-face of the same grid as "f"
  SpatFldPtr<ZFaceT> fz = SpatialFieldStore::get<ZFaceT>(f); // field on z-face of the same grid as "f"

  return 0;
}


