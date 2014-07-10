/**
 *  \file   Example_field_type_inference.cpp
 *  \date   Jul 10, 2014
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

#include <spatialops/structured/FVStaggered.h> // field type definitions, etc.
#include <spatialops/SpatialOpsTools.h>        // for IsSameType utility

#include <cassert>

using namespace SpatialOps;

int main()
{
  typedef FaceTypes<SVolField>::XFace XFaceT;      // same as SSurfXField
  typedef FaceTypes<SVolField>::YFace YFaceT;      // same as SSurfYField
  typedef FaceTypes<SVolField>::ZFace ZFaceT;      // same as SSurfZField

  assert(( is_same_type< XFaceT, SSurfXField >() ));
  assert(( is_same_type< YFaceT, SSurfYField >() ));
  assert(( is_same_type< ZFaceT, SSurfZField >() ));


  typedef FaceTypes<XVolField>::XFace XMeshXFaceT; // same as XSurfXField
  typedef FaceTypes<XVolField>::YFace XMeshYFaceT; // same as XSurfYField
  typedef FaceTypes<XVolField>::ZFace XMeshZFaceT; // same as XSurfZField

  assert(( is_same_type< XMeshXFaceT, XSurfXField >() ));
  assert(( is_same_type< XMeshYFaceT, XSurfYField >() ));
  assert(( is_same_type< XMeshZFaceT, XSurfZField >() ));


  typedef FaceTypes<YVolField>::XFace YMeshXFaceT; // same as YSurfXField
  typedef FaceTypes<YVolField>::YFace YMeshYFaceT; // same as YSurfYField
  typedef FaceTypes<YVolField>::ZFace YMeshZFaceT; // same as YSurfZField

  assert(( is_same_type< YMeshXFaceT, YSurfXField >() ));
  assert(( is_same_type< YMeshYFaceT, YSurfYField >() ));
  assert(( is_same_type< YMeshZFaceT, YSurfZField >() ));


  typedef FaceTypes<ZVolField>::XFace ZMeshXFaceT; // same as ZSurfXField
  typedef FaceTypes<ZVolField>::YFace ZMeshYFaceT; // same as ZSurfYField
  typedef FaceTypes<ZVolField>::ZFace ZMeshZFaceT; // same as ZSurfZField

  assert(( is_same_type< ZMeshXFaceT, ZSurfXField >() ));
  assert(( is_same_type< ZMeshYFaceT, ZSurfYField >() ));
  assert(( is_same_type< ZMeshZFaceT, ZSurfZField >() ));

  return 0;
}
