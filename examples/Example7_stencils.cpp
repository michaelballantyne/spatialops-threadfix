/**
 *  \file   Example7_stencils.cpp
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

#include <spatialops/structured/FVStaggered.h>
#include <spatialops/structured/Grid.h>

using namespace SpatialOps;


// If we are compiling with GPU CUDA support, create fields on the device.
// Otherwise, create them on the host.
#ifdef ENABLE_CUDA
# define LOCATION GPU_INDEX
#else
# define LOCATION CPU_INDEX
#endif


#define PI 3.141592653589793

int main()
{
  //----------------------------------------------------------------------------
  // Define the domain size and number of points
  std::vector<double> domainLength(3,PI);  // a cube of length pi on each side
  const IntVec fieldDim( 10, 1, 1 );       // a 1-D problem


  //----------------------------------------------------------------------------
  // Create fields
  const bool bcx=true, bcy=true, bcz=true;
  const GhostData nghost(1);
  const BoundaryCellInfo sVolBCInfo = BoundaryCellInfo::build<SVolField>( bcx, bcy, bcz );
  const BoundaryCellInfo ssxBCInfo = BoundaryCellInfo::build<SSurfXField>( bcx, bcy, bcz );
  const MemoryWindow sVolWindow( get_window_with_ghost( fieldDim, nghost, sVolBCInfo) );
  const MemoryWindow ssxWindow( get_window_with_ghost( fieldDim, nghost, ssxBCInfo ) );

  SVolField       x( sVolWindow, sVolBCInfo, nghost, NULL, InternalStorage, LOCATION );
  SVolField       f( sVolWindow, sVolBCInfo, nghost, NULL, InternalStorage, LOCATION );
  SSurfXField  dfdx(  ssxWindow,  ssxBCInfo, nghost, NULL, InternalStorage, LOCATION );
  SSurfXField fface(  ssxWindow,  ssxBCInfo, nghost, NULL, InternalStorage, LOCATION );

  //----------------------------------------------------------------------------
  // Build a grid. This is a convenient way to set coordinate values that will
  // be used below.
  const Grid grid( fieldDim, domainLength );
  grid.set_coord<XDIR>(x);


  //----------------------------------------------------------------------------
  // Build the operators (stencils).  Here we will use predefined stencils that
  // can be built easily.
  OperatorDatabase opDB;         // holds stencils that can be retrieved easily
  build_stencils( grid, opDB );  // builds stencils and places them in opDB

  typedef BasicOpTypes<SVolField>::GradX        GradX;  // x-derivative operator type
  typedef BasicOpTypes<SVolField>::InterpC2FX InterpX;  // x-interpolant operator type

  const GradX&     gradx = *opDB.retrieve_operator<GradX  >(); // retrieve the GradX operator
  const InterpX& interpx = *opDB.retrieve_operator<InterpX>(); // retrieve the InterpX operator


  //----------------------------------------------------------------------------
  // perform calculations
  f     <<=     sin( x ); // Initialize the function value
  dfdx  <<=   gradx( f ); // gradient of f at the x-faces
  fface <<= interpx( f ); // interpolated value of f at the x-faces

  return 0;
}
