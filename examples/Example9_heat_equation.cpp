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

// Initialize alpha, thermal conductivity:
template<typename FieldT>
void initialize_alpha(Grid const & grid,
                      FieldT & alpha) {
  const GhostData nghost(1);
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( true, true, true );
  const MemoryWindow window( get_window_with_ghost( grid.extent(), nghost, bcInfo) );
  FieldT     x( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );
  FieldT     y( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );

  grid.set_coord<XDIR>(x);
  grid.set_coord<YDIR>(y);
  //alpha <<= 0.1 + 0.4 * (x + y + 1.0);
  //alpha <<= 0.1 + 0.4 * (x + y);
  //alpha <<= (x + y + .2) /2;
  alpha <<= 1.0;

  std::cout << "Initial alpha:" << std::endl;
  print_field(alpha, std::cout, true);
}

void initialize_mask_points(Grid const & grid,
                            std::vector<IntVec> & leftSet,
                            std::vector<IntVec> & rightSet,
                            std::vector<IntVec> & topSet,
                            std::vector<IntVec> & bottomSet) {
  for(int i = -1; i <= grid.extent()[1]; i++)
    leftSet.push_back(IntVec(-1, i, 0));

  for(int i = -1; i <= grid.extent()[1]; i++)
    rightSet.push_back(IntVec(grid.extent()[0], i, 0));

  for(int i = -1; i <= grid.extent()[0]; i++)
    topSet.push_back(IntVec(i, -1, 0));

  for(int i = -1; i <= grid.extent()[0]; i++)
    bottomSet.push_back(IntVec(i, grid.extent()[1], 0));
}

template<typename FieldT>
double find_deltaT(FieldT const & alpha,
                   Grid const & grid) {
  const double deltaX = grid.spacing<XDIR>();
  const double deltaY = grid.spacing<YDIR>();
  const double sqrdDeltaX = deltaX * deltaX;
  const double sqrdDeltaY = deltaY * deltaY;
  const double sqrdDeltaXYmult = sqrdDeltaX * sqrdDeltaY;
  const double sqrdDeltaXYplus = sqrdDeltaX + sqrdDeltaY;

  return 0.25 * sqrdDeltaXYmult / ( sqrdDeltaXYplus * nebo_min(alpha) );
}

int main()
{
  typedef SVolField FieldT;

  //----------------------------------------------------------------------------
  // Define the domain size and number of points
  std::vector<double> domainLength(3,1.0);  // a cube of unit length
  // a 6 x 6 x 1 problem with a ghost cell on each edge
  const IntVec fieldDim( 6, 6, 1 );

  //----------------------------------------------------------------------------
  // Create fields
  const GhostData nghost(1);
  const BoundaryCellInfo bcInfo = BoundaryCellInfo::build<FieldT>( true, true, true );
  const MemoryWindow window( get_window_with_ghost( fieldDim, nghost, bcInfo) );

  FieldT   phi( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );
  FieldT     Q( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );
  FieldT alpha( window, bcInfo, nghost, NULL, InternalStorage, LOCATION );

  const Grid grid( fieldDim, domainLength );

  //----------------------------------------------------------------------------
  // Initialize alpha, thermal conductivity:
  initialize_alpha(grid, alpha);

  //----------------------------------------------------------------------------
  // Build and initialize masks:
  std::vector<IntVec> leftSet;
  std::vector<IntVec> rightSet;
  std::vector<IntVec> topSet;
  std::vector<IntVec> bottomSet;

  initialize_mask_points(grid, leftSet, rightSet, topSet, bottomSet);

  SpatialMask<FieldT> left(phi, leftSet);
  SpatialMask<FieldT> right(phi, rightSet);
  SpatialMask<FieldT> top(phi, topSet);
  SpatialMask<FieldT> bottom(phi, bottomSet);

  //----------------------------------------------------------------------------
  // Build stencils:
  OperatorDatabase opDB;         // holds stencils that can be retrieved easily
  build_stencils( grid, opDB );  // builds stencils and places them in opDB

  typedef BasicOpTypes<SVolField>::DivX          DivX;  // x-divergence operator type
  typedef BasicOpTypes<SVolField>::GradX        GradX;  // x-gradient operator type
  typedef BasicOpTypes<SVolField>::InterpC2FX InterpX;  // x-interpolant operator type

  const DivX&       divX = *opDB.retrieve_operator<DivX   >(); // retrieve the DivX operator
  const GradX&     gradX = *opDB.retrieve_operator<GradX  >(); // retrieve the GradX operator
  const InterpX& interpX = *opDB.retrieve_operator<InterpX>(); // retrieve the InterpX operator

  typedef BasicOpTypes<SVolField>::DivY          DivY;  // y-divergence operator type
  typedef BasicOpTypes<SVolField>::GradY        GradY;  // y-gradient operator type
  typedef BasicOpTypes<SVolField>::InterpC2FY InterpY;  // y-interpolant operator type

  const DivY&       divY = *opDB.retrieve_operator<DivY   >(); // retrieve the DivY operator
  const GradY&     gradY = *opDB.retrieve_operator<GradY  >(); // retrieve the GradY operator
  const InterpY& interpY = *opDB.retrieve_operator<InterpY>(); // retrieve the InterpY operator

  //----------------------------------------------------------------------------
  // Determine a safe deltaT:
  const double deltaT = find_deltaT(alpha, grid);

  //----------------------------------------------------------------------------
  // Initialize phi:
  phi <<= cond(left, 10.0)
              (right, 0.0)
              (5.0);

  std::cout << "Initial phi:" << std::endl;
  print_field(phi, std::cout, true);

  //----------------------------------------------------------------------------
  // Take time steps:
  const int total = 40;
  const int print_every = 1;
  int count = 1;
  int i = 0;
  for(; i <= total; i++, count++) {
    //Calculate change in phi:
    Q <<= divX( interpX(alpha) * gradX(phi) ) +
          divY( interpY(alpha) * gradY(phi) );

    //Integrate into phi:
    phi <<= phi + deltaT * Q;

    //Reset boundaries:
    phi <<= cond(left, 10.0)
                (right, 0.0)
                (top || bottom, 5.0)
                (phi);

    //print current state:
    if(count >= print_every) {
      count = 0;
      std::cout << "phi after " << i + 1 << " time steps:" << std::endl;
      print_field(phi, std::cout, true);
    }
  }

  return 0;
}
