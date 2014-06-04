#include <spatialops/SpatialOpsTools.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <spatialops/OperatorDatabase.h>

#include <spatialops/Nebo.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/SpatialMask.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/structured/stencil/FVStaggeredBCOp.h>
#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/FVStaggeredBCTools.h>
#include <spatialops/structured/stencil/StencilBuilder.h>

#include <spatialops/structured/Grid.h>
#include <test/TestHelper.h>
#include <test/FieldHelper.h>

using namespace SpatialOps;
using namespace structured;

int main()
{
  TestHelper status(true);

  using namespace SpatialOps;
  using namespace structured;

  typedef SVolField PhiFieldT;
  typedef SSurfXField GammaFieldT;

  const IntVec dim(3, 3, 1);
  const double length = 1.5;
  const double dx = length/dim[0];
  const std::vector<bool> bcFlag(3,false);
  const GhostData  fghost(1);
  const GhostData dfghost(1);
  const BoundaryCellInfo  fbc = BoundaryCellInfo::build<PhiFieldT>( bcFlag[0], bcFlag[1], bcFlag[2] );
  const BoundaryCellInfo dfbc = BoundaryCellInfo::build<GammaFieldT>( bcFlag[0], bcFlag[1], bcFlag[2] );

  SpatFldPtr<PhiFieldT  > test  = SpatialFieldStore::get_from_window<PhiFieldT>( get_window_with_ghost(dim, fghost, fbc),  fbc,  fghost );
  SpatFldPtr<PhiFieldT  > ref   = SpatialFieldStore::get_from_window<PhiFieldT>( get_window_with_ghost(dim, fghost, fbc),  fbc,  fghost );
  SpatFldPtr<GammaFieldT> gamma = SpatialFieldStore::get_from_window<GammaFieldT>( get_window_with_ghost(dim,dfghost,dfbc), dfbc, dfghost );

  /* Invert interpolant:
   *  Cell (volume) to face:
   *  -------------
   *  |     |     |
   *  |  A  B  C  |
   *  |     |     |
   *  --------------
   *
   *  dx = Lx/nx
   *
   *  Gradient:
   *  B = A * -1/dx + C * 1/dx
   *
   *  Negative side inversion:
   *  A = (B - C * 1/dx) / (-1/dx)
   *
   *  Positive side inversion:
   *  C = (B - A * -1/dx) / (1/dx)
   *
   *  Indices:
   *  A = -1, 0, 0
   *  B =  0, 0, 0
   *  C =  0, 0, 0
   */

  GammaFieldT::iterator ig = gamma->begin();
  PhiFieldT::iterator ir = ref->begin();
  PhiFieldT::iterator it = test->begin();
  for(int j=-1; j < 4; j++)
    for( int i = -1; i < 4; i++) {
      *ig = 25 + i + j * 5;
      *it = i + j * 5;
      if( (i == -1 && j == 1) ||
          (i ==  0 && j == 2) ) {
        //Negative side inversion: *ir == A
        int B = 25 + (i + 1) + j * 5;
        int C = (i + 1) + j * 5;
        *ir = (B - C * (1/dx)) / (-1/dx);
      }
      else if( (i ==  2 && j == 1) ||
               (i ==  3 && j == 2) ) {
        //Positive side inversion: *ir == C
        int B = 25 + (i) + j * 5;
        int A = (i - 1) + j * 5;
        *ir = (B - A * (-1/dx)) / (1/dx);
      }
      else
        *ir = i + j * 5;
      ig++;
      it++;
      ir++;
    }

  //make the BC:
  OperatorDatabase opdb;
  build_stencils( dim[0], dim[1], dim[2], length, length, length, opdb );
  typedef structured::BasicOpTypes<PhiFieldT>::GradX OpT;
  const OpT* const op = opdb.retrieve_operator<OpT>();
  NeboBoundaryConditionBuilder<OpT> BC(*op);

  //make the minus mask:
  std::vector<IntVec> minusSet;
  minusSet.push_back(IntVec(0,1,0));
  minusSet.push_back(IntVec(1,2,0));
  SpatialMask<GammaFieldT> minus(*gamma, minusSet);

  // evaluate the minus BC and set it in the field.
  BC(minus, *test, *gamma, true);

  //make the plus mask:
  std::vector<IntVec> plusSet;
  plusSet.push_back(IntVec(2,1,0));
  plusSet.push_back(IntVec(3,2,0));
  //build mask without gamma:
  SpatialMask<GammaFieldT> plus = SpatialMask<GammaFieldT>::build(*test, plusSet);

  // evaluate the plus BC and set it in the field.
  BC(plus, *test, *gamma, false);

  //display differences and values:
  //display_fields_compare(*test, *ref, true, true);

  // verify that the BC was set properly.
  if( field_equal(*test, *ref) ) {
    std::cout << "ALL TESTS PASSED :)" << std::endl;
    return 0;
  } else {
    std::cout << "******************************" << std::endl
              << " At least one test FAILED! :(" << std::endl
              << "******************************" << std::endl;
    return -1;
  }

}
