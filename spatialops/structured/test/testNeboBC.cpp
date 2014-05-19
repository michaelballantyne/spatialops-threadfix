#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

using std::cout;
using std::endl;

#include <spatialops/OperatorDatabase.h>

#include <spatialops/Nebo.h>
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
  const std::vector<bool> bcFlag(3,false);
  const GhostData  fghost(1);
  const GhostData dfghost(1);
  const BoundaryCellInfo  fbc = BoundaryCellInfo::build<PhiFieldT>( bcFlag[0], bcFlag[1], bcFlag[2] );
  const BoundaryCellInfo dfbc = BoundaryCellInfo::build<GammaFieldT>( bcFlag[0], bcFlag[1], bcFlag[2] );

  SpatFldPtr<PhiFieldT  > test  = SpatialFieldStore::get_from_window<PhiFieldT>( get_window_with_ghost(dim, fghost, fbc),  fbc,  fghost );
  SpatFldPtr<PhiFieldT  > ref   = SpatialFieldStore::get_from_window<PhiFieldT>( get_window_with_ghost(dim, fghost, fbc),  fbc,  fghost );
  SpatFldPtr<GammaFieldT> gamma = SpatialFieldStore::get_from_window<GammaFieldT>( get_window_with_ghost(dim,dfghost,dfbc), dfbc, dfghost );

  GammaFieldT::iterator ig = gamma->begin();
  PhiFieldT::iterator ir = ref->begin();
  PhiFieldT::iterator it = test->begin();
  for(int j=-1; j < 4; j++)
    for( int i = -1; i < 4; i++) {
      *ig = 25 + i + j * 5;
      *it = i + j * 5;
      if( (i == -1 && j == 1) ||
          (i ==  0 && j == 2) ||
          (i ==  2 && j == 1) ||
          (i ==  3 && j == 2)
          )
        *ir = 25;
      else
        *ir = i + j * 5;
      ig++;
      it++;
      ir++;
    }

  //make the minus BC:
  NeboBoundaryConditionBuilder<IndexTriplet<-1,0,0>, IndexTriplet<0,0,0>, PhiFieldT, GammaFieldT> minusBC(1.0, 1.0);

  //make the minus mask:
  std::vector<IntVec> minusSet;
  minusSet.push_back(IntVec(0,1,0));
  minusSet.push_back(IntVec(1,2,0));
  SpatialMask<SSurfXField> minus(*gamma, minusSet);

  // evaluate the minus BC and set it in the field.
  minusBC(minus, *test, *gamma, true);

  //make the plus BC:
  NeboBoundaryConditionBuilder<IndexTriplet<0,0,0>, IndexTriplet<-1,0,0>, PhiFieldT, GammaFieldT> plusBC(1.0, 1.0);

  //make the plus mask:
  std::vector<IntVec> plusSet;
  plusSet.push_back(IntVec(2,1,0));
  plusSet.push_back(IntVec(3,2,0));
  SpatialMask<SSurfXField> plus(*gamma, plusSet);

  // evaluate the plus BC and set it in the field.
  plusBC(plus, *test, *gamma, false);

  cout << "test" << endl;
  print_field(*test);

  // verify that the BC was set properly.
  if( display_fields_compare(*test, *ref, true, true) ) {
    cout << "ALL TESTS PASSED :)" << endl;
    return 0;
  } else {
    cout << "******************************" << endl
         << " At least one test FAILED! :(" << endl
         << "******************************" << endl;
    return -1;
  }

}
