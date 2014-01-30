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

  const IntVec dim(10, 10, 1);
  const std::vector<bool> bcFlag(3,false);
  const GhostData  fghost(1);
  const GhostData dfghost(1);
  const BoundaryCellInfo  fbc = BoundaryCellInfo::build<PhiFieldT>( bcFlag[0], bcFlag[1], bcFlag[2] );
  const BoundaryCellInfo dfbc = BoundaryCellInfo::build<GammaFieldT>( bcFlag[0], bcFlag[1], bcFlag[2] );

  SpatFldPtr<PhiFieldT > f  = SpatialFieldStore::get_from_window<PhiFieldT>( get_window_with_ghost(dim, fghost, fbc),  fbc,  fghost );
  SpatFldPtr<GammaFieldT> df = SpatialFieldStore::get_from_window<GammaFieldT>( get_window_with_ghost(dim,dfghost,dfbc), dfbc, dfghost );

  PhiFieldT::iterator ifld=f->begin();
  for( int i = 0; i < 12; i++)
    for(int j=0; j < 12; j++) {
      *ifld = i;
      ifld++;
    }
  //*f <<= 3.0;

  GammaFieldT::iterator idfld=df->begin();
  for( int i = 0; i < 12; i++)
    for(int j=0; j < 12; j++) {
      *idfld = j;
      idfld++;
    }
  //*df <<= 1.0;

  //make the BC:
  NeboBoundaryConditionBuilder<IndexTriplet<-1,0,0>, IndexTriplet<0,0,0>, PhiFieldT, GammaFieldT> bc(1.0, 1.0);

  //make the mask:
  std::vector<IntVec> maskSet;
  maskSet.push_back(IntVec(0,5,0));
  //maskSet.push_back(IntVec(4,0,0));
  SpatialMask<SSurfXField> mask(*df, maskSet);

  typedef NeboMaskShift<Initial, IndexTriplet<1,0,0>, NeboMask<Initial, GammaFieldT>, PhiFieldT> ShiftedMask;
  ShiftedMask shiftedMask(NeboMask<Initial, GammaFieldT>(mask));

  cout << "Gamma" << endl;
  print_field(*df);

  cout << "Phi" << endl;
  print_field(*f);

  // evaluate the BC and set it in the field.
  bc(mask, *f, *df, true);

  // verify that the BC was set properly.
  cout << "Phi" << endl;
  print_field(*f);

  return 0;
}
