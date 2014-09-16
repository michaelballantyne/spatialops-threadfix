#include <spatialops/structured/FVStaggeredFieldTypes.h>

#include <spatialops/Nebo.h>
#include <test/TestHelper.h>
#include <spatialops/structured/FieldHelper.h>

#include <spatialops/structured/FieldComparisons.h>

#include <iostream>

using namespace SpatialOps;
using std::cout;
using std::endl;

typedef SVolField Field;

class FTest {
    public:
        double operator()(double i1, double i2, double i3, double i4) const {
            return i1 + i2 + i3 + i4;
        }
};

int main( int iarg, char* carg[] )
{
    int nx, ny, nz;
    bool bcplus[] = { false, false, false };

    nx = 11;
    ny = 11;
    nz = 11;

    const int nghost = 1;
    const GhostData ghost(nghost);
    const BoundaryCellInfo bcinfo = BoundaryCellInfo::build<Field>(bcplus[0],bcplus[1],bcplus[2]);
    const MemoryWindow window( get_window_with_ghost(IntVec(nx,ny,nz),ghost,bcinfo) );

    Field input1( window, bcinfo, ghost, NULL );
    Field input2( window, bcinfo, ghost, NULL );
    Field input3( window, bcinfo, ghost, NULL );
    Field   test( window, bcinfo, ghost, NULL );
    Field    ref( window, bcinfo, ghost, NULL );
    const int total = nx * ny * nz;

    initialize_field(input1, 0.0);
    initialize_field(input2, total);
    initialize_field(input3, total);

    TestHelper status(true);

    ref <<= (input1 + input2) + input3 + 5.0 + 6.0;

    test <<= apply_pointwise<FTest>(input1 + input2, input3, 5.0, 6.0);

    status(field_equal(ref, test, 0.0), "simple pointwise");

    if( status.ok() ) {
      cout << "ALL TESTS PASSED :)" << endl;
      return 0;
    }
    else {
        cout << "******************************" << endl
             << " At least one test FAILED! :(" << endl
             << "******************************" << endl;
        return -1;
    };
}
