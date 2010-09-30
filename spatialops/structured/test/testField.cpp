#include <spatialops/structured/FVStaggeredTypes.h>
#include <test/TestHelper.h>

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;

int main()
{
  TestHelper overall(true);

  // test iterators
  {
    TestHelper status(false);
    const IntVec npts(3,3,3);
    const MemoryWindow window(npts);
    SVolField f1( window, NULL );
    SVolField f2( window, NULL );
    f1 = 2.0;
    f2 = 1.0;

    SVolField::iterator if2=f2.begin();
    const SVolField::iterator if2e=f2.end();
    SVolField::const_iterator if1=f1.begin();
    for( ; if2!=if2e; ++if1, ++if2 ){
      *if2 += *if1;
      status( *if2 == 3.0 );
    }
    overall( status.ok(), "iterator test" );
  }

  // test basic layout and operators
  {
    TestHelper status(false);

    const int npts[3] = {10,11,12};
    const MemoryWindow window(npts);
    SVolField svol1( window, NULL, InternalStorage );
    SVolField svol2( window, NULL, InternalStorage );

    for( int k=0; k<npts[2]; ++k ){
      for( int j=0; j<npts[1]; ++j ){
        for( int i=0; i<npts[0]; ++i ){
          svol1(i,j,k) = i + j + k;
        }
      }
    }

    svol2 = 2.0;
    svol1 += svol2;
    svol1 *= svol2;
    svol1 /= svol2;
    svol1 -= svol2;

    for( int k=0; k<npts[2]; ++k ){
      for( int j=0; j<npts[1]; ++j ){
        for( int i=0; i<npts[0]; ++i ){
          const double ans = (i + j + k);
          status( ans==svol1(i,j,k) );
          if( svol1(i,j,k) != ans ){
            cout << "(i,j,k)=("<<i<<","<<j<<","<<k<<"),  "
                 << svol1(i,j,k) << ", " << ans << endl;
          }
        }
      }
    }
    overall( status.ok(), "field operations" );
  }
  if( overall.isfailed() ) return -1;
  return 0;
}
