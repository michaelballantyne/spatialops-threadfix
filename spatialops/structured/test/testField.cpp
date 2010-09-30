#include <spatialops/structured/FVStaggeredTypes.h>

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;

int main()
{
  bool isokay = true;

  // test iterators
  {
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
    }
  }

  // test basic layout and operators
  {
    const int npts[3] = {10,11,12};
    const MemoryWindow window(npts);
    std::cout << window.flat_index(window.extent()) << ", " << window.npts() << std::endl;
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
          if( svol1(i,j,k) != ans ){
            isokay = false;
            cout << "(i,j,k)=("<<i<<","<<j<<","<<k<<"),  "
                 << svol1(i,j,k) << ", " << ans << endl;
          }
        }
      }
    }

  }
  if( isokay ) return 0;
  return -1;
}
