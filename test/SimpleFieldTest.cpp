#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>

using namespace SpatialOps;
using namespace structured;

int main()
{
  const int nx=10, ny=12, nz=14;
  const int npts = nx*ny*nz;

  typedef SVolField Field;

  const MemoryWindow window( get_window_with_ghost<Field>(IntVec(nx,ny,nz),true,true,true) );

  Field a( window, NULL );
  Field b( window, NULL );
  Field c( window, NULL );

  a = 1.0;
  b = 1.0;

  // example of what we would like to do:
  //  c = b + a + sin(b);

  // what we currently must do:
  Field::const_iterator ia = a.begin();
  for( Field::iterator ic=c.begin(); ic!=c.end(); ++ic, ++ia ){
    *ic = sin(*ia);
  }
  c += a;  // operator overload
  c += b;  // operator overload

  c *= a;
  c /= a;

  return 0;
}
