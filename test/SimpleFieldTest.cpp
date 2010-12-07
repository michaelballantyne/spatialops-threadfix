#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/FieldOperations.h>
#include <spatialops/FieldOperationDefinitions.h>

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
  
  DEFINE_$_ANONYMOUS_ARGUMENTS(Field);
  
  // example of what we would like to do:
  // c = b + a + sin(b);
  
  c <<= (3 + a) - 2;
  c <<= 3 + (a + a);
  c <<= (3 - (3 * $0))(a);
  c <<= (a + 3) - a;
  c <<= (a * (a + 1)) / (a - 1);
  c <<= (a - (a * $0))(a);
  c <<= ((a + a) * ((a + a) - $0))(a);
  c <<= (($0 + 3) - 3)(a);
  c <<= (($0 + a) - a)(3.0);
  c <<= (($0 + (a + a)) - (a + a))(3.0);
  c <<= (($0 * $0) - $0)(a + a);
  c <<= ($0 - ($0 * $0))(a + a);
  c <<= (($1 * $2) - $0)(3 * a, a, 3.0);
  c <<= (($1 * $2) - $0)(3 * a, a)(3.0);
  c <<= (($1 * $2) - $0)(3 * a)(a, 3.0);
  c <<= (($1 * $2) - $0)(3 * a)(a)(3.0);
  
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
