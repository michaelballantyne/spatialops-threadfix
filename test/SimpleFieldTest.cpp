#include <spatialops/structured/FVStaggeredTypes.h>

int main()
{
  const int nx=10, ny=12, nz=14;
  const int npts = nx*ny*nz;

  std::set<int> ghostIndices;

  typedef SpatialOps::structured::SVolField Field;

  Field a( npts, ghostIndices, NULL );
  Field b( npts, ghostIndices, NULL );
  Field c( npts, ghostIndices, NULL );

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
}
