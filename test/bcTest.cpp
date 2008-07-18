#include <vector>
#include <cmath>

#include <FVStaggeredBCTools.h>
#include <FVStaggered.h>

#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>  // use if you need to bind function arguments or to bind class member functions...


//====================================================================

double tmp()
{
  return 1.234;
}

//====================================================================

double tmp2( const double t )
{
  return 2*t;
}

//====================================================================

double ttmp()
{
  return 1.1;
}

//====================================================================

void test1()
{
  using namespace SpatialOps;
  using namespace FVStaggered;

  const double ATOL = 1.0e-5;
  const int nx=3, ny=4, nz=5;

  std::vector<int> dim(3);  dim[0]=nx; dim[1]=ny; dim[2]=nz;
  const bool bcx=true,  bcy=true, bcz=true;

  typedef boost::function< double() > BCFun;
  BCFun fnotime = &tmp;

  // composite function:  tmp2( ttmp() )
  BCFun tfun = boost::lambda::bind( &ttmp );
  BCFun f = boost::lambda::bind( &tmp2, boost::lambda::bind(&ttmp) );

  // pick some points to set bcs on.  Note that they need not actually
  // reside on a boundary.
  vector<IndexTriplet> pts;
  pts.push_back( IndexTriplet( 0,    0,    0    ) );
  pts.push_back( IndexTriplet( 1,    0,    0    ) );
  pts.push_back( IndexTriplet( 2,    1,    2    ) );
  pts.push_back( IndexTriplet( 0,    0,    1    ) );
  pts.push_back( IndexTriplet( 0,    1,    0    ) );
  pts.push_back( IndexTriplet( 0,    1,    0    ) );
  pts.push_back( IndexTriplet( nx-1, ny-1, nz-1 ) );
  pts.push_back( IndexTriplet( 0,    ny-1, 0    ) );
  pts.push_back( IndexTriplet( 0,    0,    nz-1 ) );
  pts.push_back( IndexTriplet( 2,    ny-1, 1    ) );

  SVolField field( get_n_tot<SVolField>(dim,bcx,bcy,bcz),
                   get_ghost_set<SVolField>(dim,bcx,bcy,bcz),
                   NULL );
  field = 1.0;

  bool isFailed = false;

  cout << "Testing simple BC usage ... " << flush;

  // apply bcs to the field using a "time varying" function.  We use one function to get the time 
  for( vector<IndexTriplet>::const_iterator ipt=pts.begin(); ipt!=pts.end(); ++ipt ){
    BoundaryCondition<SVolField,BCFun> bc( *ipt, dim, bcx, bcy, bcz, f );
    bc(field);
    // check:
    const int ix = get_index_with_ghost<SVolField>( dim, bcx, bcy, bcz, *ipt );
    isFailed = abs( field[ix] - f() ) > ATOL;
  }

  // apply constant-time BCs
  for( vector<IndexTriplet>::const_iterator ipt=pts.begin(); ipt!=pts.end(); ++ipt ){
    BoundaryCondition<SVolField,BCFun> bc( *ipt, dim, bcx, bcy, bcz, fnotime );
    bc(field);
    // check:
    const int ix = get_index_with_ghost<SVolField>( dim, bcx, bcy, bcz, *ipt );
    isFailed = abs( field[ix] - fnotime() ) > ATOL;
  }

  if( isFailed )
    cout << "***FAIL***" << endl << endl;
  else
    cout << "PASS" << endl << endl;
}

//====================================================================

int main()
{
  test1();
}

//====================================================================
