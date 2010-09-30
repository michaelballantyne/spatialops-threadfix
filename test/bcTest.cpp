#include <vector>
#include <cmath>

#include <spatialops/structured/FVStaggeredBCTools.h>
#include <spatialops/structured/FVStaggered.h>

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

bool test1()
{
  using namespace SpatialOps;
  using namespace structured;

  const double ATOL = 1.0e-5;
  const int nx=3, ny=4, nz=5;

  IntVec dim(nx,ny,nz);
  const bool bcx=true,  bcy=true, bcz=true;

  typedef boost::function< double() > BCFun;
  BCFun fnotime = &tmp;

  // composite function:  tmp2( ttmp() )
  BCFun tfun = boost::lambda::bind( &ttmp );
  BCFun f = boost::lambda::bind( &tmp2, boost::lambda::bind(&ttmp) );

  // pick some points to set bcs on.  Note that they need not actually
  // reside on a boundary.
  std::vector<IndexTriplet> pts;
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

  SVolField field( get_window_with_ghost<SVolField>(dim,bcx,bcy,bcz), NULL );
  field = 1.0;

  bool isFailed = false;

  std::cout << "Testing simple BC usage ... " << std::flush;

  // apply bcs to the field using a "time varying" function.  We use one function to get the time 
  for( std::vector<IndexTriplet>::const_iterator ipt=pts.begin(); ipt!=pts.end(); ++ipt ){
    BoundaryCondition<SVolField,BCFun> bc( *ipt, dim, bcx, bcy, bcz, f );
    bc(field);
    // check:
    const int ix = get_index_with_ghost<SVolField>( dim, bcx, bcy, bcz, *ipt );
    isFailed = std::abs( field[ix] - f() ) > ATOL;
  }

  // apply constant-time BCs
  for( std::vector<IndexTriplet>::const_iterator ipt=pts.begin(); ipt!=pts.end(); ++ipt ){
    BoundaryCondition<SVolField,BCFun> bc( *ipt, dim, bcx, bcy, bcz, fnotime );
    bc(field);
    // check:
    const int ix = get_index_with_ghost<SVolField>( dim, bcx, bcy, bcz, *ipt );
    isFailed = std::abs( field[ix] - fnotime() ) > ATOL;
  }

  if( isFailed )
    std::cout << "***FAIL***" << std::endl << std::endl;
  else
    std::cout << "PASS" << std::endl << std::endl;

  return isFailed;
}

//====================================================================

int main()
{
  return test1();
}

//====================================================================
