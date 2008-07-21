#include <FVOneDimensional.h>

#include <boost/static_assert.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace SpatialOps{
namespace FVStaggered{

//--------------------------------------------------------------------

OneDimInterpolantAssembler::
OneDimInterpolantAssembler( const int polynomialOrder,
                            const std::vector<double>& x,
                            const std::vector<double>& xs )
  : polyOrder_( polynomialOrder ),
    xsurf_( get_xs(x,xs) ),
    lagrange_( get_x(x,xs) ),
    nx_( x.size() )
{
  BOOST_STATIC_ASSERT(   SVolField::Ghost::NM == 1 );
  BOOST_STATIC_ASSERT(   SVolField::Ghost::NP == 1 );
  BOOST_STATIC_ASSERT( SSurfXField::Ghost::NM == 1 );
  BOOST_STATIC_ASSERT( SSurfXField::Ghost::NP == 1 );

  // determine how many nonzeros we have.  In general we require one
  // more point than the polynomial order.  In the case of a uniform
  // mesh where xs is at x midpoints, we require a smaller stencil.
  numNonzero_ = polynomialOrder+1;

  bool isReduced = true;
  double dx=x[1]-x[0];
  double dxold=dx;
  typedef std::vector<double> DVec;
  for( DVec::const_iterator ix=x.begin()+1; ix!=x.end(); ++ix ){
    dx = *ix - *(ix-1);
    if( std::fabs(dx-dxold)>1.0e-3*dx ) isReduced = false;
    dxold = dx;
  }
  dx = xs[1]-xs[0];
  dxold = dx;
  for( DVec::const_iterator ix=xs.begin()+1; ix!=xs.end(); ++ix ){
    dx = *ix - *(ix-1);
    if( std::fabs(dx-dxold)>1.0e-3*dx ) isReduced = false;
    dxold = dx;
  }
  // ensure that surface (xs) is at midpoint of volume (x)
  {
    const double xmid = 0.5*(x[0]+x[1]);
    if( std::fabs( xmid-xs[1] ) > 1.0e-8 ) isReduced = false;
  }

  if( isReduced ){
    --numNonzero_;
    --polyOrder_;
    //cout << "reducing interpolant stencil because we are on a uniform mesh" << endl;
  }
}

//--------------------------------------------------------------------

int
OneDimInterpolantAssembler::
get_ncols() const
{
  return nx_+2;
}

//--------------------------------------------------------------------

int
OneDimInterpolantAssembler::
get_nrows() const
{
  return nx_+3;
}

//--------------------------------------------------------------------

void
OneDimInterpolantAssembler::
get_ghost_cols( std::set<int>& ghostCols ) const
{
  ghostCols.clear();
  ghostCols.insert( 0 );
  ghostCols.insert( nx_+1 );
}

//--------------------------------------------------------------------

void
OneDimInterpolantAssembler::
get_ghost_rows( std::set<int>& ghostRows ) const
{
  ghostRows.clear();
  ghostRows.insert( 0 );
  ghostRows.insert( nx_+2 );
}

//--------------------------------------------------------------------

void
OneDimInterpolantAssembler::
get_row_entries( const int irow,
                 std::vector<double> & vals,
                 std::vector<int> & ixs ) const
{
  // irow -> index for destination field that we want.
  lagrange_.get_interp_coefs_indices( xsurf_[irow], polyOrder_, vals, ixs );
//   cout << "row: " << irow << endl
//        << setw(10) << " ix" << setw(10) << "coefs" << setw(10) << "x" << endl
//        << "-----------------------------------------" << endl;
//   std::vector<double>::const_iterator v=vals.begin();
//   std::vector<int   >::const_iterator i=ixs.begin();
//   for( ; i!=ixs.end(); ++i, ++v )
//     cout << setw(10) << *i << setw(10) << *v << setw(10) << lagrange_.get_x()[*i] << endl;
//   cout << xsurf_[irow] << endl
//        << endl;
}

//--------------------------------------------------------------------

std::vector<double>
OneDimInterpolantAssembler::
get_x( const std::vector<double>& x,
       const std::vector<double>& xs )
{
  typedef std::vector<double> DVec;
  DVec xv( x.size()+2, 0.0 );
  std::copy( x.begin(), x.end(), xv.begin()+1 );

  xv[0] = 2*xs[0] - x[0];

  DVec::iterator ixv = xv.end()-1;
  DVec::const_iterator ix=x.end()-1, ixs=xs.end()-1;
  *ixv = 2*(*ixs) - *ix;

  return xv;
}

//--------------------------------------------------------------------

std::vector<double>
OneDimInterpolantAssembler::
get_xs( const std::vector<double>& x,
        const std::vector<double>& xs )
{
  typedef std::vector<double> DVec;
  DVec xsurf( xs.size()+2, 0.0 );

  DVec::iterator ixsurf = xsurf.begin();
  DVec::const_iterator ix=x.begin(), ixs=xs.begin();

  // left side
  *ixsurf = 3*(*ixs)-2*(*ix);

  // interior
  for( ixsurf=xsurf.begin()+1, ixs=xs.begin(); ixs!=xs.end(); ++ixsurf, ++ixs )
    *ixsurf = *ixs;

  // right side
  ix = x.end()-1;
  ixs = xs.end()-1;
  *ixsurf = 3*(*ixs) - 2*(*ix);

  return xsurf;
}

//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------


OneDimGradientAssembler::
OneDimGradientAssembler( const int polynomialOrder,
                         const std::vector<double>& x,
                         const std::vector<double>& xs )
  : polyOrder_( polynomialOrder ),
    xsurf_( get_xs(x,xs) ),
    lagrange_( get_x(x,xs) ),
    nx_( x.size() )
{
  BOOST_STATIC_ASSERT(   SVolField::Ghost::NM == 1 );
  BOOST_STATIC_ASSERT(   SVolField::Ghost::NP == 1 );
  BOOST_STATIC_ASSERT( SSurfXField::Ghost::NM == 1 );
  BOOST_STATIC_ASSERT( SSurfXField::Ghost::NP == 1 );

  // determine how many nonzeros we have.  In general we require one
  // more point than the polynomial order.  In the case of a uniform
  // mesh where xs is at x midpoints, we require a smaller stencil.
  numNonzero_ = polynomialOrder+1;

  bool isReduced = true;
  double dx=x[1]-x[0];
  double dxold=dx;
  typedef std::vector<double> DVec;
  for( DVec::const_iterator ix=x.begin()+1; ix!=x.end(); ++ix ){
    dx = *ix - *(ix-1);
    if( std::fabs(dx-dxold)>1.0e-3*dx ) isReduced = false;
    dxold = dx;
  }
  dx = xs[1]-xs[0];
  dxold = dx;
  for( DVec::const_iterator ix=xs.begin()+1; ix!=xs.end(); ++ix ){
    dx = *ix - *(ix-1);
    if( std::fabs(dx-dxold)>1.0e-3*dx ) isReduced = false;
    dxold = dx;
  }
  // ensure that surface (xs) is at midpoint of volume (x)
  {
    const double xmid = 0.5*(x[0]+x[1]);
    if( std::fabs( xmid-xs[1] ) > 1.0e-8 ) isReduced = false;
  }

//   if( isReduced ){
//     --numNonzero_;
//     --polyOrder_;
//     cout << "reducing gradient stencil because we are on a uniform mesh" << endl;
//   }
}

//--------------------------------------------------------------------

int
OneDimGradientAssembler::
get_ncols() const
{
  return nx_+2;
}

//--------------------------------------------------------------------

int
OneDimGradientAssembler::
get_nrows() const
{
  return nx_+3;
}

//--------------------------------------------------------------------

void
OneDimGradientAssembler::
get_ghost_cols( std::set<int>& ghostCols ) const
{
  ghostCols.clear();
  ghostCols.insert( 0 );
  ghostCols.insert( nx_+1 );
}

//--------------------------------------------------------------------

void
OneDimGradientAssembler::
get_ghost_rows( std::set<int>& ghostRows ) const
{
  ghostRows.clear();
  ghostRows.insert( 0 );
  ghostRows.insert( nx_+2 );
}

//--------------------------------------------------------------------

void
OneDimGradientAssembler::
get_row_entries( const int irow,
                 std::vector<double> & vals,
                 std::vector<int> & ixs ) const
{
  // irow -> index for destination field that we want.
  lagrange_.get_derivative_coefs_indices( xsurf_[irow], polyOrder_, vals, ixs );
}

//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------

OneDimDivergenceAssembler::
OneDimDivergenceAssembler( const std::vector<double>& x,
                           const std::vector<double>& xs )
  : nx_( x.size() ),
    numNonzero_( 2 )
{
  // ensure that volumes are at surface midpoints
  bool isFailed = false;
  std::vector<double>::const_iterator ix=x.begin(), ixs=xs.begin();
  for( ; ix!=x.end(); ++ix, ++ixs ){
    const double xmid = 0.5*( *(ixs+1) + *ixs );
    if( std::fabs( xmid-*ix ) > 1.0e-15 ){
      isFailed = true;
      // cout << *ix << ", " << xmid << endl;
    }
  }
  if( isFailed ){
    std::ostringstream msg;
    msg << "Error in forming divergence operator for 1D mesh." << endl
        << "This operator assumes that volume variables are at" << endl
        << "the cell centroid.  That assumption has been violated" << endl
        << "by the given mesh." << endl;
    throw std::runtime_error( msg.str() );
  }

  // get all surface locations and then generate coefficients.
  dxinv_.resize(nx_+2);
  const std::vector<double> xsurf( OneDimInterpolantAssembler::get_xs(x,xs) );
  ixs = xsurf.begin();
  std::vector<double>::iterator idx=dxinv_.begin();
  for( ; idx!=dxinv_.end(); ++ixs, ++idx ){
    *idx = 1.0/( *(ixs+1)-*ixs );
  }
}

//--------------------------------------------------------------------

void
OneDimDivergenceAssembler::
get_row_entries( const int irow,
                 std::vector<double> & vals,
                 std::vector<int> & ixs ) const
{
  const double coef = dxinv_[irow];
  vals.push_back(-coef);
  vals.push_back( coef);
  ixs.push_back(irow);
  ixs.push_back(irow+1);
}

//--------------------------------------------------------------------

void
OneDimDivergenceAssembler::
get_ghost_cols( std::set<int>& ghostCols ) const
{
  ghostCols.clear();
  ghostCols.insert( 0 );
  ghostCols.insert( nx_+2 );
}

//--------------------------------------------------------------------

void
OneDimDivergenceAssembler::
get_ghost_rows( std::set<int>& ghostRows ) const
{
  ghostRows.clear();
  ghostRows.insert( 0 );
  ghostRows.insert( nx_+1 );
}

//--------------------------------------------------------------------

} // namespace SpatialOps
} // namespace FVStaggered
