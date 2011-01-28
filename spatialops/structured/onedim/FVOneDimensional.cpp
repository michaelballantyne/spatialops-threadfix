#include "FVOneDimensional.h"

#include <boost/static_assert.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace SpatialOps{
namespace structured{


//--------------------------------------------------------------------

std::vector<double>
get_x_src( const int nGhost,
           const std::vector<double>& xsrc,
           const std::vector<double>& xdest )
{
  typedef std::vector<double> DVec;
  assert( nGhost == 1 );
  DVec x( xsrc.size()+2*nGhost, 0.0 );
  std::copy( xsrc.begin(), xsrc.end(), x.begin()+1 );

  assert( nGhost == 1 );

  if( &xdest[0] == &xsrc[0] ){
    x[0] = xsrc[0] - (xsrc[1]-xsrc[0]);
    const double xL = *(xsrc.end()-1);
    const double xnm1 = *(xsrc.end()-2);
    *(x.end()-1) = xL + (xL-xnm1);
    return x;
  }

  const bool xIsSurf = ( xdest[0] > xsrc[0] );

  const DVec& xcell = xIsSurf ? xdest : xsrc;
  const DVec& xsurf = xIsSurf ? xsrc : xdest;

  const double xc0 = 2*xsurf[0] - xcell[0];
  const double xs0 = 2*xc0 - xsurf[0];
  const double xcL = 2**(xsurf.end()-1) - *(xcell.end()-1);
  const double xsL = 2*xcL - *(xsurf.end()-1);

  DVec::iterator ix0 = x.begin();
  DVec::iterator ixL = x.end()-1;

  if( xIsSurf ){
    *ix0 = xs0;
    *ixL = xsL;
  }
  else{
    *ix0 = xc0;
    *ixL = xcL;
  }

  return x;
}

//--------------------------------------------------------------------

std::vector<double>
get_x_dest( const int nGhost,
            const std::vector<double>& xsrc,
            const std::vector<double>& xdest )
{
  typedef std::vector<double> DVec;
  assert( nGhost == 1 );
  DVec x( xdest.size()+2*nGhost, 0.0 );
  std::copy( xdest.begin(), xdest.end(), x.begin()+nGhost );

  assert( nGhost == 1 );

  if( &xdest[0] == &xsrc[0] ){
    x[0] = xsrc[0] - (xsrc[1]-xsrc[0]);
    const double xL = *(xsrc.end()-1);
    const double xnm1 = *(xsrc.end()-2);
    *(x.end()-1) = xL + (xL-xnm1);
    return x;
  }

  const bool xIsSurf = ( xdest[0] < xsrc[0] );

  const DVec& xcell = xIsSurf ? xsrc : xdest;
  const DVec& xsurf = xIsSurf ? xdest : xsrc;

  const double xc0 = 2*xsurf[0] - xcell[0];
  const double xs0 = 2*xc0 - xsurf[0];
  const double xcL = 2**(xsurf.end()-1) - *(xcell.end()-1);
  const double xsL = 2*xcL - *(xsurf.end()-1);
  
  DVec::iterator ix0 = x.begin();
  DVec::iterator ixL = x.end()-1;

  if( xIsSurf ){
    *ix0 = xs0;
    *ixL = xsL;
  }
  else{
    *ix0 = xc0;
    *ixL = xcL;
  }

  return x;
}

//--------------------------------------------------------------------

OneDimInterpolantAssembler::
OneDimInterpolantAssembler( const int polynomialOrder,
                            const int nGhost,
                            const std::vector<double>& xsrc,
                            const std::vector<double>& xdest )
  : ncol_(  xsrc.size() + 2*nGhost ),
    nrow_( xdest.size() + 2*nGhost ),
    polyOrder_( polynomialOrder ),
    xdest_( get_x_dest(nGhost,xsrc,xdest) ),
    lagrange_( get_x_src(nGhost,xsrc,xdest) ),
    nx_( xsrc.size()>xdest.size() ? xdest.size() : xsrc.size() )
{
  // determine how many nonzeros we have.  In general we require one
  // more point than the polynomial order.  In the case of a uniform
  // mesh where xs is at x midpoints, we require a smaller stencil.
  numNonzero_ = polynomialOrder+1;

//   bool isReduced = true;
//   double dx=xsrc[1]-xsrc[0];
//   double dxold=dx;
//   const double TOL = 1e-5;
//   typedef std::vector<double> DVec;
//   for( DVec::const_iterator ix=xsrc.begin()+1; ix!=xsrc.end(); ++ix ){
//     dx = *ix - *(ix-1);
//     if( std::fabs(dx-dxold)>TOL*dx ){
//       isReduced = false;
//       break;
//     }
//     dxold = dx;
//   }
//   if( isReduced ){
//     dx = xdest[1]-xdest[0];
//     dxold = dx;
//     for( DVec::const_iterator ix=xdest.begin()+1; ix!=xdest.end(); ++ix ){
//       dx = *ix - *(ix-1);
//       if( std::fabs(dx-dxold)>TOL*dx ){
//         isReduced = false;
//         break;
//       }
//       dxold = dx;
//     }
//   }

//   const double xmid = 0.5*(xsrc[0]+xsrc[1]);
//   if( xsrc[0]<xdest[0] )
//     if( std::fabs( xmid-xdest[0] ) > 1.0e-8 ) isReduced = false;
//   else
//     if( std::fabs( xmid-xdest[1] ) > 1.0e-8 ) isReduced = false;

//   if( isReduced ){
//     --numNonzero_;
//     --polyOrder_;
// //     cout << "reducing interpolant stencil because we are on a uniform mesh" << endl;
//   }
}

//--------------------------------------------------------------------

int
OneDimInterpolantAssembler::
get_ncols() const
{
  return ncol_;
}

//--------------------------------------------------------------------

int
OneDimInterpolantAssembler::
get_nrows() const
{
  return nrow_;
}

//--------------------------------------------------------------------

void
OneDimInterpolantAssembler::
get_ghost_cols( std::set<size_t>& ghostCols ) const
{
  ghostCols.clear();
  ghostCols.insert( 0 );
  ghostCols.insert( nx_+1 );
}

//--------------------------------------------------------------------

void
OneDimInterpolantAssembler::
get_ghost_rows( std::set<size_t>& ghostRows ) const
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
  if( size_t(irow) >= xdest_.size() ){
    std::cout << "problems: " << irow << ", " << xdest_.size() << std::endl;
    assert( size_t(irow) < xdest_.size() );
  }
  lagrange_.get_interp_coefs_indices( xdest_[irow], polyOrder_, vals, ixs );
//   cout << "row: " << irow << endl
//        << setw(10) << " ix" << setw(10) << "coefs" << setw(10) << "x" << endl
//        << "-----------------------------------------" << endl;
//   std::vector<double>::const_iterator v=vals.begin();
//   std::vector<int   >::const_iterator i=ixs.begin();
//   for( ; i!=ixs.end(); ++i, ++v )
//     cout << setw(10) << *i << setw(10) << *v << setw(10) << lagrange_.get_x()[*i] << endl;
//   cout << xdest_[irow] << endl
//        << endl;
}

//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------


OneDimGradientAssembler::
OneDimGradientAssembler( const int polynomialOrder,
                         const int nGhost,
                         const std::vector<double>& xsrc,
                         const std::vector<double>& xdest )
  : ncol_(  xsrc.size() + 2*nGhost ),
    nrow_( xdest.size() + 2*nGhost ),
    polyOrder_( polynomialOrder ),
    xdest_( get_x_dest(nGhost,xsrc,xdest) ),
    lagrange_( get_x_src(nGhost,xsrc,xdest) ),
    nx_( xsrc.size()>xdest.size() ? xdest.size() : xsrc.size() )
{
  BOOST_STATIC_ASSERT(   SVolField::Ghost::NGHOST == 1 );
  BOOST_STATIC_ASSERT(   SVolField::Ghost::NGHOST == 1 );
  BOOST_STATIC_ASSERT( SSurfXField::Ghost::NGHOST == 1 );
  BOOST_STATIC_ASSERT( SSurfXField::Ghost::NGHOST == 1 );

  // determine how many nonzeros we have.  In general we require one
  // more point than the polynomial order.  In the case of a uniform
  // mesh where xs is at x midpoints, we require a smaller stencil.
  numNonzero_ = polynomialOrder+1;

//   bool isReduced = true;
//   double dx=xsrc[1]-xsrc[0];
//   double dxold=dx;
//   typedef std::vector<double> DVec;
//   for( DVec::const_iterator ix=xsrc.begin()+1; ix!=xsrc.end(); ++ix ){
//     dx = *ix - *(ix-1);
//     if( std::fabs(dx-dxold)>1.0e-3*dx ) isReduced = false;
//     dxold = dx;
//   }
//   dx = xdest[1]-xdest[0];
//   dxold = dx;
//   for( DVec::const_iterator ix=xdest.begin()+1; ix!=xdest.end(); ++ix ){
//     dx = *ix - *(ix-1);
//     if( std::fabs(dx-dxold)>1.0e-3*dx ) isReduced = false;
//     dxold = dx;
//   }

//   const double xmid = 0.5*(xsrc[0]+xsrc[1]);
//   if( xsrc[0]<xdest[0] )
//     if( std::fabs( xmid-xdest[0] ) > 1.0e-8 ) isReduced = false;
//   else
//     if( std::fabs( xmid-xdest[1] ) > 1.0e-8 ) isReduced = false;

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
  return ncol_;
}

//--------------------------------------------------------------------

int
OneDimGradientAssembler::
get_nrows() const
{
  return nrow_;
}

//--------------------------------------------------------------------

void
OneDimGradientAssembler::
get_ghost_cols( std::set<size_t>& ghostCols ) const
{
  ghostCols.clear();
  ghostCols.insert( 0 );
  ghostCols.insert( nx_+1 );
}

//--------------------------------------------------------------------

void
OneDimGradientAssembler::
get_ghost_rows( std::set<size_t>& ghostRows ) const
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
  assert( size_t(irow) < xdest_.size() );
  lagrange_.get_derivative_coefs_indices( xdest_[irow], polyOrder_, vals, ixs );
}

//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------

OneDimDivergenceAssembler::
OneDimDivergenceAssembler( const int nGhost,
                           const std::vector<double>& xsrc,
                           const std::vector<double>& xdest )
  : ncol_( xdest.size() + 2*nGhost ),
    nrow_(  xsrc.size() + 2*nGhost ),
    nx_( xsrc.size() ),
    numNonzero_( 2 )
{
  // ensure that volumes are at surface midpoints
  bool isFailed = false;
  std::vector<double>::const_iterator ixs=xsrc.begin(), ixd=xdest.begin();
  for( ; ixs!=xsrc.end(); ++ixs, ++ixd ){
    const double xmid = 0.5*( *(ixd+1) + *ixd );
    if( std::fabs( xmid-*ixs ) > 1.0e-15 ){
      isFailed = true;
      // cout << *ix << ", " << xmid << endl;
    }
  }
  if( isFailed ){
    std::ostringstream msg;
    msg << "Error in forming divergence operator for 1D mesh." << std::endl
        << "This operator assumes that volume variables are at" << std::endl
        << "the cell centroid.  That assumption has been violated" << std::endl
        << "by the given mesh." << std::endl;
    throw std::runtime_error( msg.str() );
  }

  // get all surface locations and then generate coefficients.
  dxinv_.resize(nx_+2*nGhost);
  const std::vector<double> xd( get_x_dest(nGhost,xsrc,xdest) );
  std::vector<double>::const_iterator ix=xd.begin();
  std::vector<double>::iterator idx=dxinv_.begin();
  for( ; idx!=dxinv_.end(); ++ix, ++idx ){
    *idx = 1.0/( *(ix+1)-*ix );
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
get_ghost_cols( std::set<size_t>& ghostCols ) const
{
  ghostCols.clear();
  ghostCols.insert( 0 );
  ghostCols.insert( nx_+2 );
}

//--------------------------------------------------------------------

void
OneDimDivergenceAssembler::
get_ghost_rows( std::set<size_t>& ghostRows ) const
{
  ghostRows.clear();
  ghostRows.insert( 0 );
  ghostRows.insert( nx_+1 );
}

//--------------------------------------------------------------------

} // namespace structured
} // namespace SpatialOps
