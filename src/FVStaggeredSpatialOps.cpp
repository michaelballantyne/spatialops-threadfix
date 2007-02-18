#include <FVStaggeredSpatialOps.h>

#include <sstream>
#include <stdexcept>
#include <numeric>
#include <algorithm>

using std::vector;

namespace SpatialOps{
namespace FVStaggeredUniform{

//====================================================================

//--------------------------------------------------------------------
int colcount( const vector<int> & dims,
	      const vector<int> & nghost )
{
  int nrows=1;
  vector<int>::const_iterator ig = nghost.begin();
  for( vector<int>::const_iterator ii=dims.begin(); ii!=dims.end(); ++ii ){
    if( *ii > 1 ){
      int nn = *ii;
      nn += *ig++;      // add "left"  side ghosts
      nn += *ig++;      // add "right" side ghosts
      nrows *= nn;
    }
  }
  return nrows;
}
//--------------------------------------------------------------------
int rowcount( const vector<int> & dims,
	      const vector<int> & nghost )
{
  return colcount(dims,nghost);
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
ScratchOperator::ScratchOperator( const vector<int> & dimExtent,
				  const vector<int> & nghostSrc,
				  const vector<int> & nghostDest,
				  const int entriesPerRow,
				  const Direction dir )
  : SpatialOperator( rowcount( dimExtent, nghostDest ),
		     colcount( dimExtent, nghostSrc  ),
		     nghostSrc,
		     nghostDest,
		     entriesPerRow,
		     dimExtent ),

    entriesPerRow_( entriesPerRow ),
    dir_( dir )
{
  assert( entriesPerRow==3 );

  ndim_ = 0;
  for( vector<int>::const_iterator i=extent_.begin(); i<extent_.end(); ++i ){
    if( *i > 1 ) ++ndim_;
  }

  setup_matrix();
  finalize();
}
//--------------------------------------------------------------------
ScratchOperator::~ScratchOperator()
{
}
//--------------------------------------------------------------------
void
ScratchOperator::setup_matrix()
{
  vector<double> vals( entriesPerRow_, 0.0 );
  vector<int>    ixs ( entriesPerRow_, 0   );

  for( int i=0; i<nrows(); ++i ){
    vals.clear();
    ixs.clear();
    get_row_entries ( i, vals, ixs );
    insert_row_entry( i, vals, ixs );
  }
}
//--------------------------------------------------------------------
void
ScratchOperator::get_row_entries( const int irow,
				  vector<double> & vals,
				  vector<int> & ixs ) const
{
  const double zero = 0.0;

  const int nx = extent_[0] + nghost_src(XDIM,MINUS) + nghost_src(XDIM,PLUS);
  const int ny = extent_[1] + nghost_src(YDIM,MINUS) + nghost_src(YDIM,PLUS);
  const int nz = extent_[2] + nghost_src(ZDIM,MINUS) + nghost_src(ZDIM,PLUS);

  switch( dir_ ){

  case X_DIR:{
  
    const int i = irow%nx;

    //  if( i >= nghost_dest(XDIM,MINUS) && i < nghost_dest(XDIM,PLUS) ){

      if( i==0 ){
	vals.push_back( zero );  ixs.push_back( irow   );
	vals.push_back( zero );  ixs.push_back( irow+1 );
	//      vals.push_back( zero );  ixs.push_back( irow+2 );
      }
      else if( i==nx-1 ){
	//      vals.push_back( 0.0 );  ixs.push_back( irow-2 );
	vals.push_back( zero );  ixs.push_back( irow-1 );
	vals.push_back( zero );  ixs.push_back( irow   );
      }
      else{
	vals.push_back( zero ); ixs.push_back( irow-1 );
	vals.push_back( zero ); ixs.push_back( irow   );
	vals.push_back( zero ); ixs.push_back( irow+1 );
      }
      //  }
    break;
  }

  case Y_DIR:{

    assert( ndim_ > 1 );

    const int j = (irow/nx)%ny;

    if( j==0 ){
      vals.push_back( zero );  ixs.push_back( irow );
      vals.push_back( zero );  ixs.push_back( irow+nx );
    }
    else if( j==ny-1 ){
      vals.push_back( zero );  ixs.push_back( irow-nx );
      vals.push_back( zero );  ixs.push_back( irow );
    }
    else{
      vals.push_back( zero );  ixs.push_back( irow-nx );
      vals.push_back( zero );  ixs.push_back( irow );
      vals.push_back( zero );  ixs.push_back( irow+nx );
    }
    break;
  }

  case Z_DIR:{

    assert( ndim_==3 );

    const int k = irow/(nx*ny);
    if( k==0 ){
      vals.push_back( zero );  ixs.push_back( irow );
      vals.push_back( zero );  ixs.push_back( irow+nx*ny );
    }
    else if( k==nz-1 ){
      vals.push_back( zero );  ixs.push_back( irow-nx*ny );
      vals.push_back( zero );  ixs.push_back( irow );
    }
    else{
      vals.push_back( zero );  ixs.push_back( irow-nx*ny );
      vals.push_back( zero );  ixs.push_back( irow );
      vals.push_back( zero );  ixs.push_back( irow+nx*ny );
    }

    break;
  }

  default:{
    std::ostringstream errmsg;
    errmsg << "ERROR: Invalid dimension for interpolation." << std::endl;
    throw std::runtime_error( errmsg.str() );
  }
  }
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
LinearInterpolant::LinearInterpolant( const vector<int> & dimExtent,
				      const vector<int> & nghostSrc,
				      const vector<int> & nghostDest,
				      const OpType opType,
				      const Direction dir )
  : SpatialOperator( rowcount( dimExtent, nghostDest ),
		     colcount( dimExtent, nghostSrc  ),
		     nghostSrc,
		     nghostDest,
		     entries_per_row(),
		     dimExtent ),
    opType_( opType ),
    dir_( dir )
{
  ndim_ = 0;
  for( vector<int>::const_iterator i=extent_.begin(); i<extent_.end(); ++i ){
    if( *i > 1 ) ++ndim_;
  }

  setup_matrix();
  finalize();
}
//--------------------------------------------------------------------
LinearInterpolant::~LinearInterpolant()
{
}
//--------------------------------------------------------------------
void
LinearInterpolant::setup_matrix()
{
  vector<double> vals( entries_per_row(), 0.0 );
  vector<int>    ixs ( entries_per_row(), 0   );

  for( int i=0; i<nrows(); ++i ){
    vals.clear();
    ixs.clear();
    get_row_entries ( i, vals, ixs );
    insert_row_entry( i, vals, ixs );
  }
}
//--------------------------------------------------------------------
void
LinearInterpolant::get_row_entries( const int irow,
				    vector<double> & vals,
				    vector<int> & ixs ) const
{
  // For cell to face interpolant, the interpolated variable is
  // staggered in the (-) direction.  For face to cell, we stagger the
  // interpolated value to the (+) direction.
  const bool staggerLeft = opType_==CellToFace;

  //
  // all nonzero entries are equal to 0.5 since we are doing midpoint
  // interpolation on a uniform mesh
  //
  const double val = 0.5;

  const int nxd = extent_[0] + nghost_dest(XDIM,MINUS) + nghost_dest(XDIM,PLUS);
  const int nyd = extent_[1] + nghost_dest(YDIM,MINUS) + nghost_dest(YDIM,PLUS);
//const int nzd = extent_[2] + nghost_dest(ZDIM,MINUS) + nghost_dest(ZDIM,PLUS);

  const int nxs = extent_[0] + nghost_src(XDIM,MINUS) + nghost_src(XDIM,PLUS);
  const int nys = extent_[1] + nghost_src(YDIM,MINUS) + nghost_src(YDIM,PLUS);
  const int nzs = extent_[2] + nghost_src(ZDIM,MINUS) + nghost_src(ZDIM,PLUS);

  // get the (ijk) index for the dest array
  const int idest = irow%nxd;
  const int jdest = irow/nxd % nyd;
  const int kdest = irow/(nxd*nyd);

  // get the (ijk) index for the src array
  const int isrc = idest - nghost_dest(XDIM,MINUS) + nghost_src(XDIM,MINUS);
  const int jsrc = jdest - nghost_dest(YDIM,MINUS) + nghost_src(YDIM,MINUS);
  const int ksrc = kdest - nghost_dest(ZDIM,MINUS) + nghost_src(ZDIM,MINUS);

  // are we in bounds?  If not, we don't have any entries...
  const int shift = (staggerLeft ? 0 :-1);
  if( isrc < 1+shift  ||  isrc >= nxs+shift ) return;
  if( extent_[1] > 1 )  if( jsrc < 1+shift  ||  jsrc >= nys+shift ) return;
  if( extent_[2] > 1 )  if( ksrc < 1+shift  ||  ksrc >= nzs+shift ) return;

  const int icol = ksrc*(nxs*nys) + jsrc*nxs + isrc;

  switch( dir_ ){

  case X_DIR:{

    const int ishift = (staggerLeft ? -1 : 1);

    vals.push_back( val );  ixs.push_back( icol        );
    vals.push_back( val );  ixs.push_back( icol+ishift );

    break;
  }

  case Y_DIR:{

    assert( ndim_ > 1 );

    const int jshift = (staggerLeft ? -1 : 1 ) * nxs;

    vals.push_back( val );  ixs.push_back( icol );
    vals.push_back( val );  ixs.push_back( icol+jshift );

    break;
  }

  case Z_DIR:{

    assert( ndim_==3 );

    const int kshift = (staggerLeft ? -1 : 1 ) * nxs*nys;
    vals.push_back( val );  ixs.push_back( icol        );
    vals.push_back( val );  ixs.push_back( icol+kshift );

    break;
  }

  default:{
    std::ostringstream errmsg;
    errmsg << "ERROR: Invalid dimension for interpolation." << std::endl;
    throw std::runtime_error( errmsg.str() );
  }
  }
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
Gradient2ndOrder::Gradient2ndOrder( const vector<double> & meshSpacing,
				    const vector<int> & dimExtent,
				    const vector<int> & nghostSrc,
				    const vector<int> & nghostDest,
				    const OpType opType,
				    const Direction dir )
  : SpatialOperator( rowcount( dimExtent, nghostDest),
		     colcount( dimExtent, nghostSrc  ),
		     nghostSrc,
		     nghostDest,
		     entries_per_row(),
		     dimExtent ),

    spacing_( meshSpacing ),
    opType_( opType ),
    dir_( dir )
{
  ndim_ = 0;
  for( vector<int>::const_iterator i=extent_.begin(); i<extent_.end(); ++i ){
    if( *i > 1 ) ++ndim_;
  }

  setup_matrix();
  finalize();
}
//--------------------------------------------------------------------
Gradient2ndOrder::~Gradient2ndOrder()
{
}
//--------------------------------------------------------------------
void
Gradient2ndOrder::setup_matrix()
{
  vector<double> vals( entries_per_row(), 0.0 );
  vector<int>    ixs ( entries_per_row(), 0   );

  for( int i=0; i<nrows(); ++i ){
    vals.clear();
    ixs.clear();
    get_row_entries ( i, vals, ixs );
    insert_row_entry( i, vals, ixs );
  }
}
//--------------------------------------------------------------------
void
Gradient2ndOrder::get_row_entries( const int irow,
				   vector<double> & vals,
				   vector<int> & ixs ) const
{
  const bool staggerLeft = opType_==CellToFace;

  const int nxd = extent_[0] + nghost_dest(XDIM,MINUS) + nghost_dest(XDIM,PLUS);
  const int nyd = extent_[1] + nghost_dest(YDIM,MINUS) + nghost_dest(YDIM,PLUS);
//const int nzd = extent_[2] + nghost_dest(ZDIM,MINUS) + nghost_dest(ZDIM,PLUS);

  const int nxs = extent_[0] + nghost_src(XDIM,MINUS) + nghost_src(XDIM,PLUS);
  const int nys = extent_[1] + nghost_src(YDIM,MINUS) + nghost_src(YDIM,PLUS);
  const int nzs = extent_[2] + nghost_src(ZDIM,MINUS) + nghost_src(ZDIM,PLUS);

  // get the (ijk) index for the dest array
  const int idest = irow%nxd;
  const int jdest = irow/nxd % nyd;
  const int kdest = irow/(nxd*nyd);

  // get the (ijk) index for the src array
  int isrc = idest - nghost_dest(XDIM,MINUS) + nghost_src(XDIM,MINUS);
  int jsrc = jdest - nghost_dest(YDIM,MINUS) + nghost_src(YDIM,MINUS);
  int ksrc = kdest - nghost_dest(ZDIM,MINUS) + nghost_src(ZDIM,MINUS);

  // are we in bounds?  If not, we don't have any entries...
  const int shift = (staggerLeft ? 0 : -1);
  bool inBoundsLo[3], inBoundsHi[3];
  for(int i=0; i<3; ++i){ inBoundsLo[i]=inBoundsHi[i]=true; }

  if( isrc <    1+shift ){ inBoundsLo[0]=false; isrc=0;     }
  if( isrc >= nxs+shift ){ inBoundsHi[0]=false; isrc=nxs-1; }
  if( extent_[1] > 1 ){
    if( jsrc <    1+shift ){ inBoundsLo[1]=false; jsrc=0;     }
    if( jsrc >= nys+shift ){ inBoundsHi[1]=false; jsrc=nys-1; }
  }
  if( extent_[2] > 1 ){
    if( ksrc <    1+shift ){ inBoundsLo[2]=false; ksrc=0;     }
    if( ksrc >= nzs+shift ){ inBoundsHi[2]=false; ksrc=nzs-1; }
  }

  const int icol = ksrc*(nxs*nys) + jsrc*nxs + isrc;

  int ifac = (staggerLeft ? -1 : 1);

  bool inBounds[3];
  for( int i=0; i<3; ++i )   inBounds[i] = ( inBoundsLo[i] && inBoundsHi[i] );
    
  switch( dir_ ){

  case X_DIR:{

    const double fac = inBounds[0] ? 1.0/spacing_[0] : 0.0;

    if( !inBoundsLo[0] ) ifac=1;
    if( !inBoundsHi[0] ) ifac=-1;
    const int ishift = ifac;
    vals.push_back( -ifac*fac );  ixs.push_back( icol        );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+ishift );

    break;
  }

  case Y_DIR:{

    assert( ndim_ > 1 );

    const double fac = inBounds[1] ? 1.0/spacing_[1] : 0.0;

    if( !inBoundsLo[1] ) ifac=1;
    if( !inBoundsHi[1] ) ifac=-1;
    const int jshift = ifac * nxs;
    vals.push_back( -ifac*fac );  ixs.push_back( icol );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+jshift );

    break;
  }

  case Z_DIR:{
    assert( ndim_ == 3 );

    const double fac = inBounds[2] ? 1.0/spacing_[2] : 0.0;

    if( !inBoundsLo[2] ) ifac=1;
    if( !inBoundsHi[2] ) ifac=-1;
    const int kshift = ifac * nxs*nys;
    vals.push_back( -ifac*fac );  ixs.push_back( icol        );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+kshift );

    break;
  }

  default:{
    std::ostringstream errmsg;
    errmsg << "ERROR: Incompatible dimension for gradient operator." << std::endl;
    throw std::runtime_error( errmsg.str() );
  }
  } // end switch
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
Divergence2ndOrder::Divergence2ndOrder( const vector<double> & cellFaceArea,
					const double cellVolume,
					const vector<int> & dimExtent,
					const vector<int> & nghostSrc,
					const vector<int> & nghostDest,
					const OpType opType,
					const Direction dir )
  : SpatialOperator( rowcount( dimExtent, nghostDest ),
		     colcount( dimExtent, nghostSrc  ),
		     nghostSrc,
		     nghostDest,
		     entries_per_row(), 
		     dimExtent ),

    faceArea_( cellFaceArea ),
    cellVol_( cellVolume ),
    opType_( opType ),
    dir_( dir )
{
  ndim_ = 0;
  for( vector<int>::const_iterator i=extent_.begin(); i<extent_.end(); ++i ){
    if( *i > 1 ) ++ndim_;
  }

  setup_matrix();
  finalize();
}
//--------------------------------------------------------------------
Divergence2ndOrder::~Divergence2ndOrder()
{
}
//--------------------------------------------------------------------
void
Divergence2ndOrder::setup_matrix()
{
  vector<double> vals( entries_per_row(), 0.0 );
  vector<int>    ixs ( entries_per_row(), 0   );

  for( int i=0; i<nrows(); ++i ){
    vals.clear();
    ixs.clear();
    get_row_entries ( i, vals, ixs );
    insert_row_entry( i, vals, ixs );
  }
}
//--------------------------------------------------------------------
void
Divergence2ndOrder::get_row_entries( const int irow,
				     vector<double> & vals,
				     vector<int> & ixs ) const
{
  const bool staggerLeft = (opType_ == CellToFace);

  const int nxd = extent_[0] + nghost_dest(XDIM,MINUS) + nghost_dest(XDIM,PLUS);
  const int nyd = extent_[1] + nghost_dest(YDIM,MINUS) + nghost_dest(YDIM,PLUS);
//const int nzd = extent_[2] + nghost_dest(ZDIM,MINUS) + nghost_dest(ZDIM,PLUS);

  const int nxs = extent_[0] + nghost_src(XDIM,MINUS) + nghost_src(XDIM,PLUS);
  const int nys = extent_[1] + nghost_src(YDIM,MINUS) + nghost_src(YDIM,PLUS);
  const int nzs = extent_[2] + nghost_src(ZDIM,MINUS) + nghost_src(ZDIM,PLUS);

  // get the (ijk) index for the dest array
  const int idest = irow%nxd;
  const int jdest = irow/nxd % nyd;
  const int kdest = irow/(nxd*nyd);

  // get the (ijk) index for the src array
  int isrc = idest - nghost_dest(XDIM,MINUS) + nghost_src(XDIM,MINUS);
  int jsrc = jdest - nghost_dest(YDIM,MINUS) + nghost_src(YDIM,MINUS);
  int ksrc = kdest - nghost_dest(ZDIM,MINUS) + nghost_src(ZDIM,MINUS);

  // are we in bounds?  If not, we don't have any entries...
  const int shift = (staggerLeft ? 0 : -1);
  bool inBoundsLo[3], inBoundsHi[3];
  for(int i=0; i<3; ++i){ inBoundsLo[i]=inBoundsHi[i]=true; }
  if( isrc <    1+shift ){ inBoundsLo[0]=false; isrc=0;     }
  if( isrc >= nxs+shift ){ inBoundsHi[0]=false; isrc=nxs-1; }
  if( extent_[1] > 1+shift ){
    if( jsrc <    1+shift ){ inBoundsLo[1]=false; jsrc=0;     }
    if( jsrc >= nys+shift ){ inBoundsHi[1]=false; jsrc=nys-1; }
  }
  if( extent_[2] > 1+shift ){
    if( ksrc <    1+shift ){ inBoundsLo[2]=false; ksrc=0;     }
    if( ksrc >= nzs+shift ){ inBoundsHi[2]=false; ksrc=nzs-1; }
  }

  const int icol = ksrc*(nxs*nys) + jsrc*nxs + isrc;

  int ifac = (staggerLeft ? -1 : 1);

  bool inBounds[3];
  for( int i=0; i<3; ++i )   inBounds[i] = ( inBoundsLo[i] && inBoundsHi[i] );

  switch( dir_ ){

  case X_DIR:{

    const double fac = inBounds[0] ? faceArea_[0]/cellVol_ : 0.0;

    if( !inBoundsLo[0] ) ifac=1;
    if( !inBoundsHi[0] ) ifac=-1;
    const int ishift = ifac;

    vals.push_back( -ifac*fac );  ixs.push_back( icol        );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+ishift );

    break;
  }

  case Y_DIR:{

    assert( ndim_ > 1 );

    const double fac = inBounds[1] ? faceArea_[1]/cellVol_ : 0.0;

    if( !inBoundsLo[1] ) ifac=1;
    if( !inBoundsHi[1] ) ifac=-1;
    const int jshift = ifac * nxs;
    vals.push_back( -ifac*fac );  ixs.push_back( icol        );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+jshift );
    break;
  }

  case Z_DIR:{
    assert( ndim_ == 3 );

    const double fac = inBounds[2] ? faceArea_[2]/cellVol_ : 0.0;

    if( !inBoundsLo[2] ) ifac=1;
    if( !inBoundsHi[2] ) ifac=-1;
    const int kshift = ifac * nxs*nys;
    vals.push_back( -ifac*fac );  ixs.push_back( icol        );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+kshift );

    break;
  }

  default:{
    std::ostringstream errmsg;
    errmsg << "ERROR: Incompatible dimension for gradient operator." << std::endl;
    throw std::runtime_error( errmsg.str() );
  }
  } // end switch

}
//--------------------------------------------------------------------


//====================================================================


} // namespace FVStaggeredUniform
} // namespace SpatialOps
