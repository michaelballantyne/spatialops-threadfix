#include <numeric>

#include <LinearSystem.h>

#include <SpatialOperator.h>
#include <SpatialField.h>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

using std::vector;


//====================================================================


//--------------------------------------------------------------------
RHS::RHS( const std::vector<int> & domainExtent )
  : extent_( domainExtent ),
    npts_( std::accumulate(extent_.begin(), extent_.end(), 1, std::multiplies<int>() ) )
{
  field_.resize( npts_ );
}
//--------------------------------------------------------------------
RHS::~RHS()
{
}
//--------------------------------------------------------------------
void
RHS::reset( const double val )
{
  field_.assign( field_.size(), val );
}
//--------------------------------------------------------------------
void
RHS::reset_value( const int rownum, const double val )
{
  field_[rownum] = val;
}
//--------------------------------------------------------------------
void
RHS::add_contribution( const SpatialOps::SpatialField & f,
		       const double scaleFac )
{
  assert( consistency_check(f) );

  /*
   *  Procedure:
   *    1. Determine how many ghost cells are in the field
   *    2. Sum field elements into RHS, skipping ghost values
   */

  const int ng = f.nghost();
  const double * const fptr = f.get_ptr();

  const int nx=extent_[0];
  const int ny=extent_[1];
  const int nz=extent_[2];

  const vector<int> & fextent = f.get_extent();
  const int nxf = ( fextent[0] > 1 ) ? fextent[0]+2*ng : 1;
  const int nyf = ( fextent[1] > 1 ) ? fextent[1]+2*ng : 1;
  const int nzf = ( fextent[2] > 1 ) ? fextent[2]+2*ng : 1;
  const int ishift = nxf>1 ? ng : 0;
  const int jshift = nyf>1 ? ng : 0;
  const int kshift = nzf>1 ? ng : 0;

  for( int k=0; k<nz; ++k ){
    int ix  = k*nx*ny;            // index for the local field
    int ixf = (k+kshift)*nxf*nyf; // index for field we are adding (skip ghost cells)
    for( int j=0; j<ny; ++j ){
      ix  += j*nx;
      ixf += (j+jshift)*nxf + ishift;
      for( int i=0; i<extent_[0]; ++i ){
	//	ixf = (k+kshift)*nxf*nyf + (j+jshift)*nxf + (i+ishift);
	field_[ix] += scaleFac * fptr[ixf];
	++ix;
	++ixf;
      }
    }
  }
}
//--------------------------------------------------------------------
bool
RHS::consistency_check( const SpatialOps::SpatialField& f ) const
{
  bool ok = f.get_extent() == extent_;
  return ok;
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
LHS::LHS( const vector<int> & extent,
	  Epetra_CrsMatrix& A )
  : A_( A ),
    extent_( extent ),
    nrows_( A.NumMyRows() ),
    ncols_( A.NumMyCols() )
{
  rowDWork_.resize( A.NumMyNonzeros() );
  rowIWork_.resize( A.NumMyNonzeros() );
}
//--------------------------------------------------------------------
LHS::~LHS()
{
}
//--------------------------------------------------------------------
void
LHS::reset( const double val )
{
  A_.PutScalar(val);
}
//--------------------------------------------------------------------
void
LHS::add_contribution( const SpatialOps::SpatialOperator & op,
		       const double scaleFac )
{
  assert( compatibility_check(op) );

  // add non-ghost elements of the local matrix to this LHS operator
  const int nx=extent_[0];
  const int ny=extent_[1];
  const int nz=extent_[2];

  const int ng  = op.nghost();
  const vector<int> & opextent = op.get_extent();
  const int nxm = (opextent[0]>1) ? opextent[0]+2*ng : 1;
  const int nym = (opextent[1]>1) ? opextent[1]+2*ng : 1;
  const int nzm = (opextent[2]>1) ? opextent[2]+2*ng : 1;

  int destRowIndex = 0;

  for( int irow=0; irow<op.nrows(); ++irow ){

    // if we are at a ghost entry, skip it.
    const int i = irow%nxm;
    if( i<ng || i>=(nxm-ng)  ) continue;

    const int j = (irow/nxm)%nym;
    if( ny>1 && (j<ng || j>=(nym-ng)) ) continue;

    const int k = irow/(nxm*nym);
    if( nz>1 && (k<ng || k>=(nzm-ng)) ) continue;


    rowDWork_.clear();
    rowIWork_.clear();

    int ncol=0;   double*vals=0;   int*ixs=0;
    op.epetra_mat().ExtractMyRowView( irow, ncol, vals, ixs );

    for( int icol=0; icol<ncol; ++icol ){

      // determine the ijk indices for this column.
      // if we are at a ghost entry, skip it.
      const int colindex = ixs[icol];

      const int i = colindex%nxm;
      if( i<ng || i>=(nxm-ng)  ) continue;

      const int j = (colindex/nxm)%nym;
      if( ny>1 && (j<ng || j>=(nym-ng)) ) continue;

      const int k = colindex/(nxm*nym);
      if( nz>1 && (k<ng || k>=(nzm-ng)) ) continue;

      // insert this value
      rowDWork_.push_back( scaleFac*vals[icol] );

      // now determine the column index for insertion of this value
      const int ii = nxm>1 ? i-ng : 0;
      const int jj = nym>1 ? j-ng : 0;
      const int kk = nzm>1 ? k-ng : 0;
      int iflat = kk*(nx*ny) + jj*(nx) + ii;
      rowIWork_.push_back( iflat );

    } // column loop

    // insert this information into the matrix
    if( rowDWork_.size() > 0 ){
      A_.SumIntoMyValues( destRowIndex, rowDWork_.size(), &rowDWork_[0], &rowIWork_[0] );
    }
    ++destRowIndex;
 
  } // row loop
}
//--------------------------------------------------------------------
void
LHS::add_contribution( const SpatialOps::SpatialField & f,
		       const double scaleFac )
{
  assert( compatibility_check(f) );

  // add non-ghost elements of the local field to this LHS operator

  const int ng  = f.nghost();
  const vector<int> & fextent = f.get_extent();
  const int nxm = (fextent[0]>1) ? fextent[0]+2*ng : 1;
  const int nym = (fextent[1]>1) ? fextent[1]+2*ng : 1;
  const int nzm = (fextent[2]>1) ? fextent[2]+2*ng : 1;

  const double * const fvals = f.get_ptr();

  for( int irow=0; irow<f.get_ntotal(); ++irow ){

    // determine the ijk indices
    const int i = irow%nxm;
    const int j = (irow/nxm)%nym;
    const int k = irow/(nxm*nym);

    // are we at a ghost entry?  If so, go to the next entry.
    if(                  i<ng || i>=nxm-ng  ) continue;
    if( extent_[1]>1 && (j<ng || j>=nym-ng) ) continue;
    if( extent_[2]>1 && (k<ng || k>=nzm-ng) ) continue;

    // add this value to the diagonal
    double val = scaleFac * fvals[irow];
    A_.SumIntoMyValues( irow, 1, &val, &irow );
  }
}
//--------------------------------------------------------------------
void
LHS::Print( std::ostream & c ) const
{
  A_.Print(c);
}
//--------------------------------------------------------------------
bool
LHS::compatibility_check( const SpatialOps::SpatialField& f ) const
{
  return ( extent_ == f.get_extent() );  
}
//--------------------------------------------------------------------
bool
LHS::compatibility_check( const SpatialOps::SpatialOperator& op ) const
{
  return ( extent_ == op.get_extent() );
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
#ifdef HAVE_MPI
LinearSystem::LinearSystem( const vector<int> & extent,
			    MPI_Comm & comm )
#else
LinearSystem::LinearSystem( const vector<int> & extent )
#endif
  : extent_( extent ),
    npts_( std::accumulate(extent.begin(),extent.end(),1,std::multiplies<int>() ) ),
    rhs_( extent ),
    lhs_( NULL ),
    solnFieldValues_( npts_ ),

    maxIterations_  ( 100    ),
    solverTolerance_( 1.0e-12 ),

    A_(NULL), b_(NULL), x_(NULL),

    comm_( NULL ),

    linProb_( NULL ),
    aztec_  ( NULL )
{
  // Build the Epetra Maps
#ifdef HAVE_MPI
  comm_ = &comm;
#else
  comm_ = new Epetra_SerialComm();
#endif
  const Epetra_Map & emap = LinSysMapFactory::self().get_map( npts_, *comm_ );

  // for now, hard code this for 2nd order.
  int entriesPerRow = 3;
  if( extent[1]>1 ) entriesPerRow += 2;
  if( extent[2]>1 ) entriesPerRow += 2;

  // Build the matrix & imprint it with the proper sparsity pattern (for now, hard-coded)
  // Note that if we want periodic BCs, we will need to modify the imprinting.
  A_ = new Epetra_CrsMatrix( Copy, emap, entriesPerRow, true );
  imprint( extent, entriesPerRow );

  lhs_ = new LHS( extent_, *A_ );

  // Build the RHS and LHS vectors - we manage storage (not trilinos)
  x_ = new Epetra_Vector( View, emap, &solnFieldValues_[0] );
  b_ = new Epetra_Vector( View, emap, rhs_.get_ptr() );

  // Build the Linear Problem
  linProb_ = new Epetra_LinearProblem( A_, x_, b_ );

  A_->FillComplete();
  A_->OptimizeStorage();

  // Build the aztec solver
  //aztec_ = new AztecOO( *linProb_ );
}
//--------------------------------------------------------------------
LinearSystem::~LinearSystem()
{
  delete aztec_;
  delete linProb_;
  delete b_;
  delete x_;
  delete A_;
  delete lhs_;

#ifndef HAVE_MPI
  delete comm_;
#endif
}
//--------------------------------------------------------------------
void
LinearSystem::imprint( const std::vector<int> & extent,
		       const int entriesPerRow )
{
  // set the sparsity pattern

  const int nx = extent[0];
  const int ny = extent[1];
  const int nz = extent[2];

  const bool doX = nx > 1;
  const bool doY = ny > 1;
  const bool doZ = nz > 1;

  assert( doX );

  std::vector<double> vals;
  std::vector<int> ixs;

  for( int irow=0; irow<A_->NumMyRows(); ++irow ){

    vals.clear();
    ixs.clear();

    const int i = irow%nx;
    const int j = (irow/ny)%ny;
    const int k = irow/(nx*ny);

    // NOTE: STENCIL IS HARD-CODED!

    if( doX ){
      if( i>0 ){
	const int l = i-1;
	const int icol = k*(nx*ny) + j*nx + l;
	ixs.push_back( icol );
	vals.push_back(1.0);
      }
      ixs.push_back( irow );
      vals.push_back(1.0);
      if( i<nx-1 ){
	const int l = i+1;
	const int icol = k*(nx*ny) + j*nx + l;
	ixs.push_back( icol );
	vals.push_back(1.0);
      }
    }

    if( doY ){
      if( j>0 ){
	const int l = j-1;
	const int icol = k*(nx*ny) + l*nx + i;
	ixs.push_back( icol );
	vals.push_back(1.0);
      }
      if( j<ny-1 ){
	const int l = j+1;
	const int icol = k*(nx*ny) + l*nx + i;
	ixs.push_back( icol );
	vals.push_back(1.0);
      }
    }

    if( doZ ){
      if( k>0 ){
	const int l = k-1;
	const int icol = l*(nx*ny) + j*nx + i;
	ixs.push_back( icol );
	vals.push_back(1.0);
      }
      if( k<nz-1 ){
	const int l = k+1;
	const int icol = l*(nx*ny) + j*nx + i;
	ixs.push_back( icol );
	vals.push_back(1.0);
      }
    }

    A_->InsertGlobalValues( irow, ixs.size(), &vals[0], &ixs[0] );

  } // row loop

}
//--------------------------------------------------------------------
void
LinearSystem::reset()
{
  rhs_.reset();
  lhs_->reset();
}
//--------------------------------------------------------------------
void
LinearSystem::solve()
{
  if( NULL == aztec_ )  aztec_ = new AztecOO( *linProb_ );
  aztec_->Iterate( maxIterations_, solverTolerance_ );
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
LinSysMapFactory::LinSysMapFactory()
{
}
//--------------------------------------------------------------------
LinSysMapFactory::~LinSysMapFactory()
{
  for( std::map<int,Epetra_Map*>::iterator ii=emap_.begin(); ii!=emap_.end(); ++ii )
    delete ii->second;
}
//--------------------------------------------------------------------
LinSysMapFactory&
LinSysMapFactory::self()
{
  static LinSysMapFactory s;
  return s;
}
//--------------------------------------------------------------------
Epetra_Map&
LinSysMapFactory::get_map( const int npts,
			   Epetra_Comm & com )
{
  std::map<int,Epetra_Map*>::const_iterator ii = emap_.find( npts );
  if( ii == emap_.end() ){
    Epetra_Map * myMap = new Epetra_Map( npts, 0, com );
    emap_.insert( std::make_pair(npts,myMap) );
    return *myMap;
  }
  return *(ii->second);
}
//--------------------------------------------------------------------


//====================================================================
