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
RHS::add_contribution( const SpatialOps::SpatialField & f )
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
	field_[ix] += fptr[ixf];
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
LHS::add_contribution( const SpatialOps::SpatialOperator & op )
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

  for( int irow=0; irow<op.nrows(); ++irow ){

    rowDWork_.clear();
    rowIWork_.clear();

    int ncol=0;   double*vals=0;   int*ixs=0;
    op.epetra_mat().ExtractMyRowView( irow, ncol, vals, ixs );

    for( int icol=0; icol<ncol; ++icol ){

      // determine the ijk indices for this column.
      const int colindex = ixs[icol];
      const int i = colindex%nxm;
      const int j = (colindex/nxm)%nym;
      const int k = colindex/(nxm*nym);

      // are we at a ghost entry?  If so, go to the next column entry.
      if(          i<ng || i>=(nxm-ng)  ) continue;
      if( ny>1 && (j<ng || j>=(nym-ng)) ) continue;
      if( nz>1 && (k<ng || k>=(nzm-ng)) ) continue;

      // insert this value
      rowDWork_.push_back( vals[icol] );

      // now determine the column index for insertion of this value
      const int ii = i-ng;
      const int jj = j-ng;
      const int kk = k-ng;
      int iflat = kk*(nx*ny) + jj*(nx) + ii;
      rowIWork_.push_back( iflat );

    } // column loop

    // insert this information into the matrix
    if( rowDWork_.size() > 0 ){
      A_.SumIntoMyValues( irow, rowDWork_.size(), &rowDWork_[0], &rowIWork_[0] );
    }
 
  } // row loop
}
//--------------------------------------------------------------------
void
LHS::add_contribution( const SpatialOps::SpatialField & f )
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
    double val = fvals[irow];
    A_.SumIntoMyValues( irow, 1, &val, &irow );
  }

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
    solnFieldValues_( new double[npts_] )
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

  // Build the matrix
  A_ = new Epetra_CrsMatrix( Copy, emap, entriesPerRow, true );
  lhs_ = new LHS( extent_, *A_ );

  // Build the RHS and LHS vectors - we manage storage (not trilinos)
  x_ = new Epetra_Vector( View, emap, solnFieldValues_ );
  b_ = new Epetra_Vector( View, emap, rhs_.get_ptr() );

  // Build the Linear Problem
  linProb_ = new Epetra_LinearProblem( A_, x_, b_ );

  // Build the aztec solver
  aztec_ = new AztecOO( *linProb_ );
}
//--------------------------------------------------------------------
LinearSystem::~LinearSystem()
{
  delete aztec_;
  delete linProb_;
  delete b_;
  delete x_;
  delete A_;

  delete [] solnFieldValues_;

#ifndef HAVE_MPI
  delete comm_;
#endif
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
