#include <numeric>
#include <sstream>
#include <stdexcept>

#include <spatialops/LinearSystem.h>
#include <spatialops/SpatialOperator.h>
#include <spatialops/LinAlgTrilinos.h>

#include <Epetra_Map.h>
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

/**
 *  @class  LinSysMapFactory
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Support for creating Epetra_Map objects for use with the
 *  LinearSystem class.
 */
class LinSysMapFactory
{
public:

  static LinSysMapFactory& self()
  {
    static LinSysMapFactory s;
    return s;
  }

  Epetra_Map& get_map( const int npts,
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

private:

  std::map<int,Epetra_Map*> emap_;

  LinSysMapFactory()
  {
  }

  ~LinSysMapFactory()
  {
    for( std::map<int,Epetra_Map*>::iterator ii=emap_.begin(); ii!=emap_.end(); ++ii )
      delete ii->second;
  }

};

//====================================================================


//--------------------------------------------------------------------
RHS::RHS( const int npts )
  : npts_( npts )
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
RHS::reset( const int rownum, const double val )
{
  field_[rownum] = val;
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
LHS::LHS( const SpatialOps::structured::IntVec& extent,
          const bool hasPlusXSideFaces,
          const bool hasPlusYSideFaces,
          const bool hasPlusZSideFaces,
          Epetra_CrsMatrix& A )
  : A_( A ),
    extent_( extent ),
    hasPlusXSideFaces_( hasPlusXSideFaces ),
    hasPlusYSideFaces_( hasPlusYSideFaces ),
    hasPlusZSideFaces_( hasPlusZSideFaces ),
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
LHS::unit_diagonal_zero_else( const int irow )
{
  int rownum = irow; // trilinos wants a non-const int.  Stupid.
  int n;
  double * vals;
  int * ixs;
  A_.ExtractMyRowView( rownum, n, vals, ixs );

  for( int i=0; i<n; ++i ){
    vals[i] = 0.0;
    if( ixs[i] == irow ) vals[i] = 1.0;
  }
}
//--------------------------------------------------------------------
void
LHS::add_diag_contribution( double x )
{
  for( int irow=0; irow<nrows_; ++irow ){
    A_.SumIntoMyValues( irow, 1, &x, &irow );
  }
}
//--------------------------------------------------------------------
void
LHS::Print( std::ostream & c ) const
{
  A_.Print(c);
}
//--------------------------------------------------------------------
void
LHS::write_matlab( const std::string prefix ) const
{
  const std::string fname = "load_"+prefix +".m";
  std::ofstream fout( fname.c_str() );
  fout << "function A = load_" << prefix << "()" << std::endl;

  // first time through count nonzeros for preallocation in matlab.
  // second time through, write the entries.
  for( int writeLoop=0; writeLoop<=1; ++writeLoop ){
    int i=0;
    for( int irow=0; irow<nrows_; ++irow ){
      SpatialOps::LinAlgTrilinos::MatrixRow row( irow, A_ );
      for( SpatialOps::LinAlgTrilinos::column_iterator icol = row.begin();
           icol!=row.end();
           ++icol, ++i )
        {
          if( writeLoop==1 ){
            fout << "row(" << i+1 << ") = " << irow+1 << ";  "
                 << "col(" << i+1 << ") = " << icol.index()+1 << ";  "
                 << "val(" << i+1 << ") = " << *icol << ";"
                 << std::endl;
          }
        }
    }
    if( writeLoop==0 ){
      fout << "row = zeros(" << i << ",1);  col=row;  val=row;" << std::endl;
    }
  }
  fout << "A = sparse( "
       << " row, col, val, "
       << nrows_ << ", " << ncols_
       << ");" << std::endl;
  fout.close();
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
#ifdef HAVE_MPI
LinearSystem::LinearSystem( const SpatialOps::structured::IntVec & extent,
                            const bool hasPlusXSideFaces,
                            const bool hasPlusYSideFaces,
                            const bool hasPlusZSideFaces,
                            MPI_Comm & comm )
  : npts_( get_global_npts(extent,comm) ),
#else
LinearSystem::LinearSystem( const SpatialOps::structured::IntVec & extent,
                            const bool hasPlusXSideFaces,
                            const bool hasPlusYSideFaces,
                            const bool hasPlusZSideFaces )
  : npts_( get_global_npts(extent) ),
#endif

    extent_( extent ),

    rhs_( npts_ ),  // construct the rhs
    lhs_( NULL  ),
    solnFieldValues_( npts_ ),  // construct the solution

    maxIterations_  ( 100    ),
    solverTolerance_( 1.0e-12 ),

    A_(NULL), b_(NULL), x_(NULL),

#ifdef HAVE_MPI
    comm_( new Epetra_MpiComm(comm) ),
#else
    comm_( new Epetra_SerialComm() ),
#endif

    linProb_( NULL ),
    aztec_  ( NULL )
{
  // Build the Epetra Maps
  const Epetra_Map & emap = LinSysMapFactory::self().get_map( npts_, *comm_ );

  // for now, hard code this for 2nd order.
  int entriesPerRow = 1;
  if( extent[0]>1 ) entriesPerRow += 2;
  if( extent[1]>1 ) entriesPerRow += 2;
  if( extent[2]>1 ) entriesPerRow += 2;

  // Build the matrix & imprint it with the proper sparsity pattern (for now, hard-coded)
  // Note that if we want periodic BCs, we will need to modify the imprinting.
  A_ = new Epetra_CrsMatrix( Copy, emap, entriesPerRow, true );
  imprint( extent, entriesPerRow );

  lhs_ = new LHS( extent_, hasPlusXSideFaces, hasPlusYSideFaces, hasPlusZSideFaces, *A_ );

  // Build the RHS and LHS vectors - we manage storage (not trilinos)
  x_ = new Epetra_Vector( View, emap, solnFieldValues_.begin() );
  b_ = new Epetra_Vector( View, emap, rhs_.begin() );

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
#ifdef HAVE_MPI
int
LinearSystem::get_global_npts( const SpatialOps::structured::IntVec & extent,
                               MPI_Comm & comm )
{
  // determine the total number of points in each dimension.
  SpatialOps::structured::IntVec tmpPts = extent;
  SpatialOps::structured::IntVec globPoints( extent.size(), 0 );
  MPI_Allreduce( &tmpPts[0], &globPoints[0], extent.size(), MPI_INT, MPI_SUM, comm );

  return std::accumulate( globPoints.begin(), globPoints.end(), 1, std::multiplies<int>() );
}
#else
int
LinearSystem::get_global_npts( const SpatialOps::structured::IntVec & extent )
{
  return extent[0]*extent[1]*extent[2];
}
#endif
//--------------------------------------------------------------------
void
LinearSystem::imprint( const SpatialOps::structured::IntVec& extent,
                       const int entriesPerRow )
{
  // set the sparsity pattern

  const int nx = extent[0];
  const int ny = extent[1];
  const int nz = extent[2];

  const bool doX = nx > 1;
  const bool doY = ny > 1;
  const bool doZ = nz > 1;

  std::vector<double> vals;
  std::vector<int> ixs;

  for( int irow=0; irow<A_->NumMyRows(); ++irow ){

    vals.clear();
    ixs.clear();

    const int i = irow%nx;
    const int j = (irow/nx)%ny;
    const int k = irow/(nx*ny);

    // NOTE: STENCIL IS HARD-CODED!

    ixs.push_back( irow );
    vals.push_back(1.0);

    if( doX ){
      if( i>0 ){
        const int l = i-1;
        const int icol = k*(nx*ny) + j*nx + l;
        ixs.push_back( icol );
        vals.push_back(1.0);
      }
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
  if( NULL == aztec_ ){
    aztec_ = new AztecOO( *linProb_ );
    aztec_->SetAztecOption(AZ_output, AZ_none  );
    aztec_->SetAztecOption(AZ_precond,AZ_Jacobi);
    aztec_->SetAztecOption(AZ_solver, AZ_gmres ); 
  }
  aztec_->Iterate( maxIterations_, solverTolerance_ );
}
//--------------------------------------------------------------------
void
LinearSystem::set_dirichlet_condition( const int irow,
                                       const double rhsVal )
{
  if( irow>npts_ ){
    std::ostringstream msg;
    msg << "ERROR: LinearSystem::set_dirichlet_condition()" << std::endl
        << "       Invalid position to set dirichlet BC on (row=" << irow << ")." << std::endl
        << "       There are only " << npts_ << " rows in the matrix." << std::endl;
    throw std::runtime_error( msg.str() );
  }

  // set the LHS row.
  lhs_->unit_diagonal_zero_else( irow );

  // set the rhs value
  rhs_.reset( irow, rhsVal );
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
#ifdef HAVE_MPI
LinSysInfo::LinSysInfo( const SpatialOps::structured::IntVec& npts,
                        const bool hasPlusXSideFaces,
                        const bool hasPlusYSideFaces,
                        const bool hasPlusZSideFaces,
                        MPI_Comm & communicator )
  : comm( communicator ),
#else
LinSysInfo::LinSysInfo( const SpatialOps::structured::IntVec& npts,
                        const bool hasPlusXSideFaces,
                        const bool hasPlusYSideFaces,
                        const bool hasPlusZSideFaces )
  :
#endif
    dimExtent( npts ),
    hasPlusXFaces( hasPlusXSideFaces ),
    hasPlusYFaces( hasPlusYSideFaces ),
    hasPlusZFaces( hasPlusZSideFaces )
{
  solverPackage  = TRILINOS;
  preconditioner = DEFAULT;
}
//--------------------------------------------------------------------
LinSysInfo::~LinSysInfo()
{
}
//--------------------------------------------------------------------
bool
LinSysInfo::operator==(const LinSysInfo& s) const
{
  return( s.solverPackage  == solverPackage &&
          s.preconditioner == preconditioner &&
          s.dimExtent == dimExtent
#ifdef HAVE_MPI
          && s.comm == comm
#endif
    );
}
//--------------------------------------------------------------------
bool
LinSysInfo::operator<(const LinSysInfo& s) const
{
  return( s.solverPackage < solverPackage &&
          s.preconditioner< preconditioner &&
          s.dimExtent < dimExtent
#ifdef HAVE_MPI
          && s.comm < comm
#endif
    );
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
LinSysFactory::LinSysFactory()
{
}
//--------------------------------------------------------------------
LinSysFactory::~LinSysFactory()
{
  for( InfoMap::iterator ii=infoMap_.begin();
       ii!=infoMap_.end();
       ++ii )
    {
      delete ii->second;
    }
}
//--------------------------------------------------------------------
LinSysFactory&
LinSysFactory::self()
{
  static LinSysFactory s;
  return s;
}
//--------------------------------------------------------------------
LinearSystem&
LinSysFactory::get_linsys( const LinSysInfo& info )
{
  // do we have a match?
  InfoMap::iterator ii = infoMap_.find( info );

  if( ii==infoMap_.end() ){

    // build a new linear system.

#ifdef HAVE_MPI
    LinearSystem * linSys = new LinearSystem( info.dimExtent,
                                              info.hasPlusXFaces,
                                              info.hasPlusYFaces,
                                              info.hasPlusZFaces,
                                              info.comm );
#else
    LinearSystem * linSys = new LinearSystem( info.dimExtent,
                                              info.hasPlusXFaces,
                                              info.hasPlusYFaces,
                                              info.hasPlusZFaces );
#endif

    std::pair<InfoMap::iterator,bool> result =
      infoMap_.insert( std::make_pair(info,linSys) );

    assert( result.second );

    ii = result.first;
  }
  return *(ii->second);
}
//--------------------------------------------------------------------
void
LinSysFactory::bind_name_to_info( const LinSysInfo& info,
                                  const std::string & name )
{
  nameInfoMap_.insert( make_pair(name,info) );
}
//--------------------------------------------------------------------
LinearSystem&
LinSysFactory::get_linsys( const std::string & name )
{
  NameInfoMap::iterator ii=nameInfoMap_.find( name );
  if( ii==nameInfoMap_.end() ){
    std::ostringstream msg;
    msg << "ERROR: No linear system has been assigned to '" << name << "'." << std::endl
        << "       You must first bind the name to the linear system information." << std::endl;
    throw std::runtime_error( msg.str() );
  }
  return get_linsys( ii->second );
}

//====================================================================
