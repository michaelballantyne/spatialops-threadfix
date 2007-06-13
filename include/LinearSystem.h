#ifndef UT_LinearSystem_h
#define UT_LinearSystem_h

#include <vector>
#include <map>
#include <string>

#include <SpatialOpsDefs.h>

//---------------------------------
// trilinos includes
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
// trilinos includes
//---------------------------------


// trilinos foward declarations
class Epetra_Map;
class Epetra_LinearProblem;
class AztecOO;


// local forward declarations
class LinearSystem;


//====================================================================

/**
 *  @struct LinSysInfo
 *  @author James C. Sutherland
 *  @date   January, 2007
 *
 *  Information required to obtain a LinearSystem from the
 *  LinSysFactory.
 */
struct LinSysInfo
{

  enum SolverPackage{
    TRILINOS  // currently only TRILINOS is supported.
  };

  enum Preconditioner{
    NONE,
    DEFAULT
  };

#ifdef HAVE_MPI
  LinSysInfo( const std::vector<int> & npts,
	      MPI_Comm & communicator );
#else
  LinSysInfo( const std::vector<int> & npts );
#endif

  ~LinSysInfo();

  SolverPackage solverPackage;
  Preconditioner preconditioner;

  //
  // these are required to build the linear system.
  //

  const std::vector<int> dimExtent;

#ifdef HAVE_MPI
  MPI_Comm & comm;
#endif

  bool operator ==(const LinSysInfo&) const;
  bool operator < (const LinSysInfo&) const;
};


//====================================================================


/**
 *  @class  LinSysFactory
 *  @author James C. Sutherland
 *  @date   January, 2007
 *
 *  Factory to produce linear systems.
 */
class LinSysFactory
{
public:
  static LinSysFactory& self();

  void bind_name_to_info( const LinSysInfo& info,
			  const std::string & name );

  LinearSystem & get_linsys( const LinSysInfo& info );

  LinearSystem & get_linsys( const std::string & name );

private:
  LinSysFactory();
  ~LinSysFactory();

  typedef std::map<LinSysInfo,LinearSystem*> InfoMap;
  InfoMap infoMap_;

  typedef std::map<std::string,LinSysInfo> NameInfoMap;
  NameInfoMap nameInfoMap_;
};


//====================================================================


/**
 *  @class  RHS
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Abstraction for the RHS field.  This is intended for use in a
 *  segregated implicit scheme where the RHS is used in conjunction
 *  with a LHS operator, or in an explicit scheme where the RHS is
 *  simply used to obtain an update for the solution field directly.
 */
class RHS
{
public:

  RHS( const std::vector<int> & domainExtent );

  ~RHS();

  void reset( const double val = 0.0 );

  void reset( const int rownum, const double value=0.0 );

  template< typename FieldType >
  void add_field_contribution( const FieldType & localField,
			       const double scaleFac = 1.0 );

  inline const std::vector<double> & get_field() const{return field_;}

  inline const std::vector<int> & get_extent() const{return extent_;}


  inline RHS& operator=(const double val){ reset(val); return *this; }


  typedef double* iterator;
  typedef double const* const_iterator;

  inline iterator begin(){ return &field_[0]; }
  inline iterator end()  { return &field_[npts_]; }

  inline const_iterator begin() const{ return &field_[0]; }
  inline const_iterator end()   const{ return &field_[npts_]; }

private:

  template<typename FieldType>
  bool consistency_check( const FieldType& f ) const;

  const std::vector<int> extent_;
  const int npts_;

  std::vector<double> field_;

  RHS();
  RHS(const RHS&);
};

//====================================================================

/**
 *  @class  LHS
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Abstraction for the LHS operator.  This is a distributed LHS
 *  operator for use in a linear system.
 */
class LHS
{
public:

  LHS( const std::vector<int> & extent,
       Epetra_CrsMatrix& A );

  ~LHS();

  int nrows() const{ return nrows_; }
  int ncols() const{ return ncols_; }

  /** reset all values in the matrix */
  void reset( const double val = 0 );

  /**
   *  For setting dirichlet conditions.  Places unity on the diagonal
   *  entry of this row and zeros elsewhere.
   */
  void unit_diagonal_zero_else( const int irow );

  /** add non-ghost elements of the local matrix to this LHS operator */
  template<typename OpType>
  void add_op_contribution( const OpType & localMat,
			    const double scaleFac = 1.0 );

  /** add non-ghost elements of the local field to the diagonal of LHS operator */
  template<typename FieldType>
  void add_field_contribution( const FieldType & localField,
			       const double scaleFac = 1.0 );


  /** add a constant to the diagonal of the matrix */
  void add_diag_contribution( double x );

        Epetra_CrsMatrix& epetra_mat()      { return A_; }
  const Epetra_CrsMatrix& epetra_mat() const{ return A_; }

  void Print( std::ostream& c ) const;

private:

  template<typename FieldType>
  bool field_compatibility_check( const FieldType& t ) const;

  template<typename OpType>
  bool op_compatibility_check( const OpType& op ) const;

  Epetra_CrsMatrix & A_;
  const std::vector<int> extent_;
  const int nrows_;
  const int ncols_;

  LHS(const LHS& );
  LHS();

  std::vector<double> rowDWork_;
  std::vector<int>    rowIWork_;
};

//====================================================================

typedef RHS SOLN;

/**
 *  @class  LinearSystem
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Basic support for a linear system distributed in parallel.
 *
 *  LinearSystem objects are constructed via the LinSysFactory.
 */
class LinearSystem
{
  friend class LinSysFactory;

public:

  /** Zero the lhs and rhs */
  void reset();

  void solve();

        RHS & get_rhs()      { return rhs_; }
  const RHS & get_rhs() const{ return rhs_; }

        LHS & get_lhs()      { return *lhs_; }
  const LHS & get_lhs() const{ return *lhs_; }

  const SOLN & get_soln_field() const{ return solnFieldValues_; }
        SOLN & get_soln_field()      { return solnFieldValues_; }

  const Epetra_Vector& get_soln_field_epetra_vec() const{ return *x_; }

  void set_tolerance( const double tol ){ solverTolerance_=tol; }

  void set_maxiter( const int maxit ){ maxIterations_=maxit; }



  /**
   *  Dirichlet conditions are accomplished by zeroing all columns of
   *  the given row and placing a "1" on the diagonal. The rhs vector
   *  is also set to accomplish the desired solution value.
   *
   *  In the case where this linear system is being applied to a
   *  residual update, the lhs value should be zero, which is set as a
   *  default.
   */
  void set_dirichlet_condition( const int irow,
				const double rhsVal = 0 );


protected:


#ifdef HAVE_MPI

  LinearSystem( const std::vector<int> & dimExtent,
		MPI_Comm & comm );

  int LinearSystem::get_global_npts( const std::vector<int> & extent,
				     MPI_Comm & comm );

#else

  LinearSystem( const std::vector<int> & dimExtent );

  int LinearSystem::get_global_npts( const std::vector<int> & extent );

#endif

  ~LinearSystem();


  void imprint( const std::vector<int> &, const int );

  const std::vector<int> extent_;
  const int npts_;
  RHS   rhs_;
  LHS * lhs_;
  SOLN  solnFieldValues_;

  int maxIterations_;
  double solverTolerance_;

private:

  Epetra_CrsMatrix * A_;  // the distributed matrix
  Epetra_Vector    * b_;  // the distributed RHS vector
  Epetra_Vector    * x_;  // the distributed solution vector

  Epetra_Comm      * comm_;

  Epetra_LinearProblem * linProb_;
  AztecOO * aztec_;

};

//====================================================================






//====================================================================

template<typename FieldType>
bool
RHS::consistency_check( const FieldType& f ) const
{
  bool ok = f.get_extent() == extent_;
  return ok;
}
//--------------------------------------------------------------------
template<typename FieldType>
void
RHS::add_field_contribution( const FieldType& f,
			     const double scaleFac )
{
  assert( consistency_check(f) );

  using namespace SpatialOps;

  /*
   *  Procedure:
   *    1. Determine how many ghost cells are in the field
   *    2. Sum field elements into RHS, skipping ghost values
   */

  const int nx=extent_[0];
  const int ny=extent_[1];
  const int nz=extent_[2];

  const int ngxl = f.template nghost<XDIR,SideMinus>();
  const int ngxr = f.template nghost<XDIR,SidePlus >();
  const int ngyl = f.template nghost<YDIR,SideMinus>();
  const int ngyr = f.template nghost<YDIR,SidePlus >();
  const int ngzl = f.template nghost<ZDIR,SideMinus>();
  const int ngzr = f.template nghost<ZDIR,SidePlus >();

  // get the dimensions of the field
  const int nxf= (extent_[0]>1) ? extent_[0] + ngxl + ngxr : 1;
  const int nyf= (extent_[1]>1) ? extent_[1] + ngyl + ngyr : 1;
  const int nzf= (extent_[2]>1) ? extent_[2] + ngzl + ngzr : 1;

  const int yskip = ngxl+ngxr;
  const int zskip = nxf * (ngyl+ngyr);

  int ixf=ngxl;
  if( nyf>1 ) ixf += nxf;
  if( nzf>1 ) ixf += nxf*nyf;

  typename FieldType::const_iterator ifld = f.begin() + ixf;
  std::vector<double>::iterator irhs = field_.begin();

  if( scaleFac==1.0 ){
    for( int k=0; k<nz; ++k ){
      for( int j=0; j<ny; ++j ){
	for( int k=0; k<nx; ++k ){
	  *irhs++ += *ifld++;
	}
	ifld += yskip;
      }
      ifld += zskip;
    }
  }
  else{
    for( int k=0; k<nz; ++k ){
      for( int j=0; j<ny; ++j ){
	for( int k=0; k<nx; ++k ){
	  *irhs++ += scaleFac * (*ifld++);
	}
	ifld += yskip;
      }
      ifld += zskip;
    }
  }
}
//--------------------------------------------------------------------

//====================================================================

//--------------------------------------------------------------------
template<typename OpType>
void
LHS::add_op_contribution( const OpType & op,
			  const double scaleFac )
{
  assert( op_compatibility_check(op) );

  // add non-ghost elements of the local matrix to this LHS operator
  const int nx=extent_[0];
  const int ny=extent_[1];
  const int nz=extent_[2];

  typename OpType::IndexTriplet ix;

  int irow = 0;
  for( int ioprow=0; ioprow<op.nrows(); ++ioprow ){

    if( op.is_row_ghost(ioprow) ) continue;

    rowDWork_.clear();
    rowIWork_.clear();

    int ncol=0;   double*vals=0;   int*ixs=0;
    op.get_linalg_mat().ExtractMyRowView( ioprow, ncol, vals, ixs );

    for( int icol=0; icol<ncol; ++icol ){

      // determine the ijk indices for this column.
      // if we are at a ghost entry, skip it.
      const int colindex = ixs[icol];

      // check to see if we are at a ghost entry - get the "interior" indices.
      if( op.is_col_ghost( colindex, &ix ) ) continue;

      // insert this value
      rowDWork_.push_back( scaleFac*vals[icol] );

      // now determine the column index for insertion of this value
      const int ii = nx>1 ? ix[0] : 0;
      const int jj = ny>1 ? ix[1] : 0;
      const int kk = nz>1 ? ix[2] : 0;
      const int iflat = kk*(nx*ny) + jj*(nx) + ii;
      rowIWork_.push_back( iflat );

    } // column loop

    // insert this information into the matrix
    if( rowDWork_.size() > 0 ){
      const int errCode = A_.SumIntoMyValues( irow, rowDWork_.size(), &rowDWork_[0], &rowIWork_[0] );
      if( errCode != 0 ){
	std::cout << "ERROR: code '" << errCode << "' returned from trilinos SumIntoMyValues." << std::endl
		  << "       This is likely due to an invalid column index for summation." << std::endl;
      }
    }
    ++irow;
 
  } // row loop
}
//--------------------------------------------------------------------
template<typename FieldType>
void
LHS::add_field_contribution( const FieldType & f,
			     const double scaleFac )
{
  //  assert( compatibility_check(f) );

  using namespace SpatialOps;
  using std::vector;

  const int ngxl = f.template nghost<XDIR,SideMinus>();
  const int ngxr = f.template nghost<XDIR,SidePlus >();
  const int ngyl = f.template nghost<YDIR,SideMinus>();
  const int ngyr = f.template nghost<YDIR,SidePlus >();
  const int ngzl = f.template nghost<ZDIR,SideMinus>();
  const int ngzr = f.template nghost<ZDIR,SidePlus >();

  // add non-ghost elements of the local field to this LHS operator

  const vector<int> & fextent = f.get_extent();
  const int nxm = (fextent[0]>1) ? fextent[0] + ngxl + ngxr : 1;
  const int nym = (fextent[1]>1) ? fextent[1] + ngyl + ngyr : 1;
  const int nzm = (fextent[2]>1) ? fextent[2] + ngzl + ngzr : 1;

  int ientry=0;
  int irow=0;
  for( typename FieldType::const_iterator ifld=f.begin(); ifld!=f.end(); ++ifld, ++irow ){

    // determine the ijk indices
    const int i = irow%nxm;
    const int j = (irow/nxm)%nym;
    const int k = irow/(nxm*nym);

    // are we at a ghost entry?  If so, go to the next entry.
    if(                  i<ngxl || i>=nxm-ngxr  ) continue;
    if( extent_[1]>1 && (j<ngyl || j>=nym-ngyr) ) continue;
    if( extent_[2]>1 && (k<ngzl || k>=nzm-ngzr) ) continue;

    // add this value to the diagonal
    double val = scaleFac * (*ifld);
    A_.SumIntoMyValues( ientry, 1, &val, &ientry );
    ++ientry;
  }
}
//--------------------------------------------------------------------
template<typename FieldType>
bool
LHS::field_compatibility_check( const FieldType& f ) const
{
  return ( extent_ == f.get_extent() );  
}
//--------------------------------------------------------------------
template<typename OpType>
bool
LHS::op_compatibility_check( const OpType& op ) const
{
  return ( extent_ == op.get_extent() );
}
//--------------------------------------------------------------------

//====================================================================

#endif
