#ifndef UT_LinearSystem_h
#define UT_LinearSystem_h

#include <spatialops/SpatialOpsConfigure.h>

#include <vector>
#include <map>
#include <string>

#include <spatialops/FVToolsTemplates.h>
#include <spatialops/SpatialField.h>

//---------------------------------
// trilinos includes
#ifndef LINALG_TRILINOS
# error LinearSystem support requires Trilinos
#endif // LINALG_TRILINOS

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


/** @todo rip out the extent information here. */

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
              const bool hasPlusXSideFaces,
              const bool hasPlusYSideFaces,
              const bool hasPlusZSideFaces,
              MPI_Comm & communicator );
#else
  LinSysInfo( const std::vector<int> & npts,
              const bool hasPlusXSideFaces,
              const bool hasPlusYSideFaces,
              const bool hasPlusZSideFaces );
#endif

  ~LinSysInfo();

  SolverPackage solverPackage;
  Preconditioner preconditioner;

  //
  // these are required to build the linear system.
  //

#ifdef HAVE_MPI
  MPI_Comm & comm;
#endif

  const std::vector<int> dimExtent;
  const bool hasPlusXFaces, hasPlusYFaces, hasPlusZFaces;

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
 *
 *  @todo Implement unary operator version of add_field_contribution
 */
class RHS
{
public:

  RHS( const int npts );

  ~RHS();

  void reset( const double val = 0.0 );

  void reset( const int rownum, const double value=0.0 );

  template< typename FieldType >
  void add_field_contribution( const FieldType & localField,
                               const double scaleFac = 1.0 );

  inline const std::vector<double> & get_field() const{return field_;}

  inline double& operator[](const int i){return field_[i];}

  inline RHS& operator=(const double val){ reset(val); return *this; }


  /** @name unary operators */
  //@{
  inline RHS& operator+=(const double val){ for(iterator i=begin(); i!=end(); ++i) *i+=val;  return *this; }
  inline RHS& operator-=(const double val){ for(iterator i=begin(); i!=end(); ++i) *i-=val;  return *this; }
  inline RHS& operator*=(const double val){ for(iterator i=begin(); i!=end(); ++i) *i*=val;  return *this; }
  inline RHS& operator/=(const double val){ return(*this *= (1.0/val)); }
  //@}


  typedef double*             iterator;
  typedef double const* const_iterator;

  inline iterator begin(){ return &field_[0]; }
  inline iterator end()  { return &field_[npts_]; }

  inline const_iterator begin() const{ return &field_[0]; }
  inline const_iterator end()   const{ return &field_[npts_]; }


  typedef iterator interior_iterator;
  typedef const_iterator const_interior_iterator;
  inline const_interior_iterator interior_begin() const{return begin();}
  inline       interior_iterator interior_begin()      {return begin();}
  inline const_interior_iterator interior_end()   const{return end();}
  inline       interior_iterator interior_end()        {return end();}

private:

  template<typename FieldType>
  bool consistency_check( const FieldType& f ) const;

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
       const bool hasPlusXSideFaces,
       const bool hasPlusYSideFaces,
       const bool hasPlusZSideFaces,
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

  /**
   *  Add non-ghost elements of the local matrix to this LHS operator.
   *
   *  @todo This WILL NOT WORK in PARALLEL because we are eliminating
   *  important entries!  Need to figure out how to make this general.
   *  This will require knowing if we are on a processor boundary or a
   *  domain boundary for each ghost point.  Ouch.
   */
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

  void write_matlab( const std::string prefix ) const;

private:

  Epetra_CrsMatrix & A_;
  const std::vector<int> extent_;
  const bool hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_;
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
 *  @todo Implement ghost support for global distributed system?  This
 *  means that the RHS field would have ghosting as well.  Somewhat
 *  strange, but this is how we will do BCs...
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
  const Epetra_Vector& get_rhs_field_epetra_vec()  const{ return *b_; }

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
                const bool hasPlusXSideFaces,
                const bool hasPlusYSideFaces,
                const bool hasPlusZSideFaces,
                MPI_Comm & comm );

  int get_global_npts( const std::vector<int> & extent,
                       MPI_Comm & comm );

#else

  LinearSystem( const std::vector<int> & dimExtent,
                const bool hasPlusXSideFaces,
                const bool hasPlusYSideFaces,
                const bool hasPlusZSideFaces );

  int get_global_npts( const std::vector<int> & extent );

#endif

  ~LinearSystem();


  void imprint( const std::vector<int> &, const int );

  const int npts_;
  const std::vector<int> extent_;
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

//--------------------------------------------------------------------
template<typename FieldType>
void
RHS::add_field_contribution( const FieldType& f,
                             const double scaleFac )
{
  typename FieldType::const_interior_iterator ifld = f.interior_begin();
  const typename FieldType::const_interior_iterator iflde = f.interior_end();

  iterator irhs = begin();
#ifndef NDEBUG
  const iterator irhse = end();
#endif

  if( scaleFac==1 ){
    for( ; ifld!=iflde; ++ifld, ++irhs ){
      assert( irhs != irhse );
      *irhs += *ifld;
    }
  }
  else if( scaleFac==-1 )
    for( ; ifld!=iflde; ++ifld, ++irhs ){
      *irhs -= *ifld;
    }
  else{
    for( ; ifld!=iflde; ++ifld, ++irhs ){
      assert( irhs != irhse );
      *irhs += *ifld * scaleFac;
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
  using namespace SpatialOps;

  // add non-ghost elements of the local matrix to this LHS operator
  const int nx=extent_[0];
  const int ny=extent_[1];
  const int nz=extent_[2];

  FVStaggered::IndexTriplet t;

  int irow = 0;
  for( int ioprow=0; ioprow<op.nrows(); ++ioprow ){

    if( op.is_row_ghost(ioprow) ) continue;

    rowDWork_.clear();
    rowIWork_.clear();

    const typename OpType::MatrixRow row( op.get_row(ioprow) );
    typename OpType::const_column_iterator icol = row.begin();
    typename OpType::const_column_iterator icole= row.end();

    for( ; icol!=icole; ++icol ){

      // determine the ijk indices for this column.
      // if we are at a ghost entry, skip it.
      const int colindex = icol.index();

      // check to see if we are at a ghost entry - get the "interior" indices.
      if( op.is_col_ghost( colindex ) ) continue;

      // insert this value
      rowDWork_.push_back( scaleFac * *icol );

      // now determine the column index for insertion of this value
      typedef typename OpType::SrcFieldType SrcField;
      t = FVStaggered::flat2ijk<SrcField>::value( extent_, colindex );
      if( nx>1 ) t.i -= SrcField::Ghost::NM;
      if( ny>1 ) t.j -= SrcField::Ghost::NM;
      if( nz>1 ) t.k -= SrcField::Ghost::NM;

      const int iflat = t.k*(nx*ny) + t.j*(nx) + t.i;
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
  using std::vector;

  const typename FieldType::const_interior_iterator iflde = f.interior_end();
  typename FieldType::const_interior_iterator ifld = f.interior_begin();

  int ientry=0;
  for( ; ifld!=iflde; ++ifld, ++ientry ){
    // add this value to the diagonal
    double val = scaleFac * (*ifld);
    A_.SumIntoMyValues( ientry, 1, &val, &ientry );
  }

}
//--------------------------------------------------------------------

//====================================================================








//====================================================================

namespace SpatialOps{

  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator=( const RHS& rhs )
  {
    // assign the RHS field to this spatial field.
    typename RHS::const_iterator irhs = rhs.begin();
#ifndef NDEBUG
    const typename RHS::const_iterator irhse = rhs.end();
#endif
    const interior_iterator iflde = this->interior_end();
    interior_iterator ifld = this->interior_begin();
    
    for( ; ifld!=iflde; ++ifld, ++irhs ){
      assert( irhs!=irhse );
      *ifld = *irhs;
    }
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator+=( const RHS& rhs )
  {
    // assign the RHS field to this spatial field.
    typename RHS::const_iterator irhs = rhs.begin();
#ifndef NDEBUG
    const typename RHS::const_iterator irhse = rhs.end();
#endif
    const interior_iterator iflde = this->interior_end();
    interior_iterator ifld = this->interior_begin();
    
    for( ; ifld!=iflde; ++ifld, ++irhs ){
      assert( irhs!=irhse );
      *ifld += *irhs;
    }
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator-=( const RHS& rhs )
  {
    // assign the RHS field to this spatial field.
    typename RHS::const_iterator irhs = rhs.begin();
    const typename RHS::const_iterator irhse = rhs.end();
    
    const interior_iterator iflde = this->interior_end();
    interior_iterator ifld = this->interior_begin();
    
    for( ; ifld!=iflde; ++ifld, ++irhs ){
      assert( irhs!=irhse );
      *ifld -= *irhs;
    }
    return *this;
  }
}// namespace SpatialOps

#endif
