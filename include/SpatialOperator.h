#ifndef UT_SpatialOperator_h
#define UT_SpatialOperator_h

#include <boost/static_assert.hpp>

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

#include <SpatialOpsDefs.h>
#include <SpatialOpsTools.h>
#include <SpatialField.h>

namespace SpatialOps{


  //====================================================================


  /**
   *  @struct OpAssemblerSelector
   *  @author James C. Sutherland
   *
   *  @brief Defines the type of Assembler required to build a given
   *  SpatialOperator.
   *
   *  This should be specialized for each operator type to provide the
   *  type of assembler required for the operator.
   *
   *  @par Template Parameters
   *    \li \b OpType The type of operator.
   *    \li \b SrcFieldT  Specifies information for the source field.
   *    \li \b DestFieldT Specifies information for the destination field.
   *
   *  Specialized versions of the OpAssemblerSelector struct must
   *  specialize one or more of the template arguments. As this is an
   *  empty struct, it simply defines an interface.  Thus, if a
   *  specialized match is not found, the compiler will fail.
   *
   *  Specialized versions must provide a typedef that defines an \c
   *  Assembler type which is the type of Assembler to be used to
   *  construct the SpatialOperator.  This is required as an input to
   *  the constructor.  As an example, assume that we have an operator
   *  assembler to build a SpatialOperator of type \c MyOp - let's
   *  call it \c MyOpAssembler If MyOp was templated on the direction
   *  and traits of the source and destination fields, e.g.,
   *
   *  \code
   *    template< typename SrcFieldT,
   *              typename DestFieldT >
   *    class MyOpAssembler
   *    {
   *      ...
   *    };
   *  \endcode
   *
   *  Then we could define an OpAssemblerSelector for \c MyOp objects as
   *
   *  \code
   *    template< typename SrcFieldT,
   *              typename DestFieldT >
   *    OpAssemblerSelector< MyOp, SrcFieldT, DestFieldT >
   *    {
   *      typedef MyOpAssembler<SrcFieldT,DestFieldT>  Assembler;
   *    };
   *  \endcode
   *
   *
   *  @par The Assembler
   *
   *  An Assembler must provide the following methods:
   *
   *   \li A method to return the number of rows in the operator.
   *   This should have the following signature: \code int get_nrows()
   *   \endcode
   *
   *   \li A method to return the number of columns in the operator.
   *   This should have the following signature: \code int get_ncols()
   *   \endcode
   *
   *   \li A method to return the nonzero entries in a given row,
   *   together with the indices of the columns corresponding to the
   *   nonzero entries.  This method should have the following
   *   signature: \code void get_row_entries( const int irow,
   *   std::vector<double>& vals, std::vector<int>& ixs ) \endcode
   *
   *   @par More Examples
   *
   *   See FV2ndOrderTypes.h for an example of defining
   *   OpAssemblerSelector structs and see FVStaggeredSpatialOps.h for
   *   examples of defining assemblers for various operators.
   */
  template< typename OpType,
            typename SrcFieldT,
            typename DestFieldT >
  struct OpAssemblerSelector{};


  //====================================================================


  /**
   *  @class  SpatialOperator
   *  @author James C. Sutherland
   *  @date   December, 2006
   *
   *  Provides support for discrete spatial operators.
   *
   *  @par Template Parameters
   *  <ul>
   *
   *    <li> \b LinAlg The Linear Algebra type for use with this
   *    SpatialOperator.  This must define:
   *    <ul>
   *      <li> \b MatType This defines the type for the underlying Matrix operator.
   *
   *      <li> \b MatrixRow This defines the type for the
   *      representation of a Matrix row.
   *
   *      <li> \b column_iterator Iterates through columns.  The
   *      <c>column_iterator</c> must support increment and decrement
   *      operators, as well as the comparison <c>!=</c> operator and
   *      <code>double& value()</code> and <c>int index()</c> methods
   *      to obtain the index and value of the column coefficient.
   *
   *    </ul>
   *
   *    <li> \b OpType This defines the type of the operator.  This is
   *    really only used to distinguish various types of operators.
   *    In many cases, an empty struct will suffice for this type.
   *
   *    <li> \b SrcFieldT Specifies information about the field
   *    type that this operator acts on.  See SpatialField class for more information.
   *
   *    <li> \b DestFieldT Specifies information about the field
   *    type that this operator produces.  It must define the same
   *    things as the SrcFieldT type does.
   *  </ul>
   *
   *
   *  @par The OpAssemblerSelector
   * 
   *  All operators have an associated assembler object that must be
   *  instantiated and passed to the SpatialOperator constructor.  The
   *  full type of the OpAssemblerSelector is determined from the
   *  template parameters to the SpatialOperator.  However, the user
   *  must provide specializations of the OpAssemblerSelector concept
   *  that have the required functionality.  See documentation on
   *  OpAssemblerSelector for more details.
   */
  template< typename LinAlg,      // linear algebra support for this operator
            typename OpType,      // type of operator
            typename SrcFieldT,   // information on the source field (field operator acts on)
            typename DestFieldT > // information on the dest field (field operator produces)
  class SpatialOperator
  {
  public:

    typedef LinAlg                         LinAlgType;
    typedef typename LinAlg::MatType       MatType;
    typedef typename LinAlg::MatrixRow     MatrixRow;

    typedef typename LinAlg::column_iterator       column_iterator;
    typedef typename LinAlg::const_column_iterator const_column_iterator;

    typedef OpType                         Type;

    typedef SrcFieldT                      SrcFieldType;
    typedef DestFieldT                     DestFieldType;

    typedef typename SrcFieldT::Ghost      SrcGhost;
    typedef typename SrcFieldT::Location   SrcLocation;

    typedef typename DestFieldT::Ghost     DestGhost;
    typedef typename DestFieldT::Location  DestLocation;

    typedef typename OpAssemblerSelector
                       < OpType,
                         SrcFieldT,
                         DestFieldT >::Assembler        Assembler;

  public:

    /**
     *  Construct a SpatialOperator.
     *
     *  @param opAssembler The assembler is a strongly typed object
     *  that provides information required to construct this
     *  operator. The assembler type is defined by the
     *  OpAssemblerSelector, which must be specialized by the client
     *  who is building the SpatialOperator.  The assembler provides
     *  methods required to construct the operator.  See documentation
     *  on OpAssemblerSelector for more information.  Since the
     *  assembler is only required during construction, it may be a
     *  temporary object that can be destroyed after a call to this
     *  constructor.
     */
    SpatialOperator( Assembler & opAssembler );

    virtual ~SpatialOperator();


    /**
     *  @brief Obtain the requested row from this operator.
     *
     *  The <c>MatrixRow</c> type is defined by the LinAlg type.
     */
    inline MatrixRow get_row(const int irow) const{ return linAlg_.get_row(irow); }


    /**
     *  @name 
     *  Return the underlying linear algebra matrix representation 
     */
    //@{
    inline       MatType& get_linalg_mat()      { return mat_; }
    inline const MatType& get_linalg_mat() const{ return mat_; }
    //@}


    /**
     *  @brief Apply this operator to the supplied source field to
     *  produce the supplied destination field.
     *
     *  Calculates the matrix-vector product (dest)=[Op](src).
     *
     *  @param src  The field to apply the operator to.
     *  @param dest The resulting field.
     */
    inline void apply_to_field( const SrcFieldT& src,
                                DestFieldT& dest ) const;


    /**
     *  @brief Apply this operator to another operator to produce a third.
     *
     *  Calculates the matrix product [D] = [Op][S], where [D]
     *  represents the destination operator, [S] represents the
     *  "source" operator, and [Op] represents this operator.
     *
     *  @param src  The operator that we act on to produce the result.
     *  @param dest The resulting operator.
     */
    template< class SrcOp, class DestOp >
    inline void apply_to_op( const SrcOp& src, DestOp& dest ) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     *  @name Matrix Operators
     *
     *  Operators to add/subtract/assign another spatial operator with
     *  the same direction, and source/destination traits.
     */
    //@{
    template< typename OT >
    inline SpatialOperator& operator=( const SpatialOperator<LinAlg,OT,SrcFieldT,DestFieldT>& );

    template< typename OT >
    inline SpatialOperator& operator+=( const SpatialOperator<LinAlg,OT,SrcFieldT,DestFieldT>& );

    template< typename OT >
    inline SpatialOperator& operator-=( const SpatialOperator<LinAlg,OT,SrcFieldT,DestFieldT>& );
    //@}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     *  @name Field Operators
     *  Operators to add/subtract a SpatialField from the diagonal of this operator.
     */
    //@{
    /** @brief sum the given field into the diagonal */
    inline SpatialOperator& operator+=( const DestFieldT& f ){ linAlg_+=f; return *this; }

    /** @brief subtract the given field from the diagonal */
    inline SpatialOperator& operator-=( const DestFieldT& f ){ linAlg_-=f; return *this; }
    //@}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     *  Scale the matrix such that A(i,j) = x(i)*A(i,j) where i denotes
     *  the row number of A and j denotes the column number of A.
     */
    inline void left_scale( const DestFieldT& s ){ linAlg_.left_scale( s.get_linalg_vec() ); }

    /**
     *  Scale the matrix such that A(i,j) = x(j)*A(i,j) where i denotes
     *  the global row number of A and j denotes the global column
     *  number of A.
     */
    inline void right_scale( const SrcFieldT& a ){ linAlg_.right_scale( a.get_linalg_vec() ); }

    /** @brief reset the coefficients in the matrix */
    inline void reset_entries( const double val = 0 ){ linAlg_.reset_entries( val ); }

    /** @brief Obtain the number of rows in this operator */
    inline int nrows() const{ return nrows_; }

    /** @brief Obtain the number of columns in this operator */
    inline int ncols() const{ return ncols_; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    bool is_col_ghost( const int colnum ) const;
    bool is_row_ghost( const int rownum ) const;

    /** @brief Print the operator's entries to the specified output stream. */
    void Print( std::ostream & ) const;

    /**
     *  @brief Write the operator to a matlab file.
     *  @param prefix The name of the operator.  The file will be called "load_prefix.m".
     */
    void write_matlab( const std::string prefix ) const;


  protected:


    /**
     *  Check compatibility of this operator applied on the one provided
     *  as an argument.  The base class implementation simply performs
     *  dimensional compatibility checks.  Derived classes may
     *  specialize this method further.
     */
    template< typename T >
    inline bool compatibility_check( const T& op,
                                     const bool isResultOp ) const;

    /**
     * Check compatability of a field with this operator.  The base
     * class method simply checks for dimensional compatibility.
     * Derived classes may/should specialize this method to ensure
     * proper ghost availability.
     */
    template< typename FieldT >
    inline bool compatibility_check( const FieldT & field ) const;

    /** Insert a row into the matrix */
    inline void insert_row_entry( const int rownum,
                                  std::vector<double> & rowValues,
                                  std::vector<int> & rowIndices );

    /** Sum a row into the matrix */
    inline void sum_into_row( const int rownum,
                              std::vector<double> & rowValues,
                              std::vector<int> & rowIndices );

  private:

    const int nrows_, ncols_;
    std::set<int> ghostCols_, ghostRows_;
    LinAlg linAlg_;
    const int nNonZero_;
    MatType & mat_;
  };


  //====================================================================


  /**
   *  @class  SpatialOpDatabase
   *  @author James C. Sutherland
   *  @date   December, 2006
   *
   *  Factory for SpatialOperator objects.  These objects are registered
   *  here (created externally, ownership is transferred) and can be
   *  recalled for later use.
   *
   *  Note that one can have multiple operators registered for a given
   *  type, and activate them as needed.  This allows the potential for
   *  dynamic operator switching.
   *
   *  @todo Need to allow the ability to put multiple operators in a
   *  database and retrieve them, potentially with locking capability to prevent multiple access.
   */
  template< class SpatialOpType >
  class SpatialOpDatabase
  {
  public:

    typedef typename SpatialOpType::SrcGhost    SrcGhost;
    typedef typename SpatialOpType::DestGhost   DestGhost;


    static SpatialOpDatabase& self();

    /**
     *  Registers a new operator.
     *
     *  @param op  The operator itself.  Ownership is transfered.  This
     *  must be a heap-allocated object (build via "new")
     *
     *  @return id  A unique identifier for this operator.
     */
    int register_new_operator( SpatialOpType * const op );

    /**
     *  Obtain the spatial operator with the given id.
     *  Throws an exception if no match is found.
     *
     *  This pointer reference can change if the default operator for
     *  this type is changed via a call to
     *  <code>set_default_operator()</code> or
     *  <code>register_new_operator()</code>.
     */
    SpatialOpType* retrieve_operator( const int id=-1 ) const;


    /**
     *  determine if an operator with this id exists in the databse
     */
    bool query_operator( const int id=-1 ) const;

    /**
     *  Remove all registered operators from the database.
     */
    void empty_database();

  private:

    SpatialOpDatabase();
    ~SpatialOpDatabase();

    SpatialOpDatabase(const SpatialOpDatabase&);
    SpatialOpDatabase&  operator=(const SpatialOpDatabase&);

    typedef std::map< int, SpatialOpType* > OpMap;
    OpMap opMap_;

    int idCounter_;
  };


  //====================================================================









  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









  //==================================================================


  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  SpatialOperator( Assembler & opAssembler )
    : nrows_ ( opAssembler.get_nrows()  ),
      ncols_ ( opAssembler.get_ncols()  ),
      nNonZero_( opAssembler.num_nonzeros() ),
      mat_( linAlg_.setup_matrix( nrows_, ncols_, nNonZero_ ) )
  {
    // build the operator
    std::vector<double> vals( nNonZero_, 0.0 );
    std::vector<int>    ixs ( nNonZero_, 0   );
    for( int i=0; i<nrows_; ++i ){
      vals.clear();
      ixs.clear();
      opAssembler.get_row_entries( i, vals, ixs );
      insert_row_entry( i, vals, ixs );
    }
    opAssembler.get_ghost_cols( ghostCols_ );
    opAssembler.get_ghost_rows( ghostRows_ );
    linAlg_.finalize();
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  ~SpatialOperator()
  {
    linAlg_.destroy_matrix();
  }
  //------------------------------------------------------------------
  /**
   *  Apply this SpatialOperator to a field.  The source and
   *  destination fields must have compatible dimension and ghosting
   *  for use with this operator.
   */
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  void
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  apply_to_field( const SrcFieldT& src,
                  DestFieldT& dest ) const
  {
    linAlg_.multiply( src.get_linalg_vec(), dest.get_linalg_vec() );
  }

  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  template< class SrcOp, class DestOp >
  void
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  apply_to_op( const SrcOp& src, DestOp& dest ) const
  {
    BOOST_STATIC_ASSERT( bool( IsSameType<SrcGhost,  typename SrcOp::DestGhost>::result ) );
    BOOST_STATIC_ASSERT( bool( IsSameType<DestGhost, typename SrcOp::SrcGhost >::result ) );
    assert( compatibility_check(src, false) );
    assert( compatibility_check(dest,true ) );
    linAlg_.multiply( src.get_linalg_mat(), dest.get_linalg_mat() );
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  template< typename OT >
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>&
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  operator=( const SpatialOperator<LinAlg,OT,SrcFieldT,DestFieldT>& op )
  {
    linAlg_ = op;
    return *this;
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  template< typename OT >
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>&
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  operator+=( const SpatialOperator<LinAlg,OT,SrcFieldT,DestFieldT>& op )
  {
    linAlg_ += op;
    return *this;
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  template< typename OT >
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>&
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  operator-=( const SpatialOperator<LinAlg,OT,SrcFieldT,DestFieldT>& op )
  {
    linAlg_ -= op;
    return *this;
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  bool
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  is_row_ghost( const int irow ) const
  {
    return ( ghostRows_.find( irow ) != ghostRows_.end() );
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  bool
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  is_col_ghost( const int irow ) const
  {
    return ( ghostCols_.find( irow ) != ghostCols_.end() );
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  template< typename T >
  bool
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  compatibility_check( const T& op, const bool isResultOp ) const
  {
    using std::endl;
    using std::cout;

    if( isResultOp ){
      if( nrows_ != op.nrows() ){
        cout << "ERROR: Destination matrix must have same number of rows as operator." << endl
             << "  Dest [nr,nc] = [" << op.nrows() << "," << op.ncols() << "]" << endl
             << "  Op   [nr,nc] = [" << nrows_     << "," << ncols_     << "]" << endl;
        return false;
      }
    }
    else{
      if( ncols_ != op.nrows() ){
        cout << "ERROR: Source matrix must have same number of rows as operator has colums." << endl
             << "  Dest [nr,nc] = [" << op.nrows() << "," << op.ncols() << "]" << endl
             << "  Op   [nr,nc] = [" << nrows_     << "," << ncols_     << "]" << endl;
        return false;
      }
    }
    return true;
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  void
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  insert_row_entry( const int rownum,
                    std::vector<double> & rowValues,
                    std::vector<int> & rowIndices )
  {
    linAlg_.insert_row_values( rownum,
                               rowValues,
                               rowIndices );
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  void
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  sum_into_row( const int rownum,
                std::vector<double> & rowValues,
                std::vector<int> & rowIndices )
  {
    linAlg_.sum_in_row_values( rownum,
                               rowValues.size(),
                               &rowValues[0],
                               &rowIndices[0] );
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  void
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  Print( std::ostream & s ) const
  {
    linAlg_.print_mat( s );
  }
  //------------------------------------------------------------------
  template< typename LinAlg, typename OpType, typename SrcFieldT, typename DestFieldT >
  void
  SpatialOperator<LinAlg,OpType,SrcFieldT,DestFieldT>::
  write_matlab( const std::string prefix ) const
  {
    const std::string fname = "load_"+prefix +".m";
    std::ofstream fout( fname.c_str() );
    fout << "function A = load_" << prefix << "()" << std::endl;

    // first time through count nonzeros for preallocation in matlab.
    // second time through, write the entries.
    for( int writeLoop=0; writeLoop<=1; ++writeLoop ){
      int i=0;
      for( int irow=0; irow<nrows_; ++irow ){
        typename LinAlg::MatrixRow row = linAlg_.get_row( irow );
        for( typename LinAlg::column_iterator icol = row.begin();
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
  //------------------------------------------------------------------

  //==================================================================


  //------------------------------------------------------------------
  template< class SpatialOpType >
  SpatialOpDatabase<SpatialOpType>&
  SpatialOpDatabase<SpatialOpType>::self()
  {
    static SpatialOpDatabase<SpatialOpType> s;
    return s;
  }
//--------------------------------------------------------------------
  template< class SpatialOpType >
  int
  SpatialOpDatabase<SpatialOpType>::
  register_new_operator( SpatialOpType * const op )
  {
    const int id = ++idCounter_;
    std::pair< typename OpMap::const_iterator, bool > result
      = opMap_.insert( std::make_pair(id,op) );
    assert( result.second );
    return id;
  }
  //------------------------------------------------------------------
  template< class SpatialOpType >
  SpatialOpType*
  SpatialOpDatabase<SpatialOpType>::
  retrieve_operator( const int id ) const
  {
    using std::endl;

    typename OpMap::const_iterator iop = opMap_.begin();

    if( id==-1 ){
      if( opMap_.size() > 1 ){
        std::ostringstream msg;
        msg << "ERROR!  You must provide a unique identifier, since multiple operators have been registered." << endl
            << "        registered ids:" << endl;
        for( iop=opMap_.begin(); iop!=opMap_.end(); ++iop ){
          msg << "     " << iop->first << std::endl;
        }
        throw std::runtime_error(msg.str());
      }
    }
    else{
      iop = opMap_.find( id );
    }
    if( iop == opMap_.end() ){
      std::ostringstream msg;
      msg << "ERROR!  Attempted to retrieve an operator that does not exist." << std::endl
          << "        Operator type name: " << typeid(SpatialOpType).name() << std::endl
          << "   Registered ids:" << std::endl;
      for( iop=opMap_.begin(); iop!=opMap_.end(); ++iop ){
        msg << "     " << iop->first <<  std::endl;
      }
      throw std::runtime_error( msg.str() );
    }
    return iop->second;
  }
  //------------------------------------------------------------------
  template< class SpatialOpType >
  void
  SpatialOpDatabase<SpatialOpType>::
  empty_database()
  {
    for( typename OpMap::iterator jj=opMap_.begin(); jj!=opMap_.end(); ++jj ){
      delete jj->second;
    }
  }
  //------------------------------------------------------------------
  template< class SpatialOpType >
  bool
  SpatialOpDatabase<SpatialOpType>::
  query_operator( const int id ) const
  {
    if( id==-1 ){
      if( opMap_.size() > 1 ){
        std::ostringstream msg;
        msg << "ERROR!  You must provide a unique identifier, since multiple operators have been registered." << std::endl
            << "        registered ids:" << std::endl;
        for( typename OpMap::const_iterator iop=opMap_.begin(); iop!=opMap_.end(); ++iop ){
          msg << "     " << iop->first <<  std::endl;
        }
        throw std::runtime_error(msg.str());
      }
      return opMap_.size() == 1;
    }
    return( opMap_.find(id) != opMap_.end() );
  }
  //------------------------------------------------------------------
  template< class SpatialOpType >
  SpatialOpDatabase<SpatialOpType>::
  SpatialOpDatabase()
  {
    idCounter_=0;
  }
  //------------------------------------------------------------------
  template< class SpatialOpType >
  SpatialOpDatabase<SpatialOpType>::
  ~SpatialOpDatabase()
  {
    empty_database();
  }
  //------------------------------------------------------------------


  //==================================================================


} // namespace SpatialOps


#endif
