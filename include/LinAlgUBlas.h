#ifndef LinAlgUBlasBoost_h
#define LinAlgUBlasBoost_h

#include <map>
#include <vector>

#include <SpatialOperator.h>
#include <SpatialField.h>

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace SpatialOps{

  /**
   *  @class  LinAlgUBlas
   *  @author James C. Sutherland
   *  @date   Jan, 2008
   *  @brief  Provides interface between SpatialOps and Boost's uBlas library.
   */
  class LinAlgUBlas
  {
  public:

    typedef boost::numeric::ublas::compressed_matrix< double, boost::numeric::ublas::row_major > MatType;
    typedef boost::numeric::ublas::shallow_array_adaptor<double> ArrayAdaptor;
    typedef boost::numeric::ublas::vector<double,ArrayAdaptor>  VecType;

    typedef boost::numeric::ublas::matrix_row<MatType> MatrixRow;
    typedef MatrixRow::iterator column_iterator;
    typedef MatrixRow::const_iterator const_column_iterator;

    LinAlgUBlas();

    ~LinAlgUBlas();

    VecType& setup_vector( const int npts,
			   double* const values );
    
    MatType& setup_matrix( const int nrows,
			   const int ncols,
			   const int entriesPerRow );

    void insert_row_values( const int rownum,
			    std::vector<double> & rowValues,
			    std::vector<int> & rowIndices );

    inline MatrixRow get_row(const int irow) const{ return MatrixRow(*mat_,irow); }

    void destroy_matrix();
    void destroy_vector();


    /**
     *  @name Matrix Assembly Operators
     *
     *  Operators to assemble matrices into this matrix
     */
    //@{

    template< typename OpType,
	      typename SrcFieldT,
	      typename DestFieldT >
    inline LinAlgUBlas& operator= ( const SpatialOperator< LinAlgUBlas, OpType, SrcFieldT, DestFieldT > & m ){
      *mat_ = m.get_linalg_mat();
      return *this;
    }

    template< typename OpType,
	      typename SrcFieldT,
	      typename DestFieldT >
    inline LinAlgUBlas& operator+=( const SpatialOperator< LinAlgUBlas, OpType, SrcFieldT, DestFieldT > & m ){
      *mat_ += m.get_linalg_mat();
      return *this;
    }

    template< typename OpType,
	      typename SrcFieldT,
	      typename DestFieldT >
    inline LinAlgUBlas& operator-=( const  SpatialOperator< LinAlgUBlas, OpType, SrcFieldT, DestFieldT > & m ){
      *mat_ -= m.get_linalg_mat();
      return *this;
    }

    //@}


    /**
     *  @name Vector Assembly Operators
     *
     *  Operators to assemble vectors into this matrix.  We only allow
     *  assembly of fields with consistent spatial storage location.
     */
    //@{
    template<typename FieldT>
    inline LinAlgUBlas& operator+=( const FieldT & f );

    template<typename FieldT>
    inline LinAlgUBlas& operator-=( const FieldT& f );
    //@}

    void left_scale( const VecType& v );
    void right_scale( const VecType& v );

    // reset entries in the matrix
    void reset_entries( const double val );


    inline void multiply( const MatType& B, MatType& C ) const;
    inline void multiply( const VecType& b, VecType& c ) const;

    void finalize();

    void print_vec( std::ostream& s ) const;
    void print_mat( std::ostream& s ) const;

  private:

    VecType* vec_;
    MatType* mat_;
  };


  //====================================================================







  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementations
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







  //====================================================================


  //--------------------------------------------------------------------

  template<typename FieldT>
  LinAlgUBlas&
  LinAlgUBlas::operator+=( const FieldT & f )
  {
    const int nrows = mat_->size1();
    assert( nrows == f.get_ntotal() );
    int i=0;
    for( typename FieldT::const_iterator ifld = f.begin(); ifld!=f.end(); ++ifld, ++i ){
      (*mat_)(i,i) += *ifld;
    }
    return *this;
  }
  //--------------------------------------------------------------------

  template<typename FieldT>
  LinAlgUBlas&
  LinAlgUBlas::operator-=( const FieldT & f )
  {
    const int nrows = mat_->size1();
    assert( nrows == f.get_ntotal() );
    int i=0;
    for( typename FieldT::const_iterator ifld = f.begin(); ifld!=f.end(); ++ifld, ++i ){
      (*mat_)(i,i) -= *ifld;
    }
    return *this;
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::multiply( const VecType& src, VecType& dest ) const
  {
    // "noalias" improves performance.  This requires that src and dest are different arrays.
    //noalias(dest) = prod(*mat_,src);

    // alternatively:
    boost::numeric::ublas::axpy_prod( *mat_, src, dest, true );
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::multiply( const MatType& B, MatType& C ) const
  {
    noalias(C) = prod(*mat_,B);

    // alternatively:
    //boost::numeric::ublas::axpy_prod( *mat_, B, C, true );
  }

  //--------------------------------------------------------------------


} // namespace SpatialOps


#endif
