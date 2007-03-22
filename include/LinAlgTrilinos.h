#ifndef LinAlgTrilinos_h
#define LinAlgTrilinos_h

#include <vector>

#include <SpatialOperator.h>
#include <SpatialField.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>

namespace SpatialOps{

  //====================================================================

  class LinAlgTrilinos
  {

  public:

    typedef Epetra_Vector    VecType;
    typedef Epetra_CrsMatrix MatType;

    LinAlgTrilinos();
    ~LinAlgTrilinos();

    VecType& setup_vector( const int npts,
			   double* const values );
    
    MatType& setup_matrix( const int nrows,
			   const int ncols,
			   const int entriesPerRow );

    void insert_row_values( const int rownum,
			    std::vector<double> & rowValues,
			    std::vector<int> & rowIndices );

    void destroy_matrix();
    void destroy_vector();

    //@{  /** Operators to assemble matrices into this matrix */

    template< typename OpType,
	      typename Direction,
	      typename SrcFieldTraits,
	      typename DestFieldTraits >
    inline LinAlgTrilinos& operator= ( const SpatialOperator< LinAlgTrilinos, OpType, Direction, SrcFieldTraits, DestFieldTraits > & m );

    template< typename OpType,
	      typename Direction,
	      typename SrcFieldTraits,
	      typename DestFieldTraits >
    inline LinAlgTrilinos& operator+=( const SpatialOperator< LinAlgTrilinos, OpType, Direction, SrcFieldTraits, DestFieldTraits > & m );

    template< typename OpType,
	      typename Direction,
	      typename SrcFieldTraits,
	      typename DestFieldTraits >
    inline LinAlgTrilinos& operator-=( const  SpatialOperator< LinAlgTrilinos, OpType, Direction, SrcFieldTraits, DestFieldTraits > & m );

    //}@


    //@{ /** Operators to assemble vectors into this matrix.  We only allow assembly of fields with consistent spatial storage location. */

    template<typename FieldT>
    inline LinAlgTrilinos& operator+=( const FieldT & f );

    template<typename FieldT>
    inline LinAlgTrilinos& operator-=( const FieldT& f );

    //}@

    void  left_scale( const VecType& v );
    void right_scale( const VecType& v );

    // reset entries in the matrix
    void reset_entries( const double val );


    void multiply( const MatType& B, MatType& C ) const;
    void multiply( const VecType& B, VecType& C ) const;

    void finalize();

  private:

    const Epetra_LocalMap *colMap_, *rowMap_;

    VecType* vec_;
    MatType* mat_;

  };




  //====================================================================





  template< typename OpType, typename Direction, typename SrcFieldTraits, typename DestFieldTraits >
  LinAlgTrilinos&
  LinAlgTrilinos::operator=( const SpatialOperator<LinAlgTrilinos,OpType,Direction,SrcFieldTraits,DestFieldTraits>& m )
  {
    // jcs this is INEFFICIENT.  It would be better to do this directly.
    reset_entries(0.0);
    return ( (*this)+=m );
  }
  //------------------------------------------------------------------
  template< typename OpType, typename Direction, typename SrcFieldTraits, typename DestFieldTraits >
  LinAlgTrilinos&
  LinAlgTrilinos::operator+=( const SpatialOperator<LinAlgTrilinos,OpType,Direction,SrcFieldTraits,DestFieldTraits>& m )
  {
    EpetraExt::MatrixMatrix::Add( m.get_linalg_mat(), false, 1.0, *mat_, 1.0 );
    return *this;
  }
  //------------------------------------------------------------------
  template< typename OpType, typename Direction, typename SrcFieldTraits, typename DestFieldTraits >
  LinAlgTrilinos&
  LinAlgTrilinos::operator-=( const SpatialOperator<LinAlgTrilinos,OpType,Direction,SrcFieldTraits,DestFieldTraits>& m )
  {
    EpetraExt::MatrixMatrix::Add( m.get_linalg_mat(), false, -1.0, *mat_, 1.0 );
    return *this;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  LinAlgTrilinos&
  LinAlgTrilinos::operator+=( const FieldT& f )
  {
    const int nrows = mat_->NumGlobalRows();
    assert( nrows == f.get_ntotal() );
    const double * const fptr = f.get_ptr();
    for( int i=0; i<nrows; ++i ){
      double val = fptr[i];
      mat_->SumIntoMyValues( i, 1, &val, &i );
    }
    return *this;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  LinAlgTrilinos&
  LinAlgTrilinos::operator-=( const FieldT & f )
  {
    const int nrows = mat_->NumGlobalRows();
    assert( nrows == f.get_ntotal() );
    const double * const fptr = f.get_ptr();
    for( int i=0; i<nrows; ++i ){
      double val = -fptr[i];
      mat_->SumIntoMyValues( i, 1, &val, &i );
    }
    return *this;
  }
  //------------------------------------------------------------------


} // namespace SpatialOps

#endif
