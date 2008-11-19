#ifndef LinAlgTrilinos_h
#define LinAlgTrilinos_h

#include <vector>

#include <spatialops/SpatialOperator.h>
#include <spatialops/SpatialField.h>

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <Epetra_Vector.h>

class Epetra_LocalMap;  // forward declaration


namespace SpatialOps{

  //====================================================================

  class LinAlgTrilinos
  {

  public:

    typedef Epetra_Vector    VecType;
    typedef Epetra_CrsMatrix MatType;

    class column_iterator;
    class const_column_iterator;

    class MatrixRow
    {
    public:
      MatrixRow( int rowNum, MatType& mat )
      {
        assert(rowNum<mat.NumMyRows());
        mat.ExtractMyRowView( rowNum, ncols_, values_, indices_ );
      }
      ~MatrixRow(){}

      double& operator()(const int i)      { return values_[i]; }
      double  operator()(const int i) const{ return values_[i]; }

      double& operator[](const int i)      { return values_[i]; }
      double  operator[](const int i) const{ return values_[i]; }


      column_iterator       begin()      { return column_iterator(*this,ncols_,0); }
      const_column_iterator begin() const{ return const_column_iterator(*this,ncols_,0); }

      column_iterator end()      { return column_iterator(*this,ncols_,ncols_); }
      const_column_iterator end() const{ return const_column_iterator(*this,ncols_,ncols_); }

      inline double*& get_column_vals()      { return values_;  }
      inline double*  get_column_vals() const{ return values_;}

      inline int*     get_column_ixs() const{ return indices_; }
      
    private:
      int ncols_;
      double* values_;
      int* indices_;
    };

    class column_iterator{
    public:
      column_iterator( MatrixRow& mr, const int ncol, const int i=0)
        : mr_(mr), ncol_(ncol), i_(i)
      {}

      inline column_iterator& operator++(){assert(i_<ncol_); ++i_; return *this;}
      inline column_iterator& operator--(){assert(i_>0);     --i_; return *this; }

      inline bool operator!=(const column_iterator c) const{return c.i_!=i_;}
      inline bool operator==(const column_iterator c) const{return c.i_==i_;}

      inline unsigned int index() const{return mr_.get_column_ixs()[i_]; }

      double& operator*()      {return mr_.get_column_vals()[i_]; }

    private:
      MatrixRow& mr_;
      const int ncol_;
      int i_;
    };

    class const_column_iterator{
    public:
      const_column_iterator( const MatrixRow& mr, const int ncol, const int i=0)
        : mr_(mr), ncol_(ncol), i_(i)
      {}

      inline const_column_iterator& operator++(){assert(i_<ncol_); ++i_; return *this;}
      inline const_column_iterator& operator--(){assert(i_>0);     --i_; return *this; }

      inline bool operator!=(const const_column_iterator c) const{return c.i_!=i_;}

      inline unsigned int index() const{return mr_.get_column_ixs()[i_]; }
      double operator*() const{return mr_.get_column_vals()[i_]; }

    private:
      const MatrixRow& mr_;
      const int ncol_;
      int i_;

      friend class column_iterator;
    };

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

    inline MatrixRow get_row( const int irow ) const{ return MatrixRow(irow,*mat_); }

    /**
     *  @name Matrix Assembly Operators
     *
     *  Operators to assemble matrices into this matrix
     */
    //@{

    template< typename OpType,
              typename SrcFieldT,
              typename DestFieldT >
    inline LinAlgTrilinos& operator= ( const SpatialOperator< LinAlgTrilinos, OpType, SrcFieldT, DestFieldT > & m );

    template< typename OpType,
              typename SrcFieldT,
              typename DestFieldT >
    inline LinAlgTrilinos& operator+=( const SpatialOperator< LinAlgTrilinos, OpType, SrcFieldT, DestFieldT > & m );

    template< typename OpType,
              typename SrcFieldT,
              typename DestFieldT >
    inline LinAlgTrilinos& operator-=( const  SpatialOperator< LinAlgTrilinos, OpType, SrcFieldT, DestFieldT > & m );

    //@}


    /**
     *  @name Vector Assembly Operators
     *
     *  Operators to assemble vectors into this matrix.  We only allow
     *  assembly of fields with consistent spatial storage location.
     */
    //@{
    template<typename FieldT>
    inline LinAlgTrilinos& operator+=( const FieldT & f );

    template<typename FieldT>
    inline LinAlgTrilinos& operator-=( const FieldT& f );
    //@}

    void  left_scale( const VecType& v );
    void right_scale( const VecType& v );

    // reset entries in the matrix
    void reset_entries( const double val );


    void multiply( const MatType& B, MatType& C ) const;
    void multiply( const VecType& B, VecType& C ) const;

    void finalize();

    void print_vec( std::ostream& s ) const;
    void print_mat( std::ostream& s ) const;

  private:

    const Epetra_LocalMap *colMap_, *rowMap_;

    VecType* vec_;
    MatType* mat_;

  };




  //====================================================================




  template< typename OpType, typename SrcFieldT, typename DestFieldT >
  LinAlgTrilinos&
  LinAlgTrilinos::operator=( const SpatialOperator<LinAlgTrilinos,OpType,SrcFieldT,DestFieldT>& m )
  {
    /**
     *  @todo improve efficiency of operator assignment.  See \c LinAlgTrilinos::operator=
     */
    // jcs this is INEFFICIENT.  It would be better to do this directly.
    reset_entries(0.0);
    return ( (*this)+=m );
  }
  //------------------------------------------------------------------
  template< typename OpType, typename SrcFieldT, typename DestFieldT >
  LinAlgTrilinos&
  LinAlgTrilinos::operator+=( const SpatialOperator<LinAlgTrilinos,OpType,SrcFieldT,DestFieldT>& m )
  {
    EpetraExt::MatrixMatrix::Add( m.get_linalg_mat(), false, 1.0, *mat_, 1.0 );
    return *this;
  }
  //------------------------------------------------------------------
  template< typename OpType, typename SrcFieldT, typename DestFieldT >
  LinAlgTrilinos&
  LinAlgTrilinos::operator-=( const SpatialOperator<LinAlgTrilinos,OpType,SrcFieldT,DestFieldT>& m )
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
    int i=0;
    for( typename FieldT::const_iterator ifld = f.begin(); ifld!=f.end(); ++ifld, ++i ){
      double val = *ifld;
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
    int i=0;
    for( typename FieldT::const_iterator ifld = f.begin(); ifld!=f.end(); ++ifld, ++i ){
      double val = -*ifld;
      mat_->SumIntoMyValues( i, 1, &val, &i );
    }
    return *this;
  }
  //------------------------------------------------------------------


} // namespace SpatialOps

#endif
