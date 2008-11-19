#include <spatialops/LinAlgUBlas.h>

#include <boost/numeric/ublas/io.hpp>

namespace SpatialOps{


  LinAlgUBlas::LinAlgUBlas()
  {
    vec_ = NULL;
    mat_ = NULL;
  }

  //--------------------------------------------------------------------

  LinAlgUBlas::~LinAlgUBlas()
  {
    destroy_matrix();
    destroy_vector();
  }
  
  //--------------------------------------------------------------------

  void
  LinAlgUBlas::destroy_vector()
  {
    if( vec_ != NULL ){
      delete vec_;
      vec_ = NULL;
    }
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::destroy_matrix()
  {
    if( mat_ != NULL ){
      delete mat_;
      mat_ = NULL;
    }
  }

  //--------------------------------------------------------------------

  LinAlgUBlas::VecType&
  LinAlgUBlas::setup_vector( const int npts,
                             double* const fieldValues )
  {
    assert( vec_ == NULL );
    vec_ = new VecType( npts, ArrayAdaptor(npts,fieldValues) );
    const double* ival = fieldValues;
    for( VecType::iterator i=vec_->begin(); i!=vec_->end(); ++i, ++ival )
      *i = *ival;
    return *vec_;
  }

  //--------------------------------------------------------------------

  LinAlgUBlas::MatType&
  LinAlgUBlas::setup_matrix( const int nrows, const int ncols, const int entriesPerRow )
  {
    assert( mat_ == NULL );
    mat_ = new MatType( nrows, ncols, entriesPerRow );
    return *mat_;
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::insert_row_values( const int rownum,
                                  std::vector<double> & rowValues,
                                  std::vector<int> & rowIndices )
  {
    MatrixRow row( *mat_, rownum );
    std::vector<double>::const_iterator ival = rowValues.begin();
    const std::vector<double>::const_iterator ivale = rowValues.end();
    std::vector<int>::const_iterator irow = rowIndices.begin();
    for( ; ival!=ivale; ++ival, ++irow )
      row(*irow) = *ival;
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::finalize(){}

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::reset_entries( const double val )
  {
    typedef MatType::iterator1 RowIter;
    typedef MatType::iterator2 ColIter;
    for( RowIter irow = mat_->begin1(); irow!=mat_->end1(); ++irow ){
      for( ColIter icol=irow.begin(); icol!=irow.end(); ++icol ){
        *icol = val;
      }
    }
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::left_scale( const VecType& v )
  {
    typedef VecType::const_iterator VecIter;
    typedef MatType::iterator1 RowIter;
    typedef MatType::iterator2 ColIter;
    VecIter ivec = v.begin();
    for( RowIter irow = mat_->begin1(); irow!=mat_->end1(); ++irow, ++ivec ){
      for( ColIter icol=irow.begin(); icol!=irow.end(); ++icol ){
        *icol *= *ivec;
      }
    }
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::right_scale( const VecType& v )
  {
    typedef VecType::const_iterator VecIter;
    typedef MatType::iterator1 RowIter;
    typedef MatType::iterator2 ColIter;
    for( RowIter irow = mat_->begin1(); irow!=mat_->end1(); ++irow ){
      VecIter ivec = v.begin();
      for( ColIter icol=irow.begin(); icol!=irow.end(); ++icol, ++ivec ){
        *icol *= *ivec;
      }
    }
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::print_vec( std::ostream& s ) const
  {
    s << *vec_ << std::endl;
  }

  //--------------------------------------------------------------------

  void
  LinAlgUBlas::print_mat( std::ostream& s ) const
  {
    s << *mat_ << std::endl;
  }

  //--------------------------------------------------------------------


}  // namespace SpatialOps
