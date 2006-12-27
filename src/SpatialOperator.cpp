#include <SpatialOperator.h>
#include <SpatialField.h>

//-----------------------------
#include <numeric>
#include <iostream>
// #include <sstream>
// #include <stdexcept>
//-----------------------------

//-----------------------------
// Trilinos includes
#include <Epetra_LocalMap.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>
// Trilinos includes
//-----------------------------


namespace SpatialOps{


//====================================================================


SpatialOperator::SpatialOperator( const int nrows,
				  const int ncols,
				  const int nghost,
				  const int entriesPerRow,
				  const std::vector<int> & extent )
  : extent_( extent ),
    nrows_ ( nrows  ),
    ncols_ ( ncols  ),
    nghost_( nghost ),

    entriesPerRow_( entriesPerRow ),


    isFinalized_( false ),

    mat_( NULL )
{
  //
  // build the appropriate trilinos objects to hold this sparse matrix.
  //

  // Step 1: Build the Epetra_Map objects to describe the matrix layout
  //         The usage of Epetra_LocalMap ensures that the entire
  //         matrix resides on a single processor.

  rowMap_ = &MapFactory::self().get_map( nrows_ );
  colMap_ = &MapFactory::self().get_map( ncols_ );

  // Step 2: Build the Epetra_CrsMatrix object
  mat_ = new Epetra_CrsMatrix( Copy, *rowMap_, *colMap_, entriesPerRow, true );
}
//--------------------------------------------------------------------
SpatialOperator::~SpatialOperator()
{
  delete mat_;
}
//--------------------------------------------------------------------
void
SpatialOperator::finalize()
{
  mat_->FillComplete( *colMap_, *rowMap_ );
  mat_->OptimizeStorage();
  isFinalized_ = true;
}
//--------------------------------------------------------------------
void
SpatialOperator::apply( const SpatialField & srcField,
			SpatialField & destField ) const
{
  assert( mat_ != NULL );
  assert( ready() );
  assert( compatibility_check(  srcField, SOURCE_FIELD ) );
  assert( compatibility_check( destField,   DEST_FIELD ) );

  mat_->Multiply( false, srcField.epetra_vec(), destField.epetra_vec() );
}
//--------------------------------------------------------------------
void
SpatialOperator::apply( const SpatialOperator & B,
			SpatialOperator & C ) const
{
  assert( ready() && B.ready() && C.ready() );

  assert( compatibility_check(B) );
  assert( compatibility_check(C) );

  const int flag =
    EpetraExt::MatrixMatrix::Multiply( epetra_mat(), false,
				       B.epetra_mat(), false,
				       C.epetra_mat() );
  if( flag!=0 ) std::cout << flag << std::endl;
  assert( flag==0 );
}
//--------------------------------------------------------------------
SpatialOperator&
SpatialOperator::operator += ( const SpatialOperator& op )
{
  assert( op.nrows() == nrows() );
  assert( op.ncols() == ncols() );

  EpetraExt::MatrixMatrix::Add( op.epetra_mat(), false,
				1.0,
				epetra_mat(),
				1.0 );
  return *this;

  // hand-coded method - only works for operators with identical sparsity patterns:
  for( int i=0; i<op.nrows(); ++i ){
    double * destRowVals=0;
    int * destRowIx=0;
    int nvals=0;
    epetra_mat().ExtractMyRowView( i, nvals, destRowVals, destRowIx );

    double * srcRowVals=0;
    int * srcRowIx=0;
    int nsrcVals=0;
    op.epetra_mat().ExtractMyRowView( i, nsrcVals, srcRowVals, srcRowIx );

    assert( nvals==nsrcVals );

    for( int k=0; k<nsrcVals; ++k ){
      assert( srcRowIx[k] == destRowIx[k] );
      *destRowVals += *srcRowVals;
      ++destRowVals; ++srcRowVals;
    }
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialOperator&
SpatialOperator::operator -= ( const SpatialOperator& op )
{
  assert( op.nrows() == nrows() );
  assert( op.ncols() == ncols() );

  EpetraExt::MatrixMatrix::Add( op.epetra_mat(), false,
				-1.0,
				epetra_mat(),
				1.0 );
  return *this;

  // hand-coded method - only works for operators with identical sparsity patterns:
  for( int i=0; i<nrows(); ++i ){
    double * destRowVals=0;
    int * destRowIx=0;
    int nvals=0;
    epetra_mat().ExtractMyRowView( i, nvals, destRowVals, destRowIx );

    double * srcRowVals=0;
    int * srcRowIx=0;
    int nsrcVals=0;
    op.epetra_mat().ExtractMyRowView( i, nsrcVals, srcRowVals, srcRowIx );

    assert( nvals==nsrcVals );

    for( int k=0; k<nsrcVals; ++k ){
      assert( srcRowIx[k] == destRowIx[k] );
      *destRowVals -= *srcRowVals;
      ++destRowVals; ++srcRowVals;
    }
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialOperator&
SpatialOperator::operator += ( const SpatialField & f )
{
  assert( nrows() == f.get_ntotal() );

  const double * const fptr = f.get_ptr();
  for( int i=0; i<nrows(); ++i ){
    double val = fptr[i];
    mat_->SumIntoMyValues( i, 1, &val, &i );
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialOperator&
SpatialOperator::operator -= ( const SpatialField & f )
{
  assert( nrows() == f.get_ntotal() );

  const double * const fptr = f.get_ptr();
  for( int i=0; i<nrows(); ++i ){
    double val = -fptr[i];
    mat_->SumIntoMyValues( i, 1, &val, &i );
  }
  return *this;
}
//--------------------------------------------------------------------
void
SpatialOperator::left_scale( const SpatialField& f )
{
  epetra_mat().LeftScale( f.epetra_vec() );
}
//--------------------------------------------------------------------
void
SpatialOperator::right_scale( const SpatialField& f )
{
  epetra_mat().RightScale( f.epetra_vec() );
}
//--------------------------------------------------------------------
void
SpatialOperator::zero_entries()
{
  for( int i=0; i<nrows(); ++i ){
    double * destRowVals=0;
    int * destRowIx=0;
    int nvals=0;
    epetra_mat().ExtractMyRowView( i, nvals, destRowVals, destRowIx );

    for( int k=0; k<nvals; ++k ){
      *destRowVals = 0.0;
      ++destRowVals;
    }
  }
}
//--------------------------------------------------------------------
void
SpatialOperator::insert_row_entry( const int rownum,
				   std::vector<double> & rowValues,
				   std::vector<int> & rowIndices )
{
  const int flag = mat_->InsertGlobalValues( rownum,
					     rowValues.size(),
					     &rowValues[0],
					     &rowIndices[0] );
  if( flag!=0 ) std::cout << flag << std::endl;
  assert( flag==0 );
}
//--------------------------------------------------------------------
void
SpatialOperator::sum_into_row( const int rownum,
			       std::vector<double> & rowValues,
			       std::vector<int> & rowIndices )
{
  int flag = mat_->SumIntoGlobalValues( rownum,
					rowValues.size(),
					&rowValues[0],
					&rowIndices[0] );
  if( flag!=0 ) std::cout << flag << std::endl;
  assert( flag==0 );
}
//--------------------------------------------------------------------
Epetra_CrsMatrix &
SpatialOperator::epetra_mat()
{
 return *mat_;
}
//--------------------------------------------------------------------
const Epetra_CrsMatrix &
SpatialOperator::epetra_mat() const
{
  return *mat_;
}
//--------------------------------------------------------------------
bool
SpatialOperator::compatibility_check( const SpatialOperator& op  ) const
{
  if( nrows_         != op.nrows_         ) return false;
  if( ncols_         != op.ncols_         ) return false;
  //  if( entriesPerRow_ != op.entriesPerRow_ ) return false;
  return true;
}
//--------------------------------------------------------------------
bool
SpatialOperator::compatibility_check( const SpatialField & field,
				      const FieldType fldType ) const
{
  switch( fldType ){
  case SOURCE_FIELD:
    if( ncols_ != field.epetra_vec().GlobalLength() ){
      std::cout << "expecting " << ncols_ << " entries, found "
		<< field.epetra_vec().GlobalLength() << std::endl;
      return false;
    }
    break;
  case DEST_FIELD:
    if( nrows_ != field.epetra_vec().GlobalLength() ){
      std::cout << "expecting " << nrows_ << " entries, found "
		<< field.epetra_vec().GlobalLength() << std::endl;
      return false;
    }
    break;
  }
  return true;
}
//--------------------------------------------------------------------


//====================================================================


//---------------------------------------------------------
MapFactory::MapFactory()
  : com_( new Epetra_SerialComm() )
{}
//---------------------------------------------------------
MapFactory::~MapFactory()
{
  for( InfoEpetraMap::iterator ii=map_.begin(); ii!=map_.end(); ++ii ){
    delete ii->second;
  }
}
//--------------------------------------------------------------------
MapFactory&
MapFactory::self()
{
  static MapFactory s;
  return s;
}
//--------------------------------------------------------------------
const Epetra_LocalMap &
MapFactory::get_map( const int npts )
{
  InfoEpetraMap::const_iterator ii = map_.find( npts );
  if( ii == map_.end() ){
    Epetra_LocalMap * myMap = new Epetra_LocalMap( npts, 0, *com_ );
    map_.insert( std::make_pair(npts,myMap) );
    return *myMap;
  }
  return *(ii->second);
}
//--------------------------------------------------------------------


//====================================================================


} // namespace SpatialOps
