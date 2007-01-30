#include <SpatialOperator.h>
#include <SpatialField.h>

//-----------------------------
#include <numeric>
#include <iostream>
#include <sstream>
#include <stdexcept>
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
				  const std::vector<int> & nghost,
				  const int entriesPerRow,
				  const std::vector<int> & extent )
  : extent_( extent ),
    nrows_ ( nrows  ),
    ncols_ ( ncols  ),

    entriesPerRow_( entriesPerRow ),

    nghost_( nghost ),

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

  const bool useTranspose = false;
  const int flag =
    EpetraExt::MatrixMatrix::Multiply( epetra_mat(), useTranspose,
				       B.epetra_mat(), useTranspose,
				       C.epetra_mat() );
  if( flag!=0 )
    std::cout << std::endl
	      << "ERROR!  Flag=" << flag
	      << " returned from EpetraExt::MatrixMatrix::Multiply()." << std::endl
	      << "        This likely indicates incompatible matrices for multiplication." << std::endl
	      << "        Check matrix sparsity patterns and dimensions for compatibility."
	      << std::endl << std::endl;
  //  assert( flag==0 );
}
//--------------------------------------------------------------------
SpatialOperator&
SpatialOperator::operator = ( const SpatialOperator& op )
{
  assert( op.nrows() == nrows() );
  assert( op.ncols() == ncols() );

  // jcs this is INEFFICIENT.  It would be better to do this directly.
  reset_entries(0.0);
  return (*this)+=(op);
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
  /*
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
  */
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

  /*
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
  */
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
  const int flag = epetra_mat().LeftScale( f.epetra_vec() );
  assert( flag == 0 );
}
//--------------------------------------------------------------------
void
SpatialOperator::right_scale( const SpatialField& f )
{
  const int flag = epetra_mat().RightScale( f.epetra_vec() );
  assert( flag == 0 );
}
//--------------------------------------------------------------------
void
SpatialOperator::reset_entries( const double val )
{
  mat_->PutScalar(val);
  return;
}
//--------------------------------------------------------------------
void
SpatialOperator::insert_row_entry( const int rownum,
				   std::vector<double> & rowValues,
				   std::vector<int> & rowIndices )
{
  const int flag = mat_->InsertMyValues( rownum,
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
void
SpatialOperator::Print( std::ostream & c ) const
{
  epetra_mat().Print(c);
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


//--------------------------------------------------------------------
SpatialOpDatabase&
SpatialOpDatabase::self()
{
  static SpatialOpDatabase s;
  return s;
}
//--------------------------------------------------------------------
void
SpatialOpDatabase::register_new_operator( const OperatorType opType,
					 SpatialOperator * const op,
					 const std::string & name,
					 const bool makeDefault )
{
  std::pair< NameOpMap::const_iterator, bool > result =
    allOps_.insert( make_pair(name,op) );

  if( !result.second ){
    std::ostringstream msg;
    msg << "ERROR! Could not insert a spatial operator of type " << type2name( opType ) << std::endl
	<< "       and name '" << name << "' - this must be a duplicate entry!" << std::endl;
    throw std::runtime_error( msg.str() );
  }

  if( makeDefault || (activeOps_.find(opType)==activeOps_.end()) ){
    activeOps_[opType] = op;
  }
}
//--------------------------------------------------------------------
SpatialOperator*&
SpatialOpDatabase::retrieve_operator( const std::string & name )
{
  NameOpMap::iterator ii = allOps_.find( name );
  if( ii == allOps_.end() ){
    std::ostringstream msg;
    msg << "ERROR! No operator named '" << name << "' has been registered!" << std::endl;
    throw std::runtime_error( msg.str() );
  }
  return ii->second;
}
//--------------------------------------------------------------------
SpatialOperator*&
SpatialOpDatabase::retrieve_operator( const OperatorType opType )
{
  TypeOpMap::iterator ii = activeOps_.find( opType );
  if( ii == activeOps_.end() ){
    std::ostringstream msg;
    msg << "ERROR! No operator of type '" << type2name(opType) << "' has been registered!" << std::endl;
    throw std::runtime_error( msg.str() );
  }
  return ii->second;
}
//--------------------------------------------------------------------
void
SpatialOpDatabase::set_default_operator( const OperatorType opType,
					const std::string & opName )
{
  // do we have an operator with the given name?  If not, just return and do nothing.
  NameOpMap::iterator ii = allOps_.find( opName );
  if( ii == allOps_.end() ) return;

  // do we have an operator of this type?  If not, just return and do nothing.
  TypeOpMap::iterator jj = activeOps_.find( opType );
  if( jj != activeOps_.end() ) jj->second = ii->second;
}
//--------------------------------------------------------------------
const std::string
SpatialOpDatabase::type2name( const OperatorType opType ) const
{
  const static std::string DivName     = "Divergence";
  const static std::string GradName    = "Gradient";
  const static std::string InterpName  = "Interpolant";
  const static std::string ScratchName = "Scratch";

  switch( opType ){

  case CELL_DIVERGENCE_X:  return (DivName+"-X Cell");
  case FACE_DIVERGENCE_X:  return (DivName+"-X Face");

  case CELL_DIVERGENCE_Y:  return (DivName+"-Y Cell");
  case FACE_DIVERGENCE_Y:  return (DivName+"-Y Face");

  case CELL_DIVERGENCE_Z:  return (DivName+"-Z Cell");
  case FACE_DIVERGENCE_Z:  return (DivName+"-Z Face");


  case CELL_GRADIENT_X:    return (GradName+"-X Cell");
  case FACE_GRADIENT_X:    return (GradName+"-X Face");

  case CELL_GRADIENT_Y:    return (GradName+"-Y Cell");
  case FACE_GRADIENT_Y:    return (GradName+"-Y Face");

  case CELL_GRADIENT_Z:    return (GradName+"-Z Cell");
  case FACE_GRADIENT_Z:    return (GradName+"-Z Face");


  case CELL_INTERPOLANT_X: return (InterpName+"-X Cell");
  case FACE_INTERPOLANT_X: return (InterpName+"-X Face");

  case CELL_INTERPOLANT_Y: return (InterpName+"-Y Cell");
  case FACE_INTERPOLANT_Y: return (InterpName+"-Y Face");

  case CELL_INTERPOLANT_Z: return (InterpName+"-Z Cell");
  case FACE_INTERPOLANT_Z: return (InterpName+"-Z Face");


  case CELL_SCRATCH_X: return (ScratchName+"-X Cell");
  case FACE_SCRATCH_X: return (ScratchName+"-X Face");

  case CELL_SCRATCH_Y: return (ScratchName+"-Y Cell");
  case FACE_SCRATCH_Y: return (ScratchName+"-Y Face");

  case CELL_SCRATCH_Z: return (ScratchName+"-Z Cell");
  case FACE_SCRATCH_Z: return (ScratchName+"-Z Face");

  default:
    throw std::runtime_error( "ERROR!  Invalid OperatorType in SpatialOpDatabase::type2name()\n" );

  }
  const static std::string err = "ERROR";
  return err;
}
//--------------------------------------------------------------------
SpatialOpDatabase::SpatialOpDatabase()
{
}
//--------------------------------------------------------------------
SpatialOpDatabase::~SpatialOpDatabase()
{
  for( NameOpMap::iterator ii=allOps_.begin(); ii!=allOps_.end(); ++ii ){
    delete ii->second;
  }
}
//--------------------------------------------------------------------


//====================================================================


} // namespace SpatialOps
