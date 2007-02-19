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
				  const std::vector<int> & nghostSrc,
				  const std::vector<int> & nghostDest,
				  const int entriesPerRow,
				  const std::vector<int> & extent )
  : extent_( extent ),
    nrows_ ( nrows  ),
    ncols_ ( ncols  ),

    entriesPerRow_( entriesPerRow ),

    nghostSrc_ ( nghostSrc  ),
    nghostDest_( nghostDest ),

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

  assert( compatibility_check(B,false) );
  assert( compatibility_check(C,true ) );

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
bool
SpatialOperator::is_row_ghost( const int irow,
			       IndexTriplet* const ix ) const
{
  const bool getIxs = (ix != NULL);
  if( getIxs ){
    IndexTriplet & i= *ix;
    i[0] = -1;
    i[1] = -1;
    i[2] = -1;
  }

  bool isGhost = false;

  // is this row a ghost entry?  Row corresponds to entry in dest vec.
  const int nxd = extent_[0] + nghostDest_[0] + nghostDest_[1];
  const int idest = irow % nxd - nghostDest_[0];
  if( idest < 0  ||  idest >= extent_[0] ){
    isGhost = true;
    if( !getIxs ) return true;
  }
  if( getIxs ) (*ix)[0] = idest;

  const int nyd = extent_[1] + nghostDest_[2] + nghostDest_[3];
  if( extent_[1] > 1 ){
    const int jdest = irow/nxd % nyd - nghostDest_[2];
    if( jdest < 0  ||  jdest >= extent_[1] ){
      isGhost = true;
      if( !getIxs ) return true;
    }
    if( getIxs ) (*ix)[1] = jdest;
  }

  if( extent_[2] > 1 ){
    const int kdest = irow/(nxd*nyd) - nghostDest_[4];
    if( kdest < 0  ||  kdest >= extent_[2] ){
      isGhost = true;
    }
    if( getIxs ) (*ix)[2] = kdest;
  }
  return isGhost;
}
//--------------------------------------------------------------------
bool
SpatialOperator::is_col_ghost( const int icol,
			       IndexTriplet* const ix ) const
{
  const bool getIxs = (ix!=NULL);
  if( getIxs ){
    IndexTriplet & i= *ix;
    i[0] = -1;
    i[1] = -1;
    i[2] = -1;
  }

  bool isGhost = false;

  // is this column a ghost entry?  Column corresponds to entry in src vec.
  const int nxs = extent_[0] + nghostSrc_[0] + nghostSrc_[1];
  const int isrc = icol%nxs - nghostSrc_[0];
  if( isrc < 0  ||  isrc >= extent_[0] ){
    isGhost = true;
    if( !getIxs ) return true;
  }
  if( getIxs ) (*ix)[0] = isrc;

  const int nys = extent_[1] + nghostSrc_[2] + nghostSrc_[3];
  if( extent_[1] > 1 ){
    const int jsrc = (icol/nxs) % nys - nghostSrc_[2];
    if( jsrc < 0  ||  jsrc >= extent_[1] ){
      isGhost = true;
      if( !getIxs ) return true;
    }
    if( getIxs ) (*ix)[1] = jsrc;
  }
  if( extent_[2] > 1 ){
    const int ksrc = icol/(nxs*nys) - nghostSrc_[4];
    if( ksrc < 0  ||  ksrc >= extent_[2] )   isGhost = true;
    if( getIxs ) (*ix)[2] = ksrc;
  }
  return isGhost;
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
SpatialOperator::compatibility_check( const SpatialOperator& op,
				      const bool isResultOp ) const
{
  if( isResultOp ){
    if( nrows_ != op.nrows_ ){
      std::cout << "Destination matrix must have same number of rows as operator." << std::endl;
      return false;
    }
  }
  else{
    if( ncols_ != op.nrows_ ){
      std::cout << "Matrix dimensions are incompatible for multiplication." << std::endl;
      return false;
    }
  }
  return true;

  // ensure identical ghosting patterns
  const vector<int> & opgs = op.nghost_src();
  const vector<int> & opgd = op.nghost_dest();
  vector<int>::const_iterator iopgs = opgs.begin();
  vector<int>::const_iterator iopgd = opgd.begin();

  vector<int>::const_iterator igs = nghost_src().begin();
  vector<int>::const_iterator igd = nghost_dest().begin();

  for( ; igs != nghost_src().end();  ++igs, ++iopgs ){
    if( *iopgs != *igs ){
      std::cout << "ghosting incompatibility detected in operators.  Source field ghosting is not identical." << std::endl;
      return false;
    }
  }

  for( ; igd != nghost_dest().end();  ++igd, ++iopgd ){
    if( *iopgd != *igd ){
      std::cout << "ghosting incompatibility detected in operators. Destination field ghosting is not identical." << std::endl;
      return false;
    }
  }

  return true;
}
//--------------------------------------------------------------------
bool
SpatialOperator::compatibility_check( const SpatialField & field,
				      const FieldType fldType ) const
{
  switch( fldType ){

  case SOURCE_FIELD:{
    if( ncols_ != field.epetra_vec().GlobalLength() ){
      std::cout << "expecting " << ncols_ << " entries for source field, found "
		<< field.epetra_vec().GlobalLength() << std::endl;
      return false;
    }

    // verify ghosting compatibility
    const vector<int> & fg = field.nghost();
    vector<int>::const_iterator ifg = fg.begin();
    for( vector<int>::const_iterator ig = nghost_src().begin();
	 ig!=nghost_src().end();
	 ++ig, ++ifg )
      {
	if( *ifg != *ifg ){
	  std::cout << "ghost incompatibility detected between field and spatial operator." << std::endl;
	  return false;
	}
      }

    break;
  }

  case DEST_FIELD:{
    if( nrows_ != field.epetra_vec().GlobalLength() ){
      std::cout << "expecting " << nrows_ << " entries for destination field, found "
		<< field.epetra_vec().GlobalLength() << std::endl;
      return false;
    }

    const vector<int> & fg = field.nghost();
    vector<int>::const_iterator ifg = fg.begin();
    for( vector<int>::const_iterator ig = nghost_dest().begin();
	 ig!=nghost_dest().end();
	 ++ig, ++ifg )
      {
	if( *ifg != *ifg ){
	  std::cout << "ghost incompatibility detected between field and spatial operator." << std::endl;
	  return false;
	}
      }

    break;
  }

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
SpatialOpDatabase::OperatorDescriptor::OperatorDescriptor( const Direction    d,
							   const OperatorLoc  l,
							   const OperatorType t )
  : dir ( d ),
    loc ( l ),
    type( t )
{}
//--------------------------------------------------------------------
bool
SpatialOpDatabase::OperatorDescriptor::operator < ( const OperatorDescriptor& s ) const
{
  // jcs: this is REALLY SLOW!  I did this because I couldn't get the
  // straight enum compare to work properly for some reason.
  ostringstream s1, s2;
  s1 << dir << " " << loc << " " << type;
  s2 << s.dir << " "  << s.loc << " " << s.type;
  return (s1.str()<s2.str());
  /*
  bool isLess =
    (dir < s.dir ) &&
    (loc < s.loc ) &&
    (type< s.type);

  return isLess;
  */
}
//--------------------------------------------------------------------
bool
SpatialOpDatabase::OperatorDescriptor::operator== ( const OperatorDescriptor& s ) const
{
  return ( (dir == s.dir)  &&
	   (loc == s.loc)  &&
	   (type== s.type) );
}

//--------------------------------------------------------------------
SpatialOpDatabase&
SpatialOpDatabase::self()
{
  static SpatialOpDatabase s;
  return s;
}
//--------------------------------------------------------------------
void
SpatialOpDatabase::register_new_operator( const OperatorDescriptor & opType,
					 SpatialOperator * const op,
					 const std::string & name,
					 const bool makeDefault )
{
  Shape s( op->get_extent(), op->nghost_src(), op->nghost_dest() );

  TypeShapeMap::iterator itsm = typeMap_.find( opType );
  if( itsm == typeMap_.end() ){
    ShapeOpMap som;
    som.insert( make_pair(s,op) );
    typeMap_.insert( make_pair(opType,som) );
  }
  else{
    ShapeOpMap & som = itsm->second;

    ShapeOpMap::iterator iop = som.find( s );
    if( iop == som.end() ){
      pair<ShapeOpMap::const_iterator,bool> result = som.insert( make_pair( s, op ) );
    }
    else if( makeDefault ){
      iop->second = op;
    }
  }

  NameShapeMap::iterator insm = nameMap_.find( name );
  if( insm == nameMap_.end() ){
    ShapeOpMap som;
    som.insert( make_pair(s,op) );
    nameMap_[name] = som;
  }
  else{
    ShapeOpMap & som = insm->second;
    pair<ShapeOpMap::const_iterator,bool> result = som.insert( make_pair( s, op ) );
    if( !result.second ){
      std::ostringstream msg;
      msg << "ERROR!  Operator named '" << name << "' has already been registered!" << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }

}
//--------------------------------------------------------------------
SpatialOperator*&
SpatialOpDatabase::retrieve_operator( const std::string & name,
				      const std::vector<int> & nxyz,
				      const std::vector<int> & nghostSrc,
				      const std::vector<int> & nghostDest )
{
  NameShapeMap::iterator ii = nameMap_.find( name );
  if( ii == nameMap_.end() ){
    std::ostringstream msg;
    msg << "ERROR!  No operator named '" << name << "'" << std::endl
	<< "        has been registered!" << std::endl;
    throw std::runtime_error( msg.str() );
  }

  ShapeOpMap & som = ii->second;
  Shape s( nxyz, nghostSrc, nghostDest );
  ShapeOpMap::iterator iop = som.find( s );
  if( iop == som.end() ){
    std::ostringstream msg;
    msg << "ERROR!  No operator named '" << name << "'" << std::endl
	<< "        with the requested dimensions and ghosting has been registered!" << std::endl;
    throw std::runtime_error( msg.str() );
  }
  return iop->second;
}
//--------------------------------------------------------------------
SpatialOperator*&
SpatialOpDatabase::retrieve_operator( const OperatorDescriptor & opType,
				      const std::vector<int> & nxyz,
				      const std::vector<int> & nghostSrc,
				      const std::vector<int> & nghostDest )
{
  // look for an operator of this type
  TypeShapeMap::iterator itsm = typeMap_.find( opType );

  if( itsm == typeMap_.end() ){
    std::ostringstream msg;
    msg << "ERROR!  Attempted to retrieve an operator that does not exist." << std::endl
	<< "        Check the operator type and shape (nx,ny,nz and nghost)" << std::endl
	<< "        and ensure that an operator of this type has been registered" << std::endl
	<< "        in the database." << std::endl;
    throw std::runtime_error( msg.str() );
  }

  ShapeOpMap & som = itsm->second;
  Shape s( nxyz, nghostSrc, nghostDest );

  // look for one with this shape
  ShapeOpMap::iterator iop = som.find( s );
  if( iop == som.end() ){
    std::ostringstream msg;
    msg << "ERROR!  Attempted to retrieve an operator that does not exist." << std::endl
	<< "        Check the operator shape (nx,ny,nz and nghost) and ensure" << std::endl
	<< "        that an operator of this type and shape has been registered" << std::endl
	<< "        in the database." << std::endl;
    throw std::runtime_error( msg.str() );
  }

  return iop->second;
}
//--------------------------------------------------------------------
void
SpatialOpDatabase::set_default_operator( const OperatorDescriptor & opType,
					 const std::string & opName,
					 const std::vector<int> & nxyz,
					 const std::vector<int> & nghostSrc,
					 const std::vector<int> & nghostDest )
{
  // do we have an operator with the given name?  If not, just return and do nothing.
  NameShapeMap::iterator ii = nameMap_.find( opName );
  if( ii == nameMap_.end() ) return;

  // do we have an operator of this type?  If not, just return and do nothing.
  TypeShapeMap::iterator jj = typeMap_.find( opType );
  if( jj == typeMap_.end() ) return;


  // See if we have an operator with this name AND shape
  Shape s( nxyz, nghostSrc, nghostDest );
  ShapeOpMap & somName = ii->second;
  ShapeOpMap::iterator iopn = somName.find( s );

  if( iopn != somName.end() ){
    ShapeOpMap & somType = jj->second;
    ShapeOpMap::iterator iopt = somType.find(s);
    iopt->second = iopn->second;
  }
}
//--------------------------------------------------------------------
SpatialOpDatabase::SpatialOpDatabase()
{
}
//--------------------------------------------------------------------
SpatialOpDatabase::~SpatialOpDatabase()
{
  for( NameShapeMap::iterator ii=nameMap_.begin(); ii!=nameMap_.end(); ++ii ){
    ShapeOpMap & som = ii->second;
    for( ShapeOpMap::iterator jj=som.begin(); jj!=som.end(); ++jj ){
      delete jj->second;
    }
  }
}
//--------------------------------------------------------------------

//====================================================================

SpatialOpDatabase::Shape::Shape( const std::vector<int> & extent,
				 const std::vector<int> ghostSrc,
				 const std::vector<int> ghostDest )
{
  nxyz = extent;
  ngS = ghostSrc;
  ngD = ghostDest;
}
//--------------------------------------------------------------------
bool
SpatialOpDatabase::Shape::operator==( const Shape & s ) const
{
  bool isequal = true;

  if( s.nxyz.size() != nxyz.size() ) return false;

  if( s.ngS.size() != ngS.size() ) return false;
  if( s.ngD.size() != ngD.size() ) return false;

  vector<int>::const_iterator is = s.nxyz.begin();
  vector<int>::const_iterator ii = nxyz.begin();
  for( ; ii!=nxyz.end(); ++ii, ++is ){
    if( *ii != *is ) isequal = false;
  }

  is = s.ngS.begin();
  ii =   ngS.begin();
  for( ; ii!=ngS.end(); ++ii, ++is ){
    if( *ii != *is ) isequal = false;
  }

  is = s.ngD.begin();
  ii =   ngD.begin();
  for( ; ii!=ngD.end(); ++ii, ++is ){
    if( *ii != *is ) isequal = false;
  }

  return isequal;
}
//--------------------------------------------------------------------
bool
SpatialOpDatabase::Shape::operator < ( const Shape& s ) const
{
  bool isLess = true;

  if( s.nxyz.size() != nxyz.size() ) return false;
  if( s.ngS.size()  != ngS.size()  ) return false;
  if( s.ngD.size()  != ngD.size()  ) return false;

  vector<int>::const_iterator is = s.nxyz.begin();
  vector<int>::const_iterator ii = nxyz.begin();
  for( ; ii!=nxyz.end(); ++ii, ++is ){
    if( *ii >= *is ) isLess = false;
  }

  is = s.ngS.begin();
  ii =   ngS.begin();
  for( ; ii!=ngS.end(); ++ii, ++is ){
    if( *ii >= *is ) isLess = false;
  }

  is = s.ngD.begin();
  ii =   ngD.begin();
  for( ; ii!=ngD.end(); ++ii, ++is ){
    if( *ii >= *is ) isLess = false;
  }

  return isLess;
}

//====================================================================


} // namespace SpatialOps
