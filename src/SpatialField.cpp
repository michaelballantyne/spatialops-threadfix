#include <SpatialField.h>
#include <SpatialOperator.h>
#include <LinearSystem.h>

#include <numeric>

#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>

using std::vector;

namespace SpatialOps{


//====================================================================


//--------------------------------------------------------------------
SpatialField::SpatialField( const vector<int> & fieldDims,
			    const vector<int> & nghost,
			    double * const fieldValues,
			    const StorageMode mode )
  : extent_( fieldDims ),
    npts_( get_npts(fieldDims,nghost) ),

    nghost_( nghost ),

    storageMode_( mode ),

    fieldValues_( (storageMode_==ExternalStorage)
		  ? fieldValues
		  : new double[ npts_ ] )
{
  //
  // copy the values into the internal buffer if requested.
  //
  if( mode==InternalStorage )  reset_values( npts_, fieldValues );
 

 // Build the Trilinos Epetra vector.  Use the "view" option to allow
  // us to directly control the memory that Trilinos uses, rather than
  // allowing Trilinos to make copies.  In this way, modifications to
  // fieldValues_ will directly affect the values that Trilinos sees.
  const Epetra_LocalMap & epetraMap = MapFactory::self().get_map( npts_ );
  vec_ = new Epetra_Vector( View, epetraMap, fieldValues_ );
}
//--------------------------------------------------------------------
SpatialField::~SpatialField()
{
  if( storageMode_ == InternalStorage )
    delete [] fieldValues_;

  delete vec_;
}
//--------------------------------------------------------------------
int
SpatialField::get_npts( const vector<int> & extent,
			const vector<int> & nghost )
{
  int npts=1;
  vector<int>::const_iterator ig = nghost.begin();
  for( vector<int>::const_iterator ii=extent.begin(); ii!=extent.end(); ++ii ){
    if( *ii > 1 ){
      // add "left" side ghosts
      int nn = *ii;
      nn += *ig;
      ++ig;
      nn += *ig;
      npts *= nn;
    }
  }
  return npts;
}
//--------------------------------------------------------------------
void
SpatialField::reset_values( const int npts,
			    const double* const values )
{
  assert( npts == npts_ );

  if( NULL == values )
    for( int i=0; i<npts; ++i )
      fieldValues_[i] = 0.0;
  else
    for( int i=0; i<npts; ++i )
      fieldValues_[i] = values[i];
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator =(const SpatialField& s)
{
  assert( consistency_check(s) );
  for( int i=0; i<npts_; ++i )
    fieldValues_[i] = s.fieldValues_[i];
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator +=(const SpatialField& s)
{
  assert( consistency_check( s ) );
  for( int i=0; i<npts_; ++i )
    fieldValues_[i] += s.fieldValues_[i];
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator -=(const SpatialField& s)
{
  assert( consistency_check(s) );
  for( int i=0; i<npts_; ++i )
    fieldValues_[i] -= s.fieldValues_[i];
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator *=(const SpatialField& s)
{
  assert( consistency_check(s) );
  for( int i=0; i<npts_; ++i )
    fieldValues_[i] *= s.fieldValues_[i];
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator /=(const SpatialField& s)
{
  assert( consistency_check(s) );
  for( int i=0; i<npts_; ++i )
    fieldValues_[i] /= s.fieldValues_[i];
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator +=(const double s)
{
  for( int i=0; i<npts_; ++i ) fieldValues_[i] += s;
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator *=(const double s)
{
  for( int i=0; i<npts_; ++i ) fieldValues_[i] *= s;
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator =(const double s)
{
  for( int i=0; i<npts_; ++i ) fieldValues_[i] = s;
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator =(const RHS& rhs)
{
  double * const f = get_ptr();

  const std::vector<double> & r = rhs.get_field();

  // get the dimensions of the field
  const int nxf= (extent_[0]>1) ? extent_[0] + nghost_[0] + nghost_[1] : 1;
  const int nyf= (extent_[1]>1) ? extent_[1] + nghost_[2] + nghost_[3] : 1;
  const int nzf= (extent_[2]>1) ? extent_[2] + nghost_[4] + nghost_[5] : 1;

  // get the dimensions of the rhs
  const int nxr = rhs.get_extent()[0];
  const int nyr = rhs.get_extent()[1];
  const int nzr = rhs.get_extent()[2];

  //
  // RHS fields do not have ghosting.
  // Thus, we must be very careful on the indices!
  //

  int ixr = 0;
  int ixf = nghost_[0];
  if( nyf>1 ) ixf += nxf;
  if( nzf>1 ) ixf += nxf*nyf;
  for( int k=0; k<nzr; ++k ){
    for( int j=0; j<nyr; ++j ){
      for( int i=0; i<nxr; ++i ){
	f[ixf] = r[ixr];
	++ixf;
	++ixr;
      }
      ixf += nghost_[0] + nghost_[1];
    }
    ixf += nxf * ( nghost_[2]+nghost_[3] );
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator +=(const RHS& rhs)
{
  double * const f = get_ptr();

  const double * const r = rhs.get_ptr();

  // get the dimensions of the field
  const int nxf= (extent_[0]>1) ? extent_[0] + nghost_[0] + nghost_[1] : 1;
  const int nyf= (extent_[1]>1) ? extent_[1] + nghost_[2] + nghost_[3] : 1;
  const int nzf= (extent_[2]>1) ? extent_[2] + nghost_[4] + nghost_[5] : 1;

  const int nxr = rhs.get_extent()[0];
  const int nyr = rhs.get_extent()[1];
  const int nzr = rhs.get_extent()[2];

  //
  // RHS fields do not have ghosting.
  // Thus, we must be very carful on the indices!
  //

  int ixr = 0;
  int ixf = nghost_[0];
  if( nyf>1 ) ixf += nxf;
  if( nzf>1 ) ixf += nxf*nyf;
  for( int k=0; k<nzr; ++k ){
    for( int j=0; j<nyr; ++j ){
      for( int i=0; i<nxr; ++i ){
	f[ixf] += r[ixr];
	++ixf;
	++ixr;
      }
      ixf += nghost_[0]+nghost_[1];
    }
    ixf += nxf * ( nghost_[2] + nghost_[3] );
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator -=(const RHS& rhs)
{
  double * const f = get_ptr();

  const double * const r = rhs.get_ptr();

  // get the dimensions of the field
  const int nxf= (extent_[0]>1) ? extent_[0] + nghost_[0] + nghost_[1] : 1;
  const int nyf= (extent_[1]>1) ? extent_[1] + nghost_[2] + nghost_[3] : 1;
  const int nzf= (extent_[2]>1) ? extent_[2] + nghost_[4] + nghost_[5] : 1;

  const int nxr = rhs.get_extent()[0];
  const int nyr = rhs.get_extent()[1];
  const int nzr = rhs.get_extent()[2];

  //
  // RHS fields do not have ghosting.
  // Thus, we must be very carful on the indices!
  //

  int ixr = 0;
  int ixf = nghost_[0];
  if( nyf>1 ) ixf += nxf;
  if( nzf>1 ) ixf += nxf*nyf;
  for( int k=0; k<nzr; ++k ){
    for( int j=0; j<nyr; ++j ){
      for( int i=0; i<nxr; ++i ){
	f[ixf] -= r[ixr];
	++ixf;
	++ixr;
      }
      ixf += nghost_[0]+nghost_[1];
    }
    ixf += nxf * ( nghost_[2] + nghost_[3] );
  }
  return *this;
}
//--------------------------------------------------------------------
Epetra_Vector&
SpatialField::epetra_vec()
{
  return *vec_;
}
//--------------------------------------------------------------------
const Epetra_Vector&
SpatialField::epetra_vec() const
{
  return *vec_;
}
//--------------------------------------------------------------------
void
SpatialField::Print( std::ostream& c ) const
{
  epetra_vec().Print(c);
}
//--------------------------------------------------------------------
bool
SpatialField::consistency_check( const SpatialField& s ) const
{
  return ( npts_ == s.npts_ );
}
//--------------------------------------------------------------------

} // namespace SpatialOps
