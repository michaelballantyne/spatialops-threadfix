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
SpatialField::SpatialField( const std::vector<int> & fieldDims,
			    const int nghost,
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
			const int nghost )
{
  int npts = 1;
  for( vector<int>::const_iterator ii=extent.begin(); ii!=extent.end(); ++ii ){
    if( *ii>1 ) npts *= *ii+2*nghost;
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
  const int ng = nghost();
  double * const f = get_ptr();

  const std::vector<double> & r = rhs.get_field();

  const int nxf= (extent_[0]>1) ? extent_[0]+2*ng : 1;
  const int nyf= (extent_[1]>1) ? extent_[1]+2*ng : 1;

  const int nxr = rhs.get_extent()[0];
  const int nyr = rhs.get_extent()[1];
  const int nzr = rhs.get_extent()[2];

  for( int k=0; k<nzr; ++k ){
    int ixf = k*nxf*nyf + ng;
    int ixr = k*nxr*nyr;
    for( int j=0; j<nyr; ++j ){
      ixf += j*nxf + ng;
      ixr += j*nxr;
      for( int i=0; i<nxr; ++i ){
	ixf += i+ng;
	ixr += i;
	f[ixf] = r[ixr];
      }
    }
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator +=(const RHS& rhs)
{
  const int ng = nghost();
  double * const f = get_ptr();

  const double * const r = rhs.get_ptr();

  const int nxf= (extent_[0]>1) ? extent_[0]+2*ng : 1;
  const int nyf= (extent_[1]>1) ? extent_[1]+2*ng : 1;

  const int nxr = rhs.get_extent()[0];
  const int nyr = rhs.get_extent()[1];
  const int nzr = rhs.get_extent()[2];

  for( int k=0; k<nzr; ++k ){
    int ixf = k*nxf*nyf + ng;
    int ixr = k*nxr*nyr;
    for( int j=0; j<nyr; ++j ){
      ixf += j*nxf + ng;
      ixr += j*nxr;
      for( int i=0; i<nxr; ++i ){
	ixf += i+ng;
	ixr += i;
	f[ixf] += r[ixr];
      }
    }
  }
  return *this;
}
//--------------------------------------------------------------------
SpatialField&
SpatialField::operator -=(const RHS& rhs)
{
  const int ng = nghost();
  double * const f = get_ptr();

  const double * const r = rhs.get_ptr();

  const int nxf= (extent_[0]>1) ? extent_[0]+2*ng : 1;
  const int nyf= (extent_[1]>1) ? extent_[1]+2*ng : 1;

  const int nxr = rhs.get_extent()[0];
  const int nyr = rhs.get_extent()[1];
  const int nzr = rhs.get_extent()[2];

  for( int k=0; k<nzr; ++k ){
    int ixf = k*nxf*nyf + ng;
    int ixr = k*nxr*nyr;
    for( int j=0; j<nyr; ++j ){
      ixf += j*nxf + ng;
      ixr += j*nxr;
      for( int i=0; i<nxr; ++i ){
	ixf += i+ng;
	ixr += i;
	f[ixf] -= r[ixr];
      }
    }
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
bool
SpatialField::consistency_check( const SpatialField& s ) const
{
  return ( npts_ == s.npts_ );
}
//--------------------------------------------------------------------

} // namespace SpatialOps
