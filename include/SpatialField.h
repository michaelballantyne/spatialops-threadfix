#ifndef UT_SpatialField_h
#define UT_SpatialField_h

#include <vector>

#include <SpatialOpsDefs.h>
#include <LinearSystem.h>


namespace SpatialOps{


  enum StorageMode{ InternalStorage, ExternalStorage };


  //==================================================================


  /**
   *  @class SpatialField
   *  @author James C. Sutherland
   *  @date   December, 2006
   *
   *  Base class for SpatialFields defined on a logically rectangular
   *  domain.
   *
   *  Policies & traits required by this class:
   *
   *    VecOps - a policy dictating the treatment of this SpatialField
   *    by the linear algebra package.
   *
   *    FieldTraits - information about this spatial field.
   *
   *     \li StorageLocation - a trait specifying the field location
   *    type (e.g. node, cell, face, etc.)
   *
   *     \li GhostTraits - information about ghosting on each logical
   *     face.
   *
   *     
   */
  template< typename VecOps,
	    typename FieldLocation,
	    typename GhostTraits >
  class SpatialField
  {
    typedef typename VecOps::VecType VecType;

  public:

    /**
     *  Construct a SpatialField.
     *
     *  @param fieldDims: The number of points (excluding any ghost
     *  cells) for the domain in each of the three ordinal directions.
     *
     *  @param fieldValues : Pointer to the field values.  Behavior is
     *  dictated by the choice of <code>StorageMode</code>
     *
     *  @param mode : Storage options.  If InternalStorage then the
     *  fieldValues will be copied into an internal buffer.  If
     *  ExternalStorage then the fieldValues will be stored externally.
     *  Efficiency suggests that ExternalStorage is best, since it will
     *  avoid excessive copies.
     */
    SpatialField( const std::vector<int> & fieldDims,
		  double * const fieldValues,
		  const StorageMode mode );

    virtual ~SpatialField();


    /**
     *  @param npts : number of points (including ghost cells)
     *  @param values : array of values to overwrite with.
     */
    inline void reset_values( const int npts,
			      const double* const values );


    //@{  /** Operators for SpatialField objects */

    inline SpatialField& operator  =(const SpatialField&);
    inline SpatialField& operator +=(const SpatialField&);
    inline SpatialField& operator -=(const SpatialField&);
    inline SpatialField& operator *=(const SpatialField&);
    inline SpatialField& operator /=(const SpatialField&);

    inline SpatialField& operator  =(const double);
    inline SpatialField& operator +=(const double);
    inline SpatialField& operator -=(const double);
    inline SpatialField& operator *=(const double);
    inline SpatialField& operator /=(const double);

    inline SpatialField& operator =(const RHS&);
    inline SpatialField& operator+=(const RHS&);
    inline SpatialField& operator-=(const RHS&);

    //}@


    VecType & get_linalg_vec(){ return vec_; }
    const VecType & get_linalg_vec() const{ return vec_; }

    /** get the total number of points (including ghost layers) */
    inline int get_ntotal() const{ return npts_; }

    template< typename Dir, typename SideType>
    static int nghost(){ return GhostTraits::template get<Dir,SideType>(); }

    inline const std::vector<int>& get_extent() const{return extent_;}


    inline double& operator[](const int i){return fieldValues_[i];}
    inline const double& operator[](const int i) const{return fieldValues_[i];}

    /** obtain a pointer to the underlying field - this should be used carefully! */
    inline       double* get_ptr()      { return fieldValues_; }
    inline const double* get_ptr() const{return fieldValues_; }

    void Print( std::ostream& ) const;

  protected:

    bool consistency_check( const SpatialField& s ) const{ return ( npts_ == s.npts_ ); }

  private:

    static int get_npts( const std::vector<int> & extent );

    VecOps linAlg_;
    const std::vector<int> extent_;
    const int npts_;
    const StorageMode storageMode_;
    double * const fieldValues_;
    VecType & vec_;

    SpatialField( const SpatialField& );
    SpatialField();

  };


  //====================================================================








  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







  //==================================================================


  template< class VecOps,
	    typename FieldLocation,
	    typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  SpatialField( const std::vector<int> & fieldDims,
		double * const fieldValues,
		const StorageMode mode )
    : extent_( fieldDims ),
      npts_( get_npts(fieldDims) ),
      storageMode_( mode ),
      
      fieldValues_( (storageMode_==ExternalStorage)
		    ? fieldValues
		    : new double[ npts_ ] ),

      vec_( linAlg_.setup_vector( npts_, fieldValues_ ) )
  {
    if( mode==InternalStorage )  reset_values( npts_, fieldValues );
  }

  template< class VecOps,
	    typename FieldLocation,
	    typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  ~SpatialField()
  {
    if( storageMode_ == InternalStorage )  delete [] fieldValues_;
    linAlg_.destroy_vector();
  }

  template< class VecOps,
	    typename FieldLocation,
	    typename GhostTraits >
  void
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  reset_values( const int npts,
		     const double* const values )
  {
    assert( npts == npts_ );
    if( NULL == values )
      for( int i=0; i<npts; ++i ) fieldValues_[i] = 0.0;
    else
      for( int i=0; i<npts; ++i ) fieldValues_[i] = values[i];
  }

  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] = s.fieldValues_[i];
    return *this;
  }

  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator+=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] += s.fieldValues_[i];
    return *this;
  }

  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator-=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] -= s.fieldValues_[i];
    return *this;
  }
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator *=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] *= s.fieldValues_[i];
    return *this;
  }
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator /=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] /= s.fieldValues_[i];
    return *this;
  }
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator =(const double a){
    for( int i=0; i<npts_; ++i ) fieldValues_[i] = a;
    return *this;
  }
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator +=(const double a)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] += a;
    return *this;
  }
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator -=(const double a)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] -= a;
    return *this;
  }
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator *=(const double a)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] *= a;
    return *this;
  }

  template< class VecOps, typename FieldLocation, typename GhostTraits >
  int
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  get_npts( const std::vector<int> & extent )
  {
    int npts = 1;
    npts *= extent[0]
      + GhostTraits::template get<XDIR,SideMinus>()
      + GhostTraits::template get<XDIR,SidePlus >();
    if( extent[1]>1 ){
      npts *= extent[1]
	+ GhostTraits::template get<YDIR,SideMinus>()
	+ GhostTraits::template get<YDIR,SidePlus >();
    }
    if( extent[2]>1 ){
      npts *= extent[2]
	+ GhostTraits::template get<ZDIR,SideMinus>()
	+ GhostTraits::template get<ZDIR,SidePlus >();
    }
    return npts;
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator =(const RHS& rhs)
  {
    double * const f = get_ptr();

    const std::vector<double> & r = rhs.get_field();

    // get the dimensions of the field
    const int nxf = extent_[0]
      + GhostTraits::template get<XDIR,SideMinus>()
      + GhostTraits::template get<XDIR,SidePlus>();

    const int nyf= ( extent_[1] > 1 ) ?
      extent_[1]
      + GhostTraits::template get<YDIR,SideMinus>()
      + GhostTraits::template get<YDIR,SidePlus>()
      : 1;
//     const int nzf= ( extent_[2] > 1 ) ?
//       extent_[2]
//       + GhostTraits::template get<YDIR,SideMinus>()
//       + GhostTraits::template get<YDIR,SidePlus>()
//       : 1;

    // get the dimensions of the rhs
    const int nxr = rhs.get_extent()[0];
    const int nyr = rhs.get_extent()[1];
    const int nzr = rhs.get_extent()[2];

    //
    // RHS fields do not have ghosting.
    // Thus, we must be very careful on the indices!
    //
    static const int ngxm = GhostTraits::template get<XDIR,SideMinus>();
    static const int ngxp = GhostTraits::template get<XDIR,SidePlus>();
    static const int ngym = GhostTraits::template get<YDIR,SideMinus>();
    static const int ngyp = GhostTraits::template get<YDIR,SidePlus>();

    int ixr = 0;
    int ixf = GhostTraits::template get<XDIR,SideMinus>();
    if( extent_[1] > 1 )  ixf += GhostTraits::template get<YDIR,SideMinus>();
    if( extent_[2] > 1 )  ixf += ( nxf * GhostTraits::template get<YDIR,SideMinus>() )
			    * ( nyf * GhostTraits::template get<ZDIR,SideMinus>() );

    for( int k=0; k<nzr; ++k ){
      for( int j=0; j<nyr; ++j ){
	for( int i=0; i<nxr; ++i ){
	  f[ixf] = r[ixr];
	  ++ixf;
	  ++ixr;
	}
	ixf += ngxm + ngxp;
      }
      ixf += nxf * (ngym+ngyp);
    }
    return *this;
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator +=(const RHS& rhs)
  {
    double * const f = get_ptr();

    const double * const r = rhs.get_ptr();

    // get the dimensions of the field
    const int nxf = extent_[0]
      + GhostTraits::template get<XDIR,SideMinus>()
      + GhostTraits::template get<XDIR,SidePlus>();

    const int nyf= ( extent_[1] > 1 ) ?
      extent_[1]
      + GhostTraits::template get<YDIR,SideMinus>()
      + GhostTraits::template get<YDIR,SidePlus>()
      : 1;

//     const int nzf= ( extent_[2] > 1 ) ?
//       extent_[2]
//       + GhostTraits::template get<YDIR,SideMinus>()
//       + GhostTraits::template get<YDIR,SidePlus>()
//       : 1;

    const int nxr = rhs.get_extent()[0];
    const int nyr = rhs.get_extent()[1];
    const int nzr = rhs.get_extent()[2];

    //
    // RHS fields do not have ghosting.
    // Thus, we must be very carful on the indices!
    //

    static const int ngxm = GhostTraits::template get<XDIR,SideMinus>();
    static const int ngxp = GhostTraits::template get<XDIR,SidePlus>();
    static const int ngym = GhostTraits::template get<YDIR,SideMinus>();
    static const int ngyp = GhostTraits::template get<YDIR,SidePlus>();

    int ixr = 0;
    int ixf = GhostTraits::template get<XDIR,SideMinus>();

    if( extent_[1] > 1 )   ixf += GhostTraits::template get<YDIR,SideMinus>();
    if( extent_[2] > 1 )   ixf += ( nxf * GhostTraits::template get<YDIR,SideMinus>() )
			     * ( nyf * GhostTraits::template get<ZDIR,SideMinus>() );

    for( int k=0; k<nzr; ++k ){
      for( int j=0; j<nyr; ++j ){
	for( int i=0; i<nxr; ++i ){
	  f[ixf] += r[ixr];
	  ++ixf;
	  ++ixr;
	}
	ixf += ngxm+ngxp;
      }
      ixf += nxf * ( ngym+ngyp );
    }
    return *this;
  }
  //--------------------------------------------------------------------
  /*  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField< VecOps, FieldLocation, GhostTraits >&
  SpatialField< VecOps, FieldLocation, GhostTraits >::
  operator -=(const RHS& rhs)
  {
    double * const f = get_ptr();

    const double * const r = rhs.get_ptr();

    // get the dimensions of the field
    const int nxf= (extent_[0]>1)
      ? extent_[0]
      + GhostTraits::template get<XDIR,SideMinus>()
      + GhostTraits::template get<XDIR,SidePlus>()
      : 1;
    const int nyf= (extent_[1]>1)
      ? extent_[1]
      + GhostTraits::template get<YDIR,SideMinus>()
      + GhostTraits::template get<YDIR,SidePlus>()
      : 1;
    const int nzf= (extent_[2]>1)
      ? extent_[2]
      + GhostTraits::template get<YDIR,SideMinus>()
      + GhostTraits::template get<YDIR,SidePlus>()
      : 1;

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
  */
  //==================================================================


} // namespace SpatialOps

#endif
