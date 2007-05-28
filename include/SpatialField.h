#ifndef UT_SpatialField_h
#define UT_SpatialField_h

#include <vector>
#include <iostream>

// allows compile-time expansion of complex expressions involving SpatialField objects.
#include <daixtrose/Daixt.h>


#include <SpatialOpsDefs.h>
#include <LinearSystem.h>
#include <SpatialFieldStore.h>


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
    friend class SpatialFieldStore<VecOps,FieldLocation,GhostTraits>;
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

    inline SpatialField& operator =(const SpatialField&);
    inline SpatialField& operator+=(const SpatialField&);
    inline SpatialField& operator-=(const SpatialField&);
    inline SpatialField& operator*=(const SpatialField&);
    inline SpatialField& operator/=(const SpatialField&);

    inline SpatialField& operator =(const double);
    inline SpatialField& operator+=(const double);
    inline SpatialField& operator-=(const double);
    inline SpatialField& operator*=(const double);
    inline SpatialField& operator/=(const double);

    inline SpatialField& operator =(const RHS&);
    inline SpatialField& operator+=(const RHS&);
    inline SpatialField& operator-=(const RHS&);

    inline bool operator==(const SpatialField&);
    inline bool operator!=(const SpatialField&);

    //}@



    //@{ Binary operators

    /**
     *  note that these return references rather than copies because
     *  we have a temporary working vector that we use internally.
     *  However, this also means that it is not safe to get a
     *  reference to these variables, since they could change very
     *  easily!  That means that you SHOULD NOT do something like
     *
     *     SpatialField a,b;
     *     SpatialField & c = a+b;
     *
     *  rather, you should do:
     *
     *     SpatialField a,b,c;
     *     c = a+b;
     *
     *  This results in a copy from a temporary to c, but is safe.
     *
     *  NOTE: this could get us into big trouble if we have threads
     *  running concurrently, since we would not be able to guarantee
     *  that two threads didn't use the same memory.
     */

    inline SpatialField& operator+(const SpatialField&);
    inline SpatialField& operator-(const SpatialField&);
    inline SpatialField& operator*(const SpatialField&);
    inline SpatialField& operator/(const SpatialField&);


    // this provides support for Daixtrose - an expression template
    // engine to allow compund expressions involving SpatialField
    // objects to be unrolled by the compiler.
    template<class T>
    inline SpatialField& operator=(const Daixt::Expr<T>&E);

    //}@


    inline       VecType & get_linalg_vec()      { return vec_; }
    inline const VecType & get_linalg_vec() const{ return vec_; }


    /**
     *  Get the total number of points (including ghost layers) in this SpatialField
     */
    inline int get_ntotal() const{ return npts_; }

    /**
     *  Obtain the number of ghost cells in the given direction and side of the patch.
     */
    template< typename Dir, typename SideType>
    static int nghost(){ return GhostTraits::template get<Dir,SideType>(); }

    /**
     *  Obtain the domain extent.  This is a vector containing the
     *  number of points in each direction, excluding ghost cells.
     */
    inline const std::vector<int>& get_extent() const{return extent_;}


    /**
     *  Obtain a reference to the field using the [] operator.  This
     *  should not generally be used, as it is not tuned for
     *  performance.
     */
    inline double& operator[](const int i){return fieldValues_[i];}
    inline const double& operator[](const int i) const{return fieldValues_[i];}


    //@{  STL-compliant iterators for this field

    typedef double*           iterator;
    typedef double const*     const_iterator;

    inline iterator       begin()      {return fieldValues_;}
    inline const_iterator begin() const{return fieldValues_;}

    inline iterator       end()      {return fieldValues_+npts_;}
    inline const_iterator end() const{return fieldValues_+npts_;}

    //}@


    /** Dump information about the field to the given output stream. */
    void Print( std::ostream& ) const;

  protected:
    
    bool consistency_check( const SpatialField& s ) const{ return ( npts_ == s.npts_ ); }

  private:

    /**
     *  Given the patch extent, this returns the total number of points
     *  in this field (including ghost points)
     */
    static int get_npts( const std::vector<int> & extent );

    /**
     *  If the temporary work array assigned to this field is not
     *  valid, obtain a valid one.  This should only be called if the
     *  temporary is required, as it may result in allocation of new
     *  memory.
     */
    inline void check_tmp_validity();



    VecOps linAlg_;
    const std::vector<int> extent_;
    const int npts_;
    const StorageMode storageMode_;
    double * const fieldValues_;
    VecType & vec_;

    // these fields facilitate more efficient operations.  +,-,*,/
    // operators require a temporary.  Here we store a temporary for
    // this purpose rather than building one each time we hit one of
    // these operators.  The SpatialFieldStore is used to assist in
    // holding and building temporaries to minimize the number of them
    // that are created.  Basically, any time any of the above
    // operators are called, a temporary must be generated.  Note that
    // the increment operators like += *= etc do not require
    // temporaries.
    SpatialField* tmp_;  // a temporary field to use for efficient operations
    size_t tmpFieldNum_; // the id for the temporary field...

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


  //------------------------------------------------------------------
  template<typename T>
  void Evaluate( const Daixt::Expr<T>& arg )
  {
    Evaluate(arg.content());
  }
  //------------------------------------------------------------------
  template<typename LHS, typename RHS, typename OP>
  void Evaluate( const Daixt::BinOp<LHS,RHS,OP>& arg )
  {
    Evaluate(arg.lhs());
    cout << OP::Symbol() << endl;
    Evaluate( arg.rhs() );
  }
  //------------------------------------------------------------------
  template<typename ARG, typename OP>
  void Evaluate( const Daixt::UnOp<ARG,OP>& arg )
  {
    cout << OP::Symbol() << endl;
    Evaluate(arg.arg());
  }
  //------------------------------------------------------------------
//   void Evaluate( const SpatialOps::FVStaggeredUniform::CellFieldNoGhost& x )
//   {
//     x.Print(cout);
//     cout << endl;
//   }
  //------------------------------------------------------------------
  template<typename T>
  void Evaluate( const T& t )
  {
    t.put(cout);
    cout << endl;
  }
  //------------------------------------------------------------------


  //==================================================================


  template< class VecOps, typename FieldLocation, typename GhostTraits >
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

      vec_( linAlg_.setup_vector( npts_, fieldValues_ ) ),

      tmp_(NULL),
      tmpFieldNum_(0)
  {
    if( mode==InternalStorage )  reset_values( npts_, fieldValues );

    tmp_ = &(SpatialFieldStore<VecOps,FieldLocation,GhostTraits>::self().get( *this, tmpFieldNum_ ));
    assert( tmp_ != NULL );
  }
  //------------------------------------------------------------------
  template< class VecOps,
	    typename FieldLocation,
	    typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  ~SpatialField()
  {
    if( storageMode_ == InternalStorage )  delete [] fieldValues_;
    linAlg_.destroy_vector();
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>::SpatialField( const SpatialField<VecOps,FieldLocation,GhostTraits>& f )
    : extent_     ( f.extent_ ),
      npts_       ( f.npts_ ),
      storageMode_( f.storageMode_ ),
      fieldValues_( new double[npts_] ),
      vec_        ( linAlg_.setup_vector(npts_,fieldValues_) ),
      tmp_        (NULL),
      tmpFieldNum_(0)
  {
    // leave the tmp field NULL.  If we need it later, we will construct it then.
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  void
  SpatialField<VecOps,FieldLocation, GhostTraits>::check_tmp_validity()
  {
    if( tmp_ == NULL ){
      ++tmpFieldNum_;
      tmp_ = &SpatialFieldStore<VecOps,FieldLocation,GhostTraits>::self().get( *this, tmpFieldNum_ );
    }
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
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
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  template<class T>
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::operator=(const Daixt::Expr<T>&E)
  {
    for( int i=0; i<npts_; ++i )
      fieldValues_[i] = Evaluate(E);
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator+(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    check_tmp_validity();
    *tmp_=*this;
    *tmp_+=s;
    return *tmp_;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator-(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    check_tmp_validity();
    *tmp_=*this;
    *tmp_-=s;
    return *tmp_;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator*(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    check_tmp_validity();
    *tmp_=*this;
    *tmp_*=s;
    return *tmp_;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator/(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    check_tmp_validity();
    *tmp_=*this;
    *tmp_/=s;
    return *tmp_;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] = s.fieldValues_[i];
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator+=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] += s.fieldValues_[i];
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator-=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] -= s.fieldValues_[i];
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator *=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] *= s.fieldValues_[i];
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator /=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] /= s.fieldValues_[i];
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator =(const double a){
    for( int i=0; i<npts_; ++i ) fieldValues_[i] = a;
    return *this;
  } 
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator +=(const double a)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] += a;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator -=(const double a)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] -= a;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator *=(const double a)
  {
    for( int i=0; i<npts_; ++i ) fieldValues_[i] *= a;
    return *this;
  }
  //------------------------------------------------------------------
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

    typename SpatialField::iterator ifld = this->begin() + ixf;
    std::vector<double>::const_iterator irfld = r.begin() + ixr;

    const int yskip = ngxm+ngxp;
    const int zskip = nxf * (ngym+ngyp);

    for( int k=0; k<nzr; ++k ){
      for( int j=0; j<nyr; ++j ){
	for( int i=0; i<nxr; ++i ){
 	  *ifld++ = *irfld++;
	}
	ifld += yskip;
      }
      ifld += zskip;
    }
    return *this;
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator +=(const RHS& rhs)
  {
    // get the dimensions of the field
    const int nxf = extent_[0]
      + GhostTraits::template get<XDIR,SideMinus>()
      + GhostTraits::template get<XDIR,SidePlus>();

    const int nyf= ( extent_[1] > 1 ) ?
      extent_[1]
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

    static const int ngxm = GhostTraits::template get<XDIR,SideMinus>();
    static const int ngxp = GhostTraits::template get<XDIR,SidePlus>();
    static const int ngym = GhostTraits::template get<YDIR,SideMinus>();
    static const int ngyp = GhostTraits::template get<YDIR,SidePlus>();

    const int yskip = ngxm+ngxp;
    const int zskip = nxf * ( ngym+ngyp );

    int ixf = GhostTraits::template get<XDIR,SideMinus>();
    if( extent_[1] > 1 )   ixf += GhostTraits::template get<YDIR,SideMinus>();
    if( extent_[2] > 1 )   ixf += ( nxf * GhostTraits::template get<YDIR,SideMinus>() )
			     * ( nyf * GhostTraits::template get<ZDIR,SideMinus>() );

    typename SpatialField::iterator ifld = this->begin() + ixf;
    typename RHS::const_iterator irhs = rhs.begin();

    for( int k=0; k<nzr; ++k ){
      for( int j=0; j<nyr; ++j ){
	for( int i=0; i<nxr; ++i ){
	  *ifld++  +=  *irhs++;
	}
	ifld += yskip;
      }
      ifld += zskip;
    }
    return *this;
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField< VecOps, FieldLocation, GhostTraits >&
  SpatialField< VecOps, FieldLocation, GhostTraits >::
  operator -=(const RHS& rhs)
  {
    // get the dimensions of the field
    const int nxf = extent_[0]
      + GhostTraits::template get<XDIR,SideMinus>()
      + GhostTraits::template get<XDIR,SidePlus>();

    const int nyf= ( extent_[1] > 1 ) ?
      extent_[1]
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

    static const int ngxm = GhostTraits::template get<XDIR,SideMinus>();
    static const int ngxp = GhostTraits::template get<XDIR,SidePlus>();
    static const int ngym = GhostTraits::template get<YDIR,SideMinus>();
    static const int ngyp = GhostTraits::template get<YDIR,SidePlus>();

    static const int yskip = ngxm+ngxp;
    static const int zskip = nxf * ( ngym+ngyp );

    int ixf = GhostTraits::template get<XDIR,SideMinus>();
    if( extent_[1] > 1 )   ixf += GhostTraits::template get<YDIR,SideMinus>();
    if( extent_[2] > 1 )   ixf += ( nxf * GhostTraits::template get<YDIR,SideMinus>() )
			     * ( nyf * GhostTraits::template get<ZDIR,SideMinus>() );

    typename SpatialField::iterator ifld = this->begin() + ixf;
    typename RHS::const_iterator irhs = rhs.begin();

    for( int k=0; k<nzr; ++k ){
      for( int j=0; j<nyr; ++j ){
	for( int i=0; i<nxr; ++i ){
	  *ifld++  -=  *irhs++;
	}
	ifld += yskip;
      }
      ifld += zskip;
    }
    return *this;
  }
  //--------------------------------------------------------------------
  template< typename VecOps, typename FieldLocation, typename GhostTraits >
  void
  SpatialField<VecOps,FieldLocation,GhostTraits>::Print(std::ostream& s) const
  {
    vec_.Print(s);
  }
  //--------------------------------------------------------------------
  template< typename VecOps, typename FieldLocation, typename GhostTraits >
  bool
  SpatialField<VecOps,FieldLocation,GhostTraits>::operator==(const SpatialField& f)
  {
    for( int i=0; i<npts_; ++i ){
      if( f[i] != fieldValues_[i] ) return false; 
    }
    return true;
  }
  //--------------------------------------------------------------------
  template< typename VecOps, typename FieldLocation, typename GhostTraits >
  bool
  SpatialField<VecOps,FieldLocation,GhostTraits>::operator!=(const SpatialField& f)
  {
    return !(*this==f);
  }
  //--------------------------------------------------------------------


  //==================================================================


} // namespace SpatialOps

#endif
