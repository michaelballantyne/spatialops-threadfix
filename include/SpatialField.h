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


  /**
   *  @enum StorageMode
   *  @brief Enumerates options for storage of a SpatialField.
   */
  enum StorageMode
    {
      InternalStorage, ///< Selector for storing fields internal to the SpatialField object.
      ExternalStorage  ///< Selector for storing fields external to the SpatialField object.
    };


  //==================================================================


  /**
   *  @class SpatialField
   *  @author James C. Sutherland
   *  @date   December, 2006
   *
   *  @brief Class to represent spatial fields defined on a logically
   *  rectangular domain.
   *
   *  @par Template Parameters
   *
   *   <ul>
   *
   *   <li> \b VecOps a policy dictating the treatment of this
   *   SpatialField by the linear algebra package. Defines a \b
   *   VecType type that is the type of object that the underlying
   *   vector is stored in.
   *
   *   <li> \b FieldLocation a trait specifying the field location
   *   type (e.g. node, cell, face, etc.)
   *
   *   <li> \b GhostTraits Defines information about ghosting.  This
   *   should provide a templated method,
   *
   *     \code
   *       template<typename Dir, typename SideType> static int get();
   *     \endcode
   *
   *   which returns the number of ghost cells in a given direction
   *   and given side of the patch.  This must be a static method, and
   *   may be specialized to deal with different ghosting on different
   *   faces of a patch.
   *
   *   </ul>
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
     *  dictated by the choice of StorageMode.
     *
     *  @param mode : Storage options.  If InternalStorage then the
     *  fieldValues will be copied into an internal buffer.  If
     *  ExternalStorage then the fieldValues will be stored
     *  externally.  Efficiency suggests that ExternalStorage is best,
     *  since it will avoid excessive copies.  Safety suggests that
     *  InternalStorage is best, since it protects against memory
     *  corruption and inadvertant deletion of the field's underlying
     *  memory.
     */
    SpatialField( const std::vector<int> & fieldDims,
		  double * const fieldValues,
		  const StorageMode mode = InternalStorage );


    virtual ~SpatialField();


    /**
     *  @brief Overwrite the values in the SpatialField with the ones supplied.
     *  @param npts : number of points (including ghost cells)
     *  @param values : array of values to overwrite with.
     */
    inline void reset_values( const int npts,
			      const double* const values );


    /**
     *  @name Operators for SpatialField objects
     */
    //@{

    inline SpatialField& operator =(const SpatialField&);  ///< Assign a SpatialField to this one.
    inline SpatialField& operator+=(const SpatialField&);  ///< Add a SpatialField to this.
    inline SpatialField& operator-=(const SpatialField&);  ///< Subtract a SpatialField from this.
    inline SpatialField& operator*=(const SpatialField&);  ///< Multiply this by a SpatialField
    inline SpatialField& operator/=(const SpatialField&);  ///< Divide this by a SpatialField

    inline SpatialField& operator =(const double);  ///< Assign this field to a constant
    inline SpatialField& operator+=(const double);  ///< Add a constant to this field
    inline SpatialField& operator-=(const double);  ///< Subtract a constant from this field
    inline SpatialField& operator*=(const double);  ///< Multiply this field by a constant
    inline SpatialField& operator/=(const double);  ///< Divide this field by a constant

    inline SpatialField& operator =(const RHS&);  ///< Assign a RHS to this field (doesn't affect ghosts)
    inline SpatialField& operator+=(const RHS&);  ///< Add a RHS to this field (doesn't affect ghosts)
    inline SpatialField& operator-=(const RHS&);  ///< Subtract a RHS from this field (doesn't affect ghosts)

    inline bool operator==(const SpatialField&);  ///< Is this field equal to the supplied one?
    inline bool operator!=(const SpatialField&);  ///< Is this field not equal to the supplied one?

    //@}


    /**
     *  @name Binary Operators
     *
     *  note that these return references rather than copies because
     *  we have a temporary working vector that we use internally.
     *  However, this also means that it is not safe to get a
     *  reference to these variables, since they could change very
     *  easily!  That means that you SHOULD NOT do something like
     *
     *  \code
     *     SpatialField a,b;
     *     SpatialField & c = a+b;
     *  \endcode
     *
     *  rather, you should do:
     *
     *  \code
     *     SpatialField a,b,c;
     *     c = a+b;
     *  \endcode
     *
     *  This results in a copy from a temporary to c, but is safe.
     *
     *  NOTE: this could get us into big trouble if we have threads
     *  running concurrently, since we would not be able to guarantee
     *  that two threads didn't use the same memory.  We could
     *  implement some sort of locking on the tmp fields, but this
     *  would require us to obtain a new field for each call that
     *  required it - a bit slower but probably worth the price.
     */
    //@{

    inline SpatialField& operator+(const SpatialField&) const;
    inline SpatialField& operator-(const SpatialField&) const;
    inline SpatialField& operator*(const SpatialField&) const;
    inline SpatialField& operator/(const SpatialField&) const;

    /**
     * this provides support for Daixtrose - an expression template
     * engine to allow compund expressions involving SpatialField
     * objects to be unrolled by the compiler.
     */
    template<class T>
    inline SpatialField& operator=(const Daixt::Expr<T>&E);

    //@}


    //@{
    /**
     * @brief Obtain the underlying VecType object that corresponds to
     * the LinAlg strategy.
     */
    inline       VecType & get_linalg_vec()      { return vec_; }
    inline const VecType & get_linalg_vec() const{ return vec_; }
    //@}

    /**
     * @brief Get the total number of points (including ghost layers)
     * in this SpatialField.
     */
    inline int get_ntotal() const{ return npts_; }

    /**
     * @brief Obtain the number of ghost cells in the given direction
     * and side of the patch.
     */
    template< typename Dir, typename SideType>
    static int nghost(){ return GhostTraits::template get<Dir,SideType>(); }

    /**
     *  @brief Obtain the domain extent.  This is a vector containing
     *  the number of points in each direction, excluding ghost cells.
     */
    inline const std::vector<int>& get_extent() const{return extent_;}


    //@{
    /**
     *  @bfief Obtain a reference to the field using the [] operator.
     *  This should not generally be used, as it is not tuned for
     *  performance.
     */
    inline       double& operator[](const int i)      { return fieldValues_[i]; }
    inline const double& operator[](const int i) const{ return fieldValues_[i]; }
    //@}


    /**
     *  @name Iterators for SpatialField objects
     */
    //@{

    typedef double*           iterator;
    typedef double const*     const_iterator;

    inline iterator       begin()      {return fieldValues_;}
    inline const_iterator begin() const{return fieldValues_;}

    inline iterator       end()      {return fieldValues_+npts_;}
    inline const_iterator end() const{return fieldValues_+npts_;}

    //@}


    /** Dump information about the field to the given output stream. */
    void Print( std::ostream& ) const;

  protected:
    
    inline bool consistency_check( const SpatialField& s ) const{ return ( npts_ == s.npts_ ); }

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
    inline void check_tmp_validity() const;



    VecOps linAlg_;
    const std::vector<int> extent_;
    const int npts_;
    const StorageMode storageMode_;
    double * const fieldValues_;
    VecType & vec_;

    /**
     * these fields facilitate more efficient operations.  +,-,*,/
     * operators require a temporary.  Here we store a temporary for
     * this purpose rather than building one each time we hit one of
     * these operators.  The SpatialFieldStore is used to assist in
     * holding and building temporaries to minimize the number of them
     * that are created.  Basically, any time any of the above
     * operators are called, a temporary must be generated.  Note that
     * the increment operators like += *= etc do not require
     * temporaries.
     */
    //@{
    mutable SpatialField* tmp_;             ///< a temporary field to use for efficient operations
    mutable size_t tmpFieldNum_; ///< the id for the temporary field...
    //@}

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
    // leave the tmp field NULL.  If we need it later, we will construct it then.
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
  SpatialField<VecOps,FieldLocation, GhostTraits>::check_tmp_validity() const
  {
    if( tmp_ == NULL ){
      tmp_ = &SpatialFieldStore<VecOps,FieldLocation,GhostTraits>::self().get( *this, tmpFieldNum_++ );
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
  operator+(const SpatialField<VecOps,FieldLocation, GhostTraits>& s) const
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
  operator-(const SpatialField<VecOps,FieldLocation, GhostTraits>& s) const
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
  operator*(const SpatialField<VecOps,FieldLocation, GhostTraits>& s) const
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
  operator/(const SpatialField<VecOps,FieldLocation, GhostTraits>& s) const
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
    double* f = fieldValues_;
    const double* sf = s.fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ = *sf++;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator+=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    double* f = fieldValues_;
    const double* sf = s.fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ += *sf++;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator-=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    double* f = fieldValues_;
    const double* sf = s.fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ -= *sf++;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator*=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    double* f = fieldValues_;
    const double* sf = s.fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ *= *sf++;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>& 
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator/=(const SpatialField<VecOps,FieldLocation, GhostTraits>& s)
  {
    double* f = fieldValues_;
    const double* sf = s.fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ /= *sf++;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator=(const double a){
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ = a;
    return *this;
  } 
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator+=(const double a)
  {
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ += a;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator-=(const double a)
  {
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ -= a;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation, GhostTraits>&
  SpatialField<VecOps,FieldLocation, GhostTraits>::
  operator*=(const double a)
  {
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ *= a;
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
  operator=(const RHS& rhs)
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

    static const int yskip = ngxm+ngxp;
    static const int zskip = nxf * (ngym+ngyp);

    int ixf = GhostTraits::template get<XDIR,SideMinus>();
    if( extent_[1] > 1 ){
      ixf = nxf + GhostTraits::template get<YDIR,SideMinus>();
    }
    if( extent_[2] > 1 ){       // must also have extent_[1] > 1
      ixf += nxf*nyf;
    }

    typename SpatialField::iterator ifld = this->begin() + ixf;
    std::vector<double>::const_iterator irfld = r.begin();

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
  operator+=(const RHS& rhs)
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
    if( extent_[1] > 1 ){
      ixf = nxf + GhostTraits::template get<YDIR,SideMinus>();
    }
    if( extent_[2] > 1 ){       // must also have extent_[1] > 1
      ixf += nxf*nyf;
    }

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
  operator-=(const RHS& rhs)
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
    if( extent_[1] > 1 ){
      ixf = nxf + GhostTraits::template get<YDIR,SideMinus>();
    }
    if( extent_[2] > 1 ){       // must also have extent_[1] > 1
      ixf += nxf*nyf;
    }

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
