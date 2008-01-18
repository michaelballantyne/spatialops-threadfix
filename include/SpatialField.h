#ifndef UT_SpatialField_h
#define UT_SpatialField_h

#include <vector>
#include <set>
#include <iostream>

// allows compile-time expansion of complex expressions involving SpatialField objects.
#include <daixtrose/Daixt.h>

#include <SpatialFieldStore.h>
#include <SFIterator.h>

// #include <LinearSystem.h>


class RHS;

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
    typedef typename VecOps::VecType VecType;

    typedef SpatialField<VecOps,FieldLocation,GhostTraits> MyType;

    friend class SpatialFieldStore<MyType>;


  public:

    typedef GhostTraits Ghost;
    typedef FieldLocation Location;


    /**
     *  Construct a SpatialField.
     *
     *  @param npts The number of points (including ghost cells) for
     *  this field.
     *
     *  @param ghostSet A std::set<int> containing the indices
     *  representing the ghost values for this field.
     *
     *  @param fieldValues  Pointer to the field values.  Behavior is
     *  dictated by the choice of StorageMode.
     *
     *  @param mode  Storage options.  If InternalStorage then the
     *  fieldValues will be copied into an internal buffer.  If
     *  ExternalStorage then the fieldValues will be stored
     *  externally.  Efficiency suggests that ExternalStorage is best,
     *  since it will avoid excessive copies.  Safety suggests that
     *  InternalStorage is best, since it protects against memory
     *  corruption and inadvertant deletion of the field's underlying
     *  memory.
     */
    SpatialField( const int npts,
		  const std::set<int>& ghostSet,
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
     *  @name Binary Operators for SpatialField objects.
     *
     *  These all return a SpatFldPtr object.  Since binary operations
     *  require a new field to hold the result, a SpatFldPtr is used.
     *  This only results in new memory allocation when required,
     *  otherwise, a temporary is used from the SpatialFieldStore.
     *  The resulting SpatFldPtr object should NOT be dereferenced and
     *  stored as a reference to an underlying SpatialField.  This
     *  will cause severe memory corruption.
     */
    //@{
    inline SpatFldPtr<SpatialField> operator+(const SpatialField&) const;  ///< Add two fields to produce a third: A=B+C
    inline SpatFldPtr<SpatialField> operator-(const SpatialField&) const;  ///< Subtract two fields to produce a third: A=B-C
    inline SpatFldPtr<SpatialField> operator*(const SpatialField&) const;  ///< Multiply two fields to produce a third: A=B*C
    inline SpatFldPtr<SpatialField> operator/(const SpatialField&) const;  ///< Divide two fields to produce a third: A=B/C
    //@}


    /**
     *  @name Unary Operators for SpatialField objects
     *
     *  Perform basic operations on a SpatialField.
     */
    //@{

    template<typename FieldT> SpatialField& operator= (const FieldT&); ///< No default implementation...
    template<typename FieldT> SpatialField& operator*=(const FieldT&); ///< No default implementation...
    template<typename FieldT> SpatialField& operator/=(const FieldT&); ///< No default implementation...
    template<typename FieldT> SpatialField& operator+=(const FieldT&); ///< No default implementation...
    template<typename FieldT> SpatialField& operator-=(const FieldT&); ///< No default implementation...

    inline SpatialField& operator= (const SpatFldPtr<SpatialField>&);  ///< Assign a SpatialField to this one.

    inline SpatialField& operator= (const SpatialField&);  ///< Assign a SpatialField to this one.
    inline SpatialField& operator+=(const SpatialField&);  ///< Add a SpatialField to this.
    inline SpatialField& operator-=(const SpatialField&);  ///< Subtract a SpatialField from this.
    inline SpatialField& operator*=(const SpatialField&);  ///< Multiply this by a SpatialField
    inline SpatialField& operator/=(const SpatialField&);  ///< Divide this by a SpatialField

    inline SpatialField& operator =(const double);  ///< Assign this field to a constant
    inline SpatialField& operator+=(const double);  ///< Add a constant to this field
    inline SpatialField& operator-=(const double);  ///< Subtract a constant from this field
    inline SpatialField& operator*=(const double);  ///< Multiply this field by a constant
    inline SpatialField& operator/=(const double);  ///< Divide this field by a constant

    inline SpatialField& operator= (const RHS&);  ///< Assign a RHS to this field (doesn't affect ghosts)
    inline SpatialField& operator+=(const RHS&);  ///< Add a RHS to this field (doesn't affect ghosts)
    inline SpatialField& operator-=(const RHS&);  ///< Subtract a RHS from this field (doesn't affect ghosts)

    inline bool operator==(const SpatialField&);  ///< Is this field equal to the supplied one?
    inline bool operator!=(const SpatialField&);  ///< Is this field not equal to the supplied one?

    //@}


    /**
     * @brief This provides support for Daixtrose - an expression template
     * engine to allow compund expressions involving SpatialField
     * objects to be unrolled by the compiler.
     */
    template<class T>
    inline SpatialField& operator=(const Daixt::Expr<T>&E);

    //@}


    /**
     * @name
     * Obtain the underlying VecType object that corresponds to
     * the LinAlg strategy.
     */
    //@{
    inline       VecType & get_linalg_vec()      { return vec_; }
    inline const VecType & get_linalg_vec() const{ return vec_; }
    //@}

    /**
     * @brief Get the total number of points (including ghost layers)
     * in this SpatialField.
     */
    inline int get_ntotal() const{ return npts_; }

    inline int get_ninterior() const{return npts_-ghostSet_.size(); }
    /**
     *  @name
     *  Obtain a reference to the field using the [] operator.
     *  This should not generally be used, as it is not tuned for
     *  performance.
     */
    //@{
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



    typedef  InteriorIterator<      double*>       interior_iterator;
    typedef  InteriorIterator<const double*> const_interior_iterator;

    inline       interior_iterator interior_begin();
    inline const_interior_iterator interior_begin() const;

    inline       interior_iterator interior_end();
    inline const_interior_iterator interior_end() const;



    typedef GhostIterator<      double*>       ghost_iterator;
    typedef GhostIterator<const double*> const_ghost_iterator;

    inline       ghost_iterator ghost_begin();
    inline const_ghost_iterator ghost_begin() const;

    inline       ghost_iterator ghost_end();
    inline const_ghost_iterator ghost_end() const;

    //@}


    int get_entries_per_comp( const int dir ) const{return nptsPerComp_[dir];}

    /** Dump information about the field to the given output stream. */
    virtual void Print( std::ostream& ) const;

  protected:
    

  private:

    VecOps linAlg_;
    const int npts_;
    std::vector<int> nptsPerComp_;
    const std::set<int> ghostSet_;
    const StorageMode storageMode_;
    double * const fieldValues_;
    VecType & vec_;

    SpatialFieldStore<SpatialField>& sfStore_;

    SpatialField( const SpatialField& );
    SpatialField();
    
  };


  //==================================================================





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
    std::cout << OP::Symbol() << std::endl;
    Evaluate( arg.rhs() );
  }
  //------------------------------------------------------------------
  template<typename ARG, typename OP>
  void Evaluate( const Daixt::UnOp<ARG,OP>& arg )
  {
    std::cout << OP::Symbol() << std::endl;
    Evaluate(arg.arg());
  }
  //------------------------------------------------------------------
  template<typename T>
  void Evaluate( const T& t )
  {
    t.put(std::cout);
    std::cout << std::endl;
  }
  //------------------------------------------------------------------


  //==================================================================


  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::interior_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  interior_begin()
  {
    iterator iter = begin();
    int ix=0;
    std::set<int>::const_iterator ig = ghostSet_.begin();
    const std::set<int>::const_iterator ige = ghostSet_.end();
    while( *ig == ix  &&  ig!=ige ){ ++ig; ++ix; ++iter; }
    return interior_iterator(iter,ix,ig,ige);
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::const_interior_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  interior_begin() const
  {
    const_iterator iter = begin();
    int ix=0;
    std::set<int>::const_iterator ig = ghostSet_.begin();
    const std::set<int>::const_iterator ige = ghostSet_.end();
    while( *ig == ix && ig!=ige ){ ++ig; ++ix; ++iter; }
    return const_interior_iterator(iter,ix,ig,ige);
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::interior_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  interior_end()
  {
    iterator iter = end();
    std::set<int>::const_iterator ig = ghostSet_.end();
    const std::set<int>::const_iterator ige = ghostSet_.begin();
    int ix=npts_;
    if( ig != ige ){
      --ig;
      while( *ig == ix && ig!=ige ){ --ig; --ix; --iter; }
    }
    return interior_iterator(iter,ix,ig,ige);
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::const_interior_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  interior_end() const
  {
    const_iterator iter = end();
    std::set<int>::const_iterator ig = ghostSet_.end();
    const std::set<int>::const_iterator ige = ghostSet_.begin();
    int ix=npts_;
    if( ig != ige ){
      --ig;
      while( *ig == ix && ig!=ige ){ --ig; --ix; --iter; }
    }
    return const_interior_iterator(iter,ix,ig,ige);
  }

  //==================================================================


  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::ghost_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  ghost_begin()
  {
    // jcs this may be wrong
    const int ig=*(ghostSet_.begin());
    return ghost_iterator(fieldValues_+ig,ig,ghostSet_);
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::const_ghost_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  ghost_begin() const
  {
    // jcs this may be wrong
    const int ig=*(ghostSet_.begin());
    return const_ghost_iterator(fieldValues_+ig,ig,ghostSet_);
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::ghost_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  ghost_end()
  {
    // jcs this may be wrong
    const int ig=*(ghostSet_.end()-1);
    return ghost_iterator(fieldValues_+ig,ig,ghostSet_);
  }
  //--------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  typename SpatialField<VecOps,FieldLocation,GhostTraits>::const_ghost_iterator
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  ghost_end() const
  {
    // jcs this may be wrong
    const int ig=*(ghostSet_.end()-1);
    return const_ghost_iterator(fieldValues_+ig,ig,ghostSet_);
  }
  //--------------------------------------------------------------------


  //==================================================================


  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  SpatialField( const int npts,
		const std::set<int>& ghostSet,
		double * const fieldValues,
		const StorageMode mode )
    : npts_( npts ),
      nptsPerComp_(3,0),
      ghostSet_( ghostSet ),
      storageMode_( mode ),
      
      fieldValues_( (storageMode_==ExternalStorage)
		    ? fieldValues
		    : new double[npts_] ),

      vec_( linAlg_.setup_vector( npts_, fieldValues_ ) ),

      sfStore_( SpatialFieldStore<SpatialField>::self() )
  {
    nptsPerComp_[0] = npts;

    if( mode==InternalStorage )  reset_values( npts_, fieldValues );
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  ~SpatialField()
  {
    if( storageMode_ == InternalStorage )  delete [] fieldValues_;
    linAlg_.destroy_vector();
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  SpatialField( const SpatialField<VecOps,FieldLocation,GhostTraits>& f )
    : npts_       ( f.npts_ ),
      nptsPerComp_( f.nptsPerComp_ ),
      ghostSet_   ( f.ghostSet_ ),
      storageMode_( InternalStorage ),
      fieldValues_( new double[npts_] ),
      vec_        ( linAlg_.setup_vector(npts_,fieldValues_) ),
      sfStore_( SpatialFieldStore<SpatialField>::self() )
  {
    reset_values( npts_, f.fieldValues_ );
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  void
  SpatialField<VecOps,FieldLocation,GhostTraits>::
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
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::operator=(const Daixt::Expr<T>&E)
  {
    for( int i=0; i<npts_; ++i )
      fieldValues_[i] = Evaluate(E);
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator= (const SpatFldPtr<SpatialField>& s)
  {
    return *this = *s;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator=(const SpatialField& s)
  {
    assert( npts_ == s.npts_ );
    typename SpatialField::iterator ifld = this->begin();
    const typename SpatialField::const_iterator iend = this->end();
    typename SpatialField::const_iterator isrc = s.begin();
    for( ; ifld!=iend; ++ifld, ++isrc ) *ifld = *isrc;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator+=(const SpatialField& s)
  {
    typename SpatialField::iterator ifld = this->begin();
    const typename SpatialField::const_iterator iend = this->end();
    typename SpatialField::const_iterator isrc = s.begin();
    for( ; ifld!=iend; ++ifld, ++isrc ) *ifld += *isrc;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>& 
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator-=(const SpatialField& s)
  {
    typename SpatialField::iterator ifld = this->begin();
    const typename SpatialField::const_iterator iend = this->end();
    typename SpatialField::const_iterator isrc = s.begin();
    for( ; ifld!=iend; ++ifld, ++isrc ) *ifld -= *isrc;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>& 
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator*=(const SpatialField& s)
  {
    typename SpatialField::iterator ifld = this->begin();
    const typename SpatialField::const_iterator iend = this->end();
    typename SpatialField::const_iterator isrc = s.begin();
    for( ; ifld!=iend; ++ifld, ++isrc ) *ifld *= *isrc;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>& 
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator/=(const SpatialField& s)
  {
    typename SpatialField::iterator ifld = this->begin();
    const typename SpatialField::const_iterator iend = this->end();
    typename SpatialField::const_iterator isrc = s.begin();
    for( ; ifld!=iend; ++ifld, ++isrc ) *ifld /= *isrc;
    return *this;
  }

  //------------------------------------------------------------------

  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator=(const double a){
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ = a;
    return *this;
  } 
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator+=(const double a)
  {
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ += a;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator-=(const double a)
  {
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ -= a;
    return *this;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatialField<VecOps,FieldLocation,GhostTraits>&
  SpatialField<VecOps,FieldLocation,GhostTraits>::
  operator*=(const double a)
  {
    double* f = fieldValues_;
    for( int i=0; i<npts_; ++i ) *f++ *= a;
    return *this;
  }
  //------------------------------------------------------------------
  template< typename VecOps, typename FieldLocation, typename GhostTraits >
  void
  SpatialField<VecOps,FieldLocation,GhostTraits>::Print(std::ostream& s) const
  {
    linAlg_.print_vec(s);
  }
  //------------------------------------------------------------------
  template< typename VecOps, typename FieldLocation, typename GhostTraits >
  bool
  SpatialField<VecOps,FieldLocation,GhostTraits>::operator==(const SpatialField& f)
  {
    for( int i=0; i<npts_; ++i ){
      if( f[i] != fieldValues_[i] ) return false; 
    }
    return true;
  }
  //------------------------------------------------------------------
  template< typename VecOps, typename FieldLocation, typename GhostTraits >
  bool
  SpatialField<VecOps,FieldLocation,GhostTraits>::operator!=(const SpatialField& f)
  {
    return !(*this==f);
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatFldPtr< SpatialField< VecOps, FieldLocation,GhostTraits > >
  SpatialField< VecOps, FieldLocation,GhostTraits >::  
  operator+( const SpatialField& f ) const
  {
    SpatFldPtr<SpatialField> result( sfStore_.get(*this) );
    *result  = *this;
    *result += f;
    return result;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatFldPtr< SpatialField< VecOps, FieldLocation,GhostTraits > >
  SpatialField< VecOps, FieldLocation,GhostTraits >::  
  operator-( const SpatialField& f ) const
  {
    SpatFldPtr<SpatialField> result( sfStore_.get(*this) );
    *result  = *this;
    *result -= f;
    return result;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatFldPtr< SpatialField< VecOps, FieldLocation,GhostTraits > >
  SpatialField< VecOps, FieldLocation,GhostTraits >::  
  operator*( const SpatialField& f ) const
  {
    SpatFldPtr<SpatialField> result( sfStore_.get(*this) );
    *result  = *this;
    *result *= f;
    return result;
  }
  //------------------------------------------------------------------
  template< class VecOps, typename FieldLocation, typename GhostTraits >
  SpatFldPtr< SpatialField< VecOps, FieldLocation,GhostTraits > >
  SpatialField< VecOps, FieldLocation,GhostTraits >::  
  operator/( const SpatialField& f ) const
  {
    SpatFldPtr<SpatialField> result( sfStore_.get(*this) );
    *result  = *this;
    *result /= f;
    return result;
  }
  
  //------------------------------------------------------------------

  //==================================================================


} // namespace SpatialOps

#endif
