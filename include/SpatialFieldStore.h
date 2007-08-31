#ifndef UT_SpatialFieldStore_h
#define UT_SpatialFieldStore_h

#include <queue>
#include <map>


namespace SpatialOps{


  // forward declaration
  template<typename T>  class SpatialFieldStore;

  /**
   *  @class  SpatFldPtr
   *  @author James C. Sutherland
   *  @date   June, 2007
   *
   *  @brief Wrapper for pointers to SpatialField objects.  Provides
   *  reference counting.
   *
   *  The SpatFldPtr class provides reference counted pointer
   *  management for SpatialField objects.  Furthermore, it can be
   *  used in conjunction with the SpatialFieldStore class to obtain
   *  temporary (scratch) SpatialField objects.  It also supports
   *  binary operations among SpatFldPtr objects - an extension to the
   *  SpatialField class, which doesn't support binary operators since
   *  the issue of creating intermediate/temporary objects is
   *  dangerous.
   *
   *  You should NOT dereference SpatFldPtr objects and store those
   *  references.  That is VERY dangerous and can lead to memory
   *  corruption.
   *
   *  Note that the underlying SpatialField object will not be owned
   *  by this class.  The real functionality provided here is to allow
   *  an interface to the SpatialFieldStore so that when a SpatFldPtr
   *  is created from there, it may be reference-counted so that when
   *  the last one is destroyed it returns control of the memory back
   *  to the SpatialFieldStore.  You can create SpatFldPtr onjects
   *  from an existing SpatialField, but the SpatFldPtr will not
   *  assume ownership of the memory.
   */
  template<typename FieldT>
  class SpatFldPtr
  {
  public:

    /**
     *  @brief Construct a SpatFldPtr.
     *
     *  @param field The field to wrap.  This constructor should be
     *  used if you want to wrap an existing SpatialField for use as a
     *  SpatFldPtr.
     */
    SpatFldPtr( FieldT& field );

    /**
     *  @brief Constructor for use from the SpatialFieldStore class only.
     *
     *  @param field The field to wrap.
     *
     *  @param builtFromStore if true, then SpatFldPtr will return the
     *  memory it owns to the SpatialFieldStore class once the last
     *  reference is destroyed.  If false then this will simply alias
     *  a FieldT object and provide binary operations.
     */
    SpatFldPtr( FieldT& field, const bool builtFromStore );

    ~SpatFldPtr();

    /** @brief Copy constructor */
    SpatFldPtr( const SpatFldPtr<FieldT>& p );

    /** @brief Skeletal constructor */
    SpatFldPtr();

    SpatFldPtr( FieldT* const p );

    /** @brief Assignment operator */
    SpatFldPtr& operator=( const SpatFldPtr& p );

    /** @brief Assignment operator */
    SpatFldPtr& operator=( FieldT& f );


    inline       FieldT& operator*()      {return *f_;}
    inline const FieldT& operator*() const{return *f_;}

    inline       FieldT* operator->()      {return f_;}
    inline const FieldT* operator->() const{return f_;}

    inline bool isnull() const{ return f_ == NULL; }

    /**
     *  @name binary Operators
     *
     *  These operators only result in new memory allocation when
     *  required, otherwise, a temporary is used from the
     *  SpatialFieldStore.  The resulting SpatFldPtr object should NOT
     *  be dereferenced and stored as a reference to an underlying
     *  SpatialField.  This will cause severe memory corruption.
     *
     *  These operator simply call through to the ones defined on the
     *  underlying SpatialField object.
     */
    //@{
    inline SpatFldPtr operator+(const SpatFldPtr& p) const{return (*f_ + *p);}  ///< Add two fields to produce a third: A=B+C
    inline SpatFldPtr operator-(const SpatFldPtr& p) const{return (*f_ - *p);}  ///< Subtract two fields to produce a third: A=B-C
    inline SpatFldPtr operator*(const SpatFldPtr& p) const{return (*f_ * *p);}  ///< Multiply two fields to produce a third: A=B*C
    inline SpatFldPtr operator/(const SpatFldPtr& p) const{return (*f_ / *p);}  ///< Divide two fields to produce a third: A=B/C
    //@}


    /**
     *  @name unary operators
     *  these simply call through to the corresponding SpatialField unary operators.
     */
    //@{
    inline SpatFldPtr& operator+=(const SpatFldPtr& p){*f_ += *p; return *this;}  ///< Add a SpatFldPtr to this.
    inline SpatFldPtr& operator-=(const SpatFldPtr& p){*f_ -= *p; return *this;}  ///< Subtract a SpatFldPtr from this.
    inline SpatFldPtr& operator*=(const SpatFldPtr& p){*f_ *= *p; return *this;}  ///< Multiply this by a SpatFldPtr
    inline SpatFldPtr& operator/=(const SpatFldPtr& p){*f_ /= *p; return *this;}  ///< Divide this by a SpatFldPtr

    inline SpatFldPtr& operator+=(const FieldT& p){*f_ += p; return *this;}  ///< Add a FieldT to this.
    inline SpatFldPtr& operator-=(const FieldT& p){*f_ -= p; return *this;}  ///< Subtract a FieldT from this.
    inline SpatFldPtr& operator*=(const FieldT& p){*f_ *= p; return *this;}  ///< Multiply this by a FieldT
    inline SpatFldPtr& operator/=(const FieldT& p){*f_ /= p; return *this;}  ///< Divide this by a FieldT

    inline SpatFldPtr& operator =(const double x){*f_  = x; return *this;}  ///< Assign this field to a constant
    inline SpatFldPtr& operator+=(const double x){*f_ += x; return *this;}  ///< Add a constant to this field
    inline SpatFldPtr& operator-=(const double x){*f_ -= x; return *this;}  ///< Subtract a constant from this field
    inline SpatFldPtr& operator*=(const double x){*f_ *= x; return *this;}  ///< Multiply this field by a constant
    inline SpatFldPtr& operator/=(const double x){*f_ /= x; return *this;}  ///< Divide this field by a constant
    //@}

  private:
    SpatialFieldStore<FieldT>& store_;
    FieldT* f_;
    int* count_;
    bool builtFromStore_;
  };



  /**
   *  @class  SpatialFieldStore
   *  @author James C. Sutherland
   *  @date   May, 2007
   *
   *  @brief Provides a common interface to obtain temporary (work) fields.
   *
   *  The SpatialFieldStore class provides a mechanism to generate
   *  temporary SpatialField objects. This prevents multiple
   *  allocation/deallocation of such objects that would be required
   *  otherwise.  It is implemented as a singleton, and provides a
   *  method:
   *
   *  \code
   *    SpatFldPtr<FieldT> SpatialFieldStore<FieldT>::get( const FieldT& f )
   *  \endcode
   *
   *  to return a field with the asme dimensions as the provided
   *  template field.  Note that the field will not necessarily have
   *  the same values as the provided field.  The supplied field is
   *  simply used to provide information needed to construct clones.
   *
   *  Note that the returned type, <code>SpatFldPtr<FieldT></code>,
   *  should not be dereferenced and saved as a SpatialField.  Doing
   *  so can cause serious memory corruption.
   *
   *
   *  @par Thread-Parallelism Issues:
   *
   *  NOTE: this could get us into big trouble if we have threads
   *  running concurrently, since we would not be able to guarantee
   *  that two threads didn't use the same memory.
   *
   *  \todo Implement a thread-safe version of this concept.
   */
  template< typename FieldT >
  class SpatialFieldStore
  {
    friend class SpatFldPtr<FieldT>;

  public:
    static SpatialFieldStore& self();


    /**
     *  @brief Obtain a temporary field.
     *
     *  @param f  A field to model this one after.
     *
     *  Note that you should not dereference the SpatFldPtr object to
     *  store a SpatialField reference.  Doing so can cause memory
     *  corruption.
     */
    inline SpatFldPtr<FieldT> get( const FieldT& f );


    inline SpatFldPtr<FieldT> get( const int ntot,
				   const std::vector<int>& entriesPerComponent,
				   const std::set<int>& ghostSet );

  private:

    /**
     *  @brief Restores a field to the store for future use.
     *
     *  Note that this method is private to ensure it is only called
     *  by SpatFldPtr objects.  Calling it anywhere else can result in
     *  memory corruption.
     */
    void restore_field( FieldT& f );

    SpatialFieldStore(){};
    ~SpatialFieldStore();


    typedef std::queue<FieldT*> FieldQueue;
    typedef std::map<int,FieldQueue> FQMap;

    FQMap fqmap_;
  };


  //==================================================================






  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






  //=================================================================


  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::SpatFldPtr( FieldT& f )
    : store_( SpatialFieldStore<FieldT>::self() )
  {
    f_ = &f;
    count_ = new int;
    *count_ = 1;
    builtFromStore_ = false;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::SpatFldPtr( FieldT& f, const bool builtFromStore )
    : store_( SpatialFieldStore<FieldT>::self() )
  {
    f_ = &f;
    count_ = new int;
    *count_ = 1;
    builtFromStore_ = builtFromStore;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::SpatFldPtr( const SpatFldPtr<FieldT>& p )
    : store_( SpatialFieldStore<FieldT>::self() )
  {
    f_ = p.f_;
    count_ = p.count_;
    ++(*count_);
    builtFromStore_ = p.builtFromStore_;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::SpatFldPtr()
    : store_( SpatialFieldStore<FieldT>::self() )
  {
    f_     = NULL;
    count_ = NULL;
    builtFromStore_ = false;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::SpatFldPtr( FieldT* const p )
    : store_( SpatialFieldStore<FieldT>::self() )
  {
    f_ = p;
    *count_ = 1;
    builtFromStore_ = false;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>&
  SpatFldPtr<FieldT>::operator=( const SpatFldPtr& p )
  {
    // was this an active SpatFldPtr?
    if( count_ != NULL ){
      // this one is dying so decrement the count.
      --(*count_);
      // kill the old one if needed
      if( *count_ == 0 ){
	if( builtFromStore_ ) store_.restore_field( *f_ );
	delete count_;
      }
    }
    // reassign
    f_ = p.f_;
    count_ = p.count_;
    builtFromStore_ = p.builtFromStore_;
    // increment copy count
    ++(*count_);

    return *this;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>&
  SpatFldPtr<FieldT>::operator=( FieldT& f )
  {
    // was this an active SpatFldPtr?
    if( count_ != NULL ){
      // this one is dying so decrement the count.
      --(*count_);
      // kill the old one if needed
      if( *count_ == 0 ){
	if( builtFromStore_ ) store_.restore_field( *f_ );
	delete count_;
      }
    }
    // reassign
    f_ = &f;
    count_ = new int;
    *count_ = 1;
    builtFromStore_ = false;

    return *this;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::~SpatFldPtr()
  {
    if( count_ != NULL ){
      --(*count_);
      if( *count_ == 0 ){
	if( builtFromStore_ ) store_.restore_field( *f_ );
	delete count_;
      }
    }
  }
  //------------------------------------------------------------------
  

  //==================================================================


  //------------------------------------------------------------------
  template<typename FieldT>
  SpatialFieldStore<FieldT>::~SpatialFieldStore()
  {
    for( typename FQMap::iterator ii=fqmap_.begin(); ii!=fqmap_.end(); ++ii ){
      FieldQueue& q = ii->second;
      while( !q.empty() ){
	FieldT* field = q.front();
	delete field;
	q.pop();
      }
    }
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatialFieldStore<FieldT>&
  SpatialFieldStore<FieldT>::self()
  {
    static SpatialFieldStore<FieldT> s;
    return s;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>
  SpatialFieldStore<FieldT>::get( const FieldT& f )
  {
    // find the proper map
    FieldQueue& q = fqmap_[ f.get_ntotal() ];

    if( q.empty() ){
      FieldT* fnew = new FieldT( f );
      q.push( fnew );
    }

    FieldT* fnew = q.front();
    q.pop();

    return SpatFldPtr<FieldT>(*fnew,true);
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>
  SpatialFieldStore<FieldT>::get( const int ntot,
				  const std::vector<int>& entriesPerComponent,
				  const std::set<int>& ghostSet )
  {
    // find the proper map
    FieldQueue& q = fqmap_[ ntot ];

    if( q.empty() ){
      FieldT* fnew = new FieldT( ntot, entriesPerComponent, ghostSet, NULL );
      q.push( fnew );
    }

    FieldT* fnew = q.front();
    q.pop();

    return SpatFldPtr<FieldT>(*fnew,true);
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  void
  SpatialFieldStore<FieldT>::restore_field( FieldT& field )
  {
    FieldQueue& q = fqmap_[ field.get_ntotal() ];
    q.push( &field );
  }
  //------------------------------------------------------------------


  //====================================================================


} // namespace SpatialOps

#endif
