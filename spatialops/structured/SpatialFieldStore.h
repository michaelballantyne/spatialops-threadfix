#ifndef UT_SpatialFieldStore_h
#define UT_SpatialFieldStore_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/SpatialField.h>

#include <queue>
#include <map>
#include <set>

#include <boost/type_traits.hpp>

#ifdef ENABLE_THREADS
# include <boost/thread/mutex.hpp>
#endif

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
     *  SpatFldPtr.  Ownership of this pointer is transfered.
     */
    SpatFldPtr( FieldT* const field );

    /**
     *  @brief Constructor for use from the SpatialFieldStore class only.
     *
     *  @param field The field to wrap.  Ownership is transfered.
     *
     *  @param builtFromStore if true, then SpatFldPtr will return the
     *  memory it owns to the SpatialFieldStore class once the last
     *  reference is destroyed.  If false then this will simply alias
     *  a FieldT object.
     */
    SpatFldPtr( FieldT* const field, const bool builtFromStore );

    ~SpatFldPtr();

    /** @brief Copy constructor */
    SpatFldPtr( const SpatFldPtr<FieldT>& p );

    /** @brief Skeletal constructor */
    SpatFldPtr();

    /** @brief Assignment operator */
    SpatFldPtr& operator=( const SpatFldPtr& p );

    /** @brief Assignment operator */
    SpatFldPtr& operator=( FieldT* const f );


    inline       FieldT& operator*()      {return *f_;}
    inline const FieldT& operator*() const{return *f_;}

    inline       FieldT* operator->()      {return f_;}
    inline const FieldT* operator->() const{return f_;}

    inline bool isnull() const{ return f_ == NULL; }


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

    const int getCount() { return *count_; }
    const bool getBFS() { return builtFromStore_; }

    void free();
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
   *    SpatFldPtr<FieldT> field = SpatialFieldStore<FieldT>::get( const FieldT& f )
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

    inline SpatFldPtr<FieldT> get( const structured::MemoryWindow& window );

  private:

    /**
     *  @brief Restores a field to the store for future use.
     *
     *  Note that this method is private to ensure it is only called
     *  by SpatFldPtr objects.  Calling it anywhere else can result in
     *  memory corruption.
     */
    inline void restore_field( FieldT& f );

#ifdef ENABLE_THREADS
    /**
     *  Used to lock threads to prevent simultaneous access.
     */
    inline boost::mutex& get_mutex(){ static boost::mutex m; return m; }
#endif

    SpatialFieldStore(){};
    ~SpatialFieldStore();

    template< typename FT, typename IsPODT > struct ValTypeSelector;
    template< typename FT > struct ValTypeSelector<FT,boost::true_type >{ typedef FT type; };
    template< typename FT > struct ValTypeSelector<FT,boost::false_type>{ typedef typename FT::AtomicT type; };

    typedef typename ValTypeSelector< FieldT, typename boost::is_pod<FieldT>::type >::type AtomicT;
    typedef std::queue<AtomicT*> FieldQueue;
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
  SpatFldPtr<FieldT>::SpatFldPtr( FieldT* const f )
    : store_( SpatialFieldStore<FieldT>::self() ),
      f_( f ),
      count_( new int ),
      builtFromStore_( false )
  {
    *count_ = 1;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::SpatFldPtr( FieldT* const f, const bool builtFromStore )
    : store_( SpatialFieldStore<FieldT>::self() ),
      f_( f ),
      count_( new int ),
      builtFromStore_( builtFromStore )
  {
    *count_ = 1;
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
			count_=NULL;
        delete f_;
			f_=NULL;
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
  SpatFldPtr<FieldT>::operator=( FieldT* const f )
  {
    // was this an active SpatFldPtr?
    if( count_ != NULL ){
      // this one is dying so decrement the count.
      --(*count_);
      // kill the old one if needed
      if( *count_ == 0 ){
        if( builtFromStore_ ) store_.restore_field( *f_ );
        delete count_;
			count_=NULL;
        delete f_;
			f_=NULL;
      }
    }
    // reassign
    f_ = f;
    count_ = new int;
    *count_ = 1;
    builtFromStore_ = false;

    return *this;
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  void
  SpatFldPtr<FieldT>::free()
  {
    // was this an active SpatFldPtr?
    if( count_ != NULL ){
      // this one is dying so decrement the count.
      --(*count_);
      if( *count_ == 0 ){
         // kill the old one if needed
        if( builtFromStore_ ) store_.restore_field( *f_ );
        delete count_;
			count_=NULL;
        delete f_;
        count_ = NULL;
        f_ = NULL;
      }
    }
  }
  //------------------------------------------------------------------
  
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>::~SpatFldPtr() { free(); }
  //------------------------------------------------------------------
  

  //==================================================================


  //------------------------------------------------------------------
  template<typename FieldT>
  SpatialFieldStore<FieldT>::~SpatialFieldStore()
  {
    for( typename FQMap::iterator ii=fqmap_.begin(); ii!=fqmap_.end(); ++ii ){
      FieldQueue& q = ii->second;
      while( !q.empty() ){
        AtomicT* field = q.front();
        delete [] field;
			field=NULL;
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
    return get( f.window_with_ghost() );
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  SpatFldPtr<FieldT>
  SpatialFieldStore<FieldT>::get( const structured::MemoryWindow& window )
  {
#ifdef ENABLE_THREADS
    boost::mutex::scoped_lock lock( get_mutex() );
#endif
    // find the proper map
    const int npts = window.local_npts();
    FieldQueue& q = fqmap_[ npts ];

    AtomicT* fnew;
    if( q.empty() ){
      fnew = new AtomicT[ npts ];
    }
    else{
      fnew = q.front();
      q.pop();
    }

    return SpatFldPtr<FieldT>( new FieldT(window,fnew,structured::ExternalStorage), true );
  }
  //------------------------------------------------------------------
  template<typename FieldT>
  void
  SpatialFieldStore<FieldT>::restore_field( FieldT& field )
  {
#ifdef ENABLE_THREADS
    boost::mutex::scoped_lock lock( get_mutex() );
#endif
    const structured::MemoryWindow& w = field.window_with_ghost();
    FieldQueue& q = fqmap_[ w.local_npts() ];
    q.push( &field[0] );
  }
  //------------------------------------------------------------------


  //====================================================================



  // specialized for doubles masquerading as spatialfields
  template<>
  inline void
  SpatialFieldStore<double>::restore_field( double& d )
  {}

  template<>
  SpatFldPtr<double>
  inline SpatialFieldStore<double>::get( const double& d )
  {
#ifdef ENABLE_THREADS
    boost::mutex::scoped_lock lock( get_mutex() );
#endif
    return SpatFldPtr<double>( new double, true );
  }

  template<>
  inline SpatFldPtr<double>
  SpatialFieldStore<double>::get( const structured::MemoryWindow& w )
  {
#ifdef ENABLE_THREADS
    boost::mutex::scoped_lock lock( get_mutex() );
#endif
    return SpatFldPtr<double>( new double, true );
  }
  
} // namespace SpatialOps

#endif
