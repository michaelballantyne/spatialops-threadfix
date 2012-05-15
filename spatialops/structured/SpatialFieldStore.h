/*
 * Copyright (c) 2011 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef UT_SpatialFieldStore_h
#define UT_SpatialFieldStore_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/ExternalAllocators.h>
#include <spatialops/structured/IndexTriplet.h>

#include <stack>
#include <map>

#include <boost/type_traits.hpp>

#ifdef ENABLE_THREADS
# include <boost/thread/mutex.hpp>
#endif

namespace SpatialOps {

// forward declaration
class SpatialFieldStore;

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
 *  to the SpatialFieldStore.  You can create SpatFldPtr objects
 *  from an existing SpatialField, but the SpatFldPtr will not
 *  assume ownership of the memory.
 */
template<typename FieldT>
class SpatFldPtr {
public:

  /**
   *  @brief Construct a SpatFldPtr.
   *
   *  @param field The field to wrap.  This constructor should be
   *  used if you want to wrap an existing SpatialField for use as a
   *  SpatFldPtr.  Ownership of this pointer is transfered.
   */
  SpatFldPtr(FieldT* const field);

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
  SpatFldPtr(FieldT* const field, const bool builtFromStore);

  /**
   *  @brief Constructor for use from SpatialFieldStore only and ONLY for
   *  		masquerade fields. This is specifically so that we can wrap
   *  		or allocate a double on an external device while pretending
   *  		it is a full spatial field.
   *
   *  @param pointer to wrap
   *  @param do we build it
   *  @param where the memory is allocated
   *  @param device index for memory lookup
   */
  SpatFldPtr(FieldT* const field, const bool builtFromStore, const MemoryType mtype, const unsigned short int deviceIndex);

  ~SpatFldPtr();

  /** @brief Copy constructor */
  SpatFldPtr(const SpatFldPtr<FieldT>& p);

  /** @brief Skeletal constructor */
  SpatFldPtr();

  /** @brief Assignment operator */
  SpatFldPtr& operator=(const SpatFldPtr& p);

  /** @brief Assignment operator */
  SpatFldPtr& operator=(FieldT* const f);

  inline FieldT& operator*() {
    return *f_;
  }

  inline const FieldT& operator*() const {
    return *f_;
  }

  inline FieldT* operator->() {
    return f_;
  }

  inline const FieldT* operator->() const {
    return f_;
  }

  inline bool isnull() const {
    return f_ == NULL;
  }

  int count() const {
    return *count_;
  }

  bool built_from_store() const {
    return builtFromStore_;
  }

  // Devin: this gets around having to make a variety of structures to mimic
  //		  partial specializations in downstream classes... It may not be
  //		  ideal, but it works for now and is easily replaced

  // Wrap some spatial field calls to get around
  // problems when trying to call methods of de-referenced pointers to
  inline unsigned int   allocated_bytes() const;
  inline double*        field_values()    const;
  inline unsigned short device_index()    const;

  void detach();

private:
  FieldT* f_;
  int* count_;
  bool builtFromStore_;
  unsigned short deviceIndex_;
  MemoryType memType_;
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
 *    SpatFldPtr<FieldT> field = SpatialFieldStore::get<FieldT>( const FieldT& f )
 *  \endcode
 *
 *  to return a field with the same dimensions as the provided
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
class SpatialFieldStore {

  template<typename T>
  class Pool{
    typedef std::stack<T*> FieldQueue;
    typedef std::map<size_t,FieldQueue> FQSizeMap;

    FQSizeMap fqm_;
    size_t pad_;
    size_t highWater_;

    Pool();
    ~Pool();

  public:
    static Pool& self();
    inline T* get( const size_t n );
    inline void put(T*, size_t);
    inline size_t active() const;
    inline size_t total() const{ return highWater_; }
  };

  template<typename FT, typename IsPODT> struct ValTypeSelector;
  template<typename FT> struct ValTypeSelector<FT, boost::true_type> {
    typedef FT type;
  };
  template<typename FT> struct ValTypeSelector<FT, boost::false_type> {
    typedef typename FT::AtomicT type;
  };

public:

  /**
   *  @brief Obtain a temporary field.
   *
   *  @param w  A memory window describing the desired field dimensions
   *
   *  Note that you should not dereference the SpatFldPtr object to
   *  store a SpatialField reference.  Doing so can cause memory
   *  corruption.
   */
  template<typename FieldT>
  inline static SpatFldPtr<FieldT>
  get_from_window( const structured::MemoryWindow& window,
                   const MemoryType mtype = LOCAL_RAM,
                   const unsigned short int deviceIndex = 0 );

  /**
   *  @brief Obtain a temporary field.
   *
   *  @param f  A field to model this one after.
   *
   *  Note that you should not dereference the SpatFldPtr object to
   *  store a SpatialField reference.  Doing so can cause memory
   *  corruption.
   */
  template< typename FieldT, typename ProtoT >
  inline static SpatFldPtr<FieldT>
  get( const ProtoT& f,
       const MemoryType mtype = LOCAL_RAM,
       const unsigned short int deviceIndex = 0 )
  {
    using namespace structured;
    const MemoryWindow& ws = f.window_with_ghost();
    const MemoryWindow w( ws.glob_dim() + Subtract< typename FieldT::Location::BCExtra, typename ProtoT::Location::BCExtra >::result::int_vec() * ws.has_bc(),
                          ws.offset(),
                          ws.extent()   + Subtract< typename FieldT::Location::BCExtra, typename ProtoT::Location::BCExtra >::result::int_vec() * ws.has_bc(),
                          ws.has_bc(0), ws.has_bc(1), ws.has_bc(2) );
    return get_from_window<FieldT>(w,mtype,deviceIndex);
  }



  /**
   *  @brief Restores a field to the store for future use.
   *
   *  @param Field to be restored to the store
   *
   *  Note that this method is should only be called by SpatFldPtr
   *  objects.  Calling it anywhere else can result in memory corruption.
   */
  template<typename FieldT>
  inline static void restore_field(FieldT& f);

//  inline static size_t active(){ return Pool<double>::self().active(); }
//  inline static size_t total() { return Pool<double>::self().total(); }

private:

#ifdef ENABLE_THREADS
  /**
   *  Used to lock threads to prevent simultaneous access.
   */
  inline static boost::mutex& get_mutex() {static boost::mutex m; return m;}
#endif

};

//==================================================================

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//  Implementation
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//=================================================================


//TODO: This is a problem.... We really need to specialize the SpatialField class to handle singleton values efficiently
//              As it is now, there are a number of semantic problems with how we have to treat doubles and we're not providing
//              the proper meta information to effectively work around them
template<>
inline SpatFldPtr<double>::
SpatFldPtr( double* const field,
            const bool builtFromStore,
            const MemoryType mtype,
            const unsigned short int deviceIndex )
  : f_(field),
    count_(new int),
    builtFromStore_(builtFromStore),
    deviceIndex_(deviceIndex),
    memType_(mtype)
{
  if( mtype != LOCAL_RAM ){
    std::cout << "TRIED TO BUILD AN EXTRNAL DOUBLE -- NOT SUPPORTED YET\n";
  }
  *count_ = 1;
}
// specialized for doubles masquerading as spatial fields

template<>
inline unsigned int SpatFldPtr<double>::allocated_bytes() const{
  return sizeof(double);
}

template<>
inline double* SpatFldPtr<double>::field_values() const {
  return f_;
}

template<>
inline unsigned short SpatFldPtr<double>::device_index() const {
  return deviceIndex_;
}

//====================================================================

template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr( FieldT* const f )
  : f_(f),
    count_(new int),
    builtFromStore_(false),
    memType_( ( f != NULL ? f->memory_device_type() : LOCAL_RAM ) ),
    deviceIndex_( ( f != NULL ? f->device_index() : 0 ) )
{
  *count_ = 1;
}

template<>
inline SpatFldPtr<double>::
SpatFldPtr( double* const f, const bool builtFromStore )
  : f_(f),
    count_(new int),
    builtFromStore_(builtFromStore),
    memType_(LOCAL_RAM),
    deviceIndex_(0)
{
  *count_ = 1;
}

//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr(FieldT* const f, const bool builtFromStore)
  : f_(f),
    count_(new int),
    builtFromStore_(builtFromStore),
    memType_( ( f != NULL ? f->memory_device_type() : LOCAL_RAM ) ),
    deviceIndex_( ( f != NULL ? f->device_index() : 0 ) )
{
  *count_ = 1;
}

//------------------------------------------------------------------

template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr(const SpatFldPtr<FieldT>& p)
  : f_(p.f_),
    count_(p.count_),
    builtFromStore_( p.builtFromStore_ ),
    deviceIndex_( p.deviceIndex_),
    memType_( p.memType_ )
{
    ++(*count_);
}

//------------------------------------------------------------------

template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr()
  : f_( NULL ),
    count_( NULL ),
    builtFromStore_( false ),
    deviceIndex_( 0 ),
    memType_( LOCAL_RAM )
{}
//------------------------------------------------------------------

template<typename FieldT>
SpatFldPtr<FieldT>&
SpatFldPtr<FieldT>::operator=(const SpatFldPtr& p)
{
  // was this an active SpatFldPtr?
  if (count_ != NULL) {
    // this one is dying so decrement the count.
    --(*count_);
    // kill the old one if needed
    if (*count_ == 0) {
      if (builtFromStore_) {
        switch (memType_) {
        case LOCAL_RAM: {
          SpatialFieldStore::restore_field(*f_);
        }
        break;

        case EXTERNAL_CUDA_GPU:
          //Do nothing, we don't store cuda memory here
          break;

        default:
          throw(std::runtime_error("Attempt to detach an unknown field type."));
        }
      }
      delete f_;
      delete count_;
      count_ = NULL;
      f_ = NULL;
    }
  }
  // reassign
  f_ = p.f_;
  count_ = p.count_;
  builtFromStore_ = p.builtFromStore_;
  memType_ = p.memType_;
  deviceIndex_ = p.deviceIndex_;
  // increment copy count
  ++(*count_);

  return *this;
}
//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT>&
SpatFldPtr<FieldT>::operator=(FieldT* const f) {
  // was this an active SpatFldPtr?
  if (count_ != NULL) {
    // this one is dying so decrement the count.
    --(*count_);
    // kill the old one if needed
    if (*count_ == 0) {
      if (builtFromStore_) {
        switch (memType_) {
        case LOCAL_RAM: {
          SpatialFieldStore::restore_field(*f_);
        }
        break;

        case EXTERNAL_CUDA_GPU:
          break;

        default:
          throw(std::runtime_error("Attempt to detach an unknown field type."));
        }
      }

      delete f_;
      delete count_;
      count_ = NULL;
      f_ = NULL;
    }
  }
  // reassign
  f_ = f;
  count_ = new int;
  *count_ = 1;
  builtFromStore_ = false;
  memType_ = f->memory_device_type();
  deviceIndex_ = f->device_index();

  return *this;
}

//------------------------------------------------------------------

template<typename FieldT>
void SpatFldPtr<FieldT>::detach() {
  // was this an active SpatFldPtr?
  if (count_ != NULL) {
    // this one is dying so decrement the count.
    --(*count_);
    if (*count_ == 0) {
      // kill the old one if needed
      if (builtFromStore_) {
        switch (memType_) {
        case LOCAL_RAM:
          SpatialFieldStore::restore_field(*f_);
          break;
        case EXTERNAL_CUDA_GPU:
          //Currently GPU fields are never built from the store
          break;
        default:
          throw(std::runtime_error("Attempt to detach an unknown field type."));
        }
      }
      delete count_;
      count_ = NULL;
      delete f_;
      count_ = NULL;
      f_ = NULL;
    }
  }
}

//------------------------------------------------------------------

template<typename FieldT>
inline unsigned int SpatFldPtr<FieldT>::allocated_bytes() const{
  return f_->allocated_bytes();
}

//------------------------------------------------------------------

template<typename FieldT>
inline double* SpatFldPtr<FieldT>::field_values() const {
  return f_->field_values();
}

//------------------------------------------------------------------

template<typename FieldT>
inline unsigned short SpatFldPtr<FieldT>::device_index() const {
  return f_->device_index();;
}

//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT>::~SpatFldPtr() {
  detach();
}
//------------------------------------------------------------------
//Wrap a double
template<>
inline SpatFldPtr<double>::SpatFldPtr(double* const f)
  : f_(f),
    count_(new int),
    builtFromStore_(false),
    memType_(LOCAL_RAM),
    deviceIndex_(0)
{
  *count_ = 1;
}

//==================================================================


// specialized for doubles masquerading as SpatialFields

template<>
inline
void
SpatialFieldStore::restore_field<double>( double& d )
{}

template<>
inline
SpatFldPtr<double>
SpatialFieldStore::
get<double,double>( const double& d,
                    const MemoryType mtype,
                    const unsigned short int deviceIndex )
{
# ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
# endif
  return SpatFldPtr<double>( new double, true, mtype, deviceIndex );
}

//====================================================================

//------------------------------------------------------------------

template<typename FieldT>
inline
SpatFldPtr<FieldT>
SpatialFieldStore::
get_from_window( const structured::MemoryWindow& window,
                 const MemoryType mtype,
                 const unsigned short int deviceIndex )
{
  typedef typename ValTypeSelector<FieldT,typename boost::is_pod<FieldT>::type>::type AtomicT;

#ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
#endif
  switch (mtype) {
  case LOCAL_RAM: { // Allocate from a store
    const structured::MemoryWindow mw( window.extent(),
                                       structured::IntVec(0,0,0),
                                       window.extent(),
                                       window.has_bc(0), window.has_bc(1), window.has_bc(2) );
    const size_t npts = mw.glob_npts();
    AtomicT* fnew = Pool<AtomicT>::self().get(npts);
#   ifndef NDEBUG  // only zero the field for debug runs.
    for( size_t i=0; i<npts; ++i )  fnew[i] = 0.0;
#   endif
    return SpatFldPtr<FieldT>( new FieldT(mw,fnew,structured::ExternalStorage), true );
  }
#ifdef ENABLE_CUDA
  //Dvn: I'm not having the store hold GPU memory right now, as I'm not sure it would be entirely stable
  // for single GPU systems ( we could end up holding all the memory and break the display functionality )
  // So they're all allocated as internal storage, but wrapped in a SpatFldPtr
  case EXTERNAL_CUDA_GPU: {
    return SpatFldPtr<FieldT>(
        new FieldT(window, NULL, structured::InternalStorage, mtype,
            deviceIndex ), true );
  }
#endif
  default: {
    std::ostringstream msg;
    msg << "Attempt to create Spatial Field Pointer wrapping ( "
        << DeviceTypeTools::get_memory_type_description(mtype)
    << " ) field type, without supporting libraries included\n";
    msg << "\t " << __FILE__ << " : " << __LINE__;
    throw(std::runtime_error(msg.str()));
  }
  }
}

//------------------------------------------------------------------

template<typename FieldT>
inline
void SpatialFieldStore::restore_field(FieldT& field)
{
# ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
# endif
  typedef typename ValTypeSelector<FieldT,typename boost::is_pod<FieldT>::type>::type AtomicT;
  Pool<AtomicT>::self().put( field.field_values(),
                             field.window_with_ghost().glob_npts() );
}

//------------------------------------------------------------------

template< typename T >
SpatialFieldStore::Pool<T>::Pool()
{
  pad_ = 0;
  highWater_ = 0;
}

template< typename T >
SpatialFieldStore::Pool<T>::~Pool()
{
  for( typename FQSizeMap::iterator i=fqm_.begin(); i!=fqm_.end(); ++i ){
    FieldQueue& fq = i->second;
    while( !fq.empty() ){
      delete [] fq.top();
      fq.pop();
    }
  }
}

template< typename T >
SpatialFieldStore::Pool<T>&
SpatialFieldStore::Pool<T>::self()
{
  static Pool p;
  return p;
}

template< typename T >
T*
SpatialFieldStore::Pool<T>::get( const size_t n )
{
  if( pad_==0 ) pad_ = n+n/10;
  typename FQSizeMap::iterator ifq = fqm_.lower_bound( n );
  if( ifq == fqm_.end() ){
    ifq = fqm_.insert( ifq, make_pair(n+pad_,FieldQueue()) );
  }
  T* field = NULL;
  FieldQueue& fq = ifq->second;
  if( fq.empty() ){
    ++highWater_;
    field = new T[n+pad_];
  }
  else{
    field = fq.top(); fq.pop();
  }
  return field;
}

template< typename T >
void
SpatialFieldStore::Pool<T>::put( T* t, size_t n )
{
  const typename FQSizeMap::iterator ifq = fqm_.lower_bound( n );
  assert( ifq != fqm_.end() );
  FieldQueue& fq = ifq->second;
  fq.push( t );
}

template<typename T>
size_t
SpatialFieldStore::Pool<T>::active() const{
  size_t n=0;
  for( typename FQSizeMap::const_iterator ifq=fqm_.begin(); ifq!=fqm_.end(); ++ifq ){
    n += ifq->second.size();
  }
  return highWater_-n;
}

} // namespace SpatialOps

#endif
