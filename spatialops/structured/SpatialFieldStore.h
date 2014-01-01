/*
 * Copyright (c) 2014 The University of Utah
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
#include <spatialops/structured/MemoryPool.h>
#include <spatialops/structured/GhostData.h>

#include <stack>
#include <map>

#include <boost/type_traits.hpp>

#ifdef ENABLE_THREADS
# include <boost/thread/mutex.hpp>
#endif

namespace SpatialOps {


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

  /**
   * Dissociate the field that this SpatFldPtr points to from this object,
   * potentially releasing the memory that it points to as well.
   */
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


public:

  /**
   *  @brief Obtain a temporary field.
   *
   *  @param w  A memory window describing the desired field dimensions
   *  @param bc the information on boundaries
   *  @param ghost ghost information
   *  @param mtype where this field should be located
   *  @param deviceIndex for multiple GPUs, which one to locate the field on
   *
   *  Note that you should not dereference the SpatFldPtr object to
   *  store a SpatialField reference.  Doing so can cause memory
   *  corruption.
   */
  template<typename FieldT>
  inline static SpatFldPtr<FieldT>
  get_from_window( const structured::MemoryWindow& window,
                   const structured::BoundaryCellInfo& bc,
                   const structured::GhostData& ghost,
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
       MemoryType mtype = UNKNOWN,
       short int deviceIndex = -9999 )
  {
    using namespace structured;

    if( deviceIndex == -9999 ) deviceIndex = f.device_index();

    if( mtype == UNKNOWN ) mtype = f.memory_device_type();

    const MemoryWindow& ws = f.window_with_ghost();
    GhostData gs = f.get_ghost_data();

    const BoundaryCellInfo bc = BoundaryCellInfo::build<FieldT>( f.boundary_info().has_bc() );

    const IntVec inc = bc.has_bc() * Subtract< typename FieldT::Location::BCExtra, typename ProtoT::Location::BCExtra >::result::int_vec();
    const MemoryWindow w( ws.glob_dim() + inc,
                          ws.offset(),
                          ws.extent()   + inc );

    return get_from_window<FieldT>(w,bc,gs,mtype,deviceIndex);
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
  inline static void restore_field(const MemoryType mtype, FieldT& f);

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

template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr( FieldT* const f )
  : f_(f),
    count_(new int),
    builtFromStore_(false),
    deviceIndex_( ( f != NULL ? f->device_index() : 0 ) ),
    memType_( ( f != NULL ? f->memory_device_type() : LOCAL_RAM ) )
{
  *count_ = 1;
}

//------------------------------------------------------------------

template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr(FieldT* const f, const bool builtFromStore)
  : f_(f),
    count_(new int),
    builtFromStore_(builtFromStore),
    deviceIndex_( ( f != NULL ? f->device_index() : 0 ) ),
    memType_( ( f != NULL ? f->memory_device_type() : LOCAL_RAM ) )
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
	SpatialFieldStore::restore_field(memType_, *f_);
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
	SpatialFieldStore::restore_field(memType_, *f_);
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
	SpatialFieldStore::restore_field(memType_, *f_);
      }
      delete count_;
      count_ = NULL;
      delete f_;
      f_ = NULL;
    }
  }
}

//------------------------------------------------------------------

template<typename FieldT>
SpatFldPtr<FieldT>::~SpatFldPtr() {
  detach();
}

//------------------------------------------------------------------

//====================================================================

//------------------------------------------------------------------

template<typename FieldT>
inline
SpatFldPtr<FieldT>
SpatialFieldStore::
get_from_window( const structured::MemoryWindow& window,
                 const structured::BoundaryCellInfo& bc,
                 const structured::GhostData& ghost,
                 const MemoryType mtype,
                 const unsigned short int deviceIndex )
{
  typedef typename FieldT::value_type ValT;
# ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
# endif
    const structured::MemoryWindow mw( window.extent(),
                                       structured::IntVec(0,0,0),
                                       window.extent() );
    const size_t npts = mw.glob_npts();

  switch (mtype) {
  case LOCAL_RAM: { // Allocate from a store
    ValT* fnew = structured::Pool<ValT>::self().get(mtype,npts);
#   ifndef NDEBUG  // only zero the field for debug runs.
    ValT* iftmp = fnew;
    for( size_t i=0; i<npts; ++i, ++iftmp )  *iftmp = 0.0;
#   endif
    return SpatFldPtr<FieldT>( new FieldT(mw,bc,ghost,fnew,structured::ExternalStorage), true );
  }
# ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    ValT* fnew = structured::Pool<ValT>::self().get(mtype, npts);
    return SpatFldPtr<FieldT>( new FieldT(window, bc, ghost, fnew,
                                          structured::ExternalStorage,
                                          mtype, deviceIndex ),
                               true );
  }
# endif
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
void SpatialFieldStore::restore_field( const MemoryType mtype, FieldT& field )
{
# ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
# endif
  typedef typename FieldT::value_type ValT;
  ValT * values = const_cast<ValT *>((const_cast<FieldT const &>(field)).field_values(field.memory_device_type(), field.device_index()));
  structured::Pool<ValT>::self().put( mtype, values );
}

} // namespace SpatialOps

#endif
