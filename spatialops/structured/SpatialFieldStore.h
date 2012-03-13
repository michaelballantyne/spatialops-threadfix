#ifndef UT_SpatialFieldStore_h
#define UT_SpatialFieldStore_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/ExternalAllocators.h>
#include <spatialops/structured/IndexTriplet.h>

#include <queue>
#include <map>
#include <set>

#include <boost/type_traits.hpp>

#ifdef ENABLE_THREADS
# include <boost/thread/mutex.hpp>
#endif

//TODO: GPU Single element field implementation -- consume only --
namespace SpatialOps {

// forward declaration
template<typename T> class SpatialFieldStore;

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

    inline SpatFldPtr& operator =(const double x) {
      *f_ = x;
      return *this;
    } ///< Assign this field to a constant

    /**
     *  @name unary operators
     *  these simply call through to the corresponding SpatialField unary operators.
     *  WARNING: (DEPRECATED) These should not be used going forward.
     */
    //@{
    inline SpatFldPtr& operator+=(const SpatFldPtr& p) {
      *f_ += *p;
      return *this;
    } ///< Add a SpatFldPtr to this.
    inline SpatFldPtr& operator-=(const SpatFldPtr& p) {
      *f_ -= *p;
      return *this;
    } ///< Subtract a SpatFldPtr from this.
    inline SpatFldPtr& operator*=(const SpatFldPtr& p) {
      *f_ *= *p;
      return *this;
    } ///< Multiply this by a SpatFldPtr
    inline SpatFldPtr& operator/=(const SpatFldPtr& p) {
      *f_ /= *p;
      return *this;
    } ///< Divide this by a SpatFldPtr

    inline SpatFldPtr& operator+=(const FieldT& p) {
      *f_ += p;
      return *this;
    } ///< Add a FieldT to this.
    inline SpatFldPtr& operator-=(const FieldT& p) {
      *f_ -= p;
      return *this;
    } ///< Subtract a FieldT from this.
    inline SpatFldPtr& operator*=(const FieldT& p) {
      *f_ *= p;
      return *this;
    } ///< Multiply this by a FieldT
    inline SpatFldPtr& operator/=(const FieldT& p) {
      *f_ /= p;
      return *this;
    } ///< Divide this by a FieldT

    inline SpatFldPtr& operator+=(const double x) {
      *f_ += x;
      return *this;
    } ///< Add a constant to this field
    inline SpatFldPtr& operator-=(const double x) {
      *f_ -= x;
      return *this;
    } ///< Subtract a constant from this field
    inline SpatFldPtr& operator*=(const double x) {
      *f_ *= x;
      return *this;
    } ///< Multiply this field by a constant
    inline SpatFldPtr& operator/=(const double x) {
      *f_ /= x;
      return *this;
    } ///< Divide this field by a constant
    //@}

    int getCount() {
      return *count_;
    }

    bool getBFS() {
      return builtFromStore_;
    }

    // Devin: this gets around having to make a variety of structures to mimic
    //		  partial specializations in downstream classes... It may not be
    //		  ideal, but it works for now and is easily replaced
    // Wrap some spatial field calls to get around
    // problems when trying to call methods of de-referenced pointers to
    inline unsigned int allocated_bytes() 	const;
    inline double* field_values() 			const;
    inline unsigned short device_index() 	const;

    void detach();

  private:
    SpatialFieldStore<FieldT>& store_;
    MemoryType memType_;
    unsigned short deviceIndex_;
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
template<typename FieldT>
class SpatialFieldStore {
    friend class SpatFldPtr<FieldT> ;
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
    inline SpatFldPtr<FieldT> get( const FieldT& f,
    		const MemoryType mtype = LOCAL_RAM, const unsigned short int deviceIndex = 0 );

    inline SpatFldPtr<FieldT> get(const structured::MemoryWindow& window,
         const MemoryType mtype = LOCAL_RAM, const unsigned short int deviceIndex = 0 );

    template< typename ProtoT >
    inline SpatFldPtr<FieldT> get( const ProtoT& f )
    {
      using namespace structured;
      const MemoryWindow& ws = f.window_with_ghost();
      const MemoryWindow w( ws.glob_dim() + Subtract< typename FieldT::Location::BCExtra, typename ProtoT::Location::BCExtra >::result::int_vec() * ws.has_bc(),
                            ws.offset(),
                            ws.extent()   + Subtract< typename FieldT::Location::BCExtra, typename ProtoT::Location::BCExtra >::result::int_vec() * ws.has_bc(),
                            ws.has_bc(0), ws.has_bc(1), ws.has_bc(2) );
      return get(w);
    }


  private:

    /**
     *  @brief Restores a field to the store for future use.
     *
     *  @param Field to be restored to the store
     *
     *  Note that this method is private to ensure it is only called
     *  by SpatFldPtr objects.  Calling it anywhere else can result in
     *  memory corruption.
     */
    inline void restore_field(FieldT& f);

#ifdef ENABLE_THREADS
    /**
     *  Used to lock threads to prevent simultaneous access.
     */
    inline boost::mutex& get_mutex() {static boost::mutex m; return m;}
#endif

    SpatialFieldStore() {};
    ~SpatialFieldStore();

    template<typename FT, typename IsPODT> struct ValTypeSelector;
    template<typename FT> struct ValTypeSelector<FT, boost::true_type> {
        typedef FT type;
    };
    template<typename FT> struct ValTypeSelector<FT, boost::false_type> {
        typedef typename FT::AtomicT type;
    };

    typedef typename ValTypeSelector<FieldT,
        typename boost::is_pod<FieldT>::type>::type AtomicT;
    typedef std::queue<AtomicT*> FieldQueue;
    typedef std::map<int, FieldQueue> FQMap;

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
SpatFldPtr<FieldT>::SpatFldPtr(FieldT* const f) :
    store_(
    		SpatialFieldStore<FieldT>::self()), f_(f),
    		count_(new int), builtFromStore_(false),
    		memType_( ( f != NULL ? f->memory_device_type() : LOCAL_RAM ) ),
    		deviceIndex_( ( f != NULL ? f->device_index() : 0 ) ) {

		*count_ = 1;
}

template<>
inline SpatFldPtr<double>::SpatFldPtr(double* const f, const bool builtFromStore) :
    store_(
    	SpatialFieldStore<double>::self()), f_(f),
    	count_(new int), builtFromStore_(builtFromStore),
    	memType_(LOCAL_RAM), deviceIndex_(0) {

		*count_ = 1;
}

//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr(FieldT* const f, const bool builtFromStore) :
    store_(
    	SpatialFieldStore<FieldT>::self()), f_(f),
    	count_(new int), builtFromStore_(
        builtFromStore),
        memType_( ( f != NULL ? f->memory_device_type() : LOCAL_RAM ) ),
        deviceIndex_( ( f != NULL ? f->device_index() : 0 ) ) {

		*count_ = 1;
}

//------------------------------------------------------------------

template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr(const SpatFldPtr<FieldT>& p) :
    store_(SpatialFieldStore<FieldT>::self()) {
  f_ = p.f_;
  count_ = p.count_;
  ++(*count_);
  builtFromStore_ = p.builtFromStore_;
  memType_ = p.memType_;
  deviceIndex_ = p.deviceIndex_;
}

//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT>::SpatFldPtr() :
    store_(SpatialFieldStore<FieldT>::self()) {
  f_ = NULL;
  count_ = NULL;
  deviceIndex_ = 0;
  builtFromStore_ = false;
  memType_ = LOCAL_RAM;
}
//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT>&
SpatFldPtr<FieldT>::operator=(const SpatFldPtr& p) {
  // was this an active SpatFldPtr?
  if (count_ != NULL) {
    // this one is dying so decrement the count.
    --(*count_);
    // kill the old one if needed
    if (*count_ == 0) {
      if (builtFromStore_) {
        switch (memType_) {
          case LOCAL_RAM: {
            store_.restore_field(*f_);
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
            store_.restore_field(*f_);
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
          case LOCAL_RAM: {
            store_.restore_field(*f_);
          }
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

//==================================================================

//------------------------------------------------------------------
template<typename FieldT>
SpatialFieldStore<FieldT>::~SpatialFieldStore() {
  for (typename FQMap::iterator ii = fqmap_.begin(); ii != fqmap_.end(); ++ii) {
    FieldQueue& q = ii->second;
    while (!q.empty()) {
      AtomicT* field = q.front();
      delete[] field;
      field = NULL;
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
SpatFldPtr<FieldT> SpatialFieldStore<FieldT>::get(const FieldT& f,
		const MemoryType mtype, const unsigned short int deviceIndex ) {
  // jcs note that we could create a window from the parent window
  // that was the minimum size and create the field based on that.
  // This could save a lot of memory in some cases.
  return get(f.window_with_ghost(), f.memory_device_type(), f.device_index());
}

//------------------------------------------------------------------
template<typename FieldT>
SpatFldPtr<FieldT> SpatialFieldStore<FieldT>::get(
    const structured::MemoryWindow& window, MemoryType mtype,
    unsigned short int deviceIndex) {
#ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
#endif
  // find the proper map

  switch (mtype) {
    case LOCAL_RAM: { // Allocate from a store
      const size_t npts = window.glob_npts();
      FieldQueue& q = fqmap_[npts];

      AtomicT* fnew;
      if (q.empty()) {
        fnew = new AtomicT[npts];
        for (size_t i = 0; i < npts; ++i)
          fnew[i] = 0.0;
      } else {
        fnew = q.front();
        for (size_t i = 0; i < npts; ++i)
          fnew[i] = 0.0;
        q.pop();
      }
      return SpatFldPtr<FieldT>(
          new FieldT(window, fnew, structured::ExternalStorage, mtype,
              deviceIndex), true);
    }
#ifdef ENABLE_CUDA
      //Dvn: I'm not having the store hold GPU memory right now, as I'm not sure it would be entirely stable
      // for single GPU systems ( we could end up holding all the memory and break the display functionality )
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
void SpatialFieldStore<FieldT>::restore_field(FieldT& field) {
#ifdef ENABLE_THREADS
  boost::mutex::scoped_lock lock( get_mutex() );
#endif
  const structured::MemoryWindow& w = field.window_with_ghost();
  FieldQueue& q = fqmap_[w.local_npts()];
  q.push(&field[0]);
}
//------------------------------------------------------------------



  //====================================================================


  // specialized for doubles masquerading as spatialfields
  template<>
  inline void
  SpatialFieldStore<double>::restore_field( double& d )
  {}

  template<>
  template<>
  SpatFldPtr<double>
  inline SpatialFieldStore<double>::get<double>( const double& d )
  {
#ifdef ENABLE_THREADS
    boost::mutex::scoped_lock lock( get_mutex() );
#endif

    return SpatFldPtr<double>( new double, true );
  }

//------------------------------------------------------------------
//====================================================================

//Wrap a double
template<>
inline SpatFldPtr<double>::SpatFldPtr(double* const f) :
    store_(
    	SpatialFieldStore<double>::self()), f_(f),
    	count_(new int), builtFromStore_(false),
    	memType_(LOCAL_RAM), deviceIndex_(0) {

		*count_ = 1;
}
//SpatFldPtr(FieldT* const field, const bool builtFromStore, const MemoryType mtype, const unsigned short int deviceIndex);

template<>
inline SpatFldPtr<double>::SpatFldPtr(double* const field, const bool builtFromStore,
		const MemoryType mtype, const unsigned short int deviceIndex) :
		store_(SpatialFieldStore<double>::self()), f_(field),
    	count_(new int), builtFromStore_(builtFromStore),
    	memType_(mtype), deviceIndex_(deviceIndex) {

		if( mtype != LOCAL_RAM ){
			std::cout << "WE TRIED TO BUILD AN EXTRNAL DOUBLE -- NOT SUPPORTED YET\n";
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

template<>
SpatFldPtr<double>
inline SpatialFieldStore<double>::get(const double& d, const MemoryType mtype, unsigned short deviceIndex) {
	#ifdef ENABLE_THREADS
	  boost::mutex::scoped_lock lock( get_mutex() );
	#endif

	return SpatFldPtr<double>(new double, true, mtype, deviceIndex);
}

template<>
inline SpatFldPtr<double> SpatialFieldStore<double>::get(
    const structured::MemoryWindow& w, MemoryType mtype,
    	unsigned short int deviceIndex ) {
	#ifdef ENABLE_THREADS
	  boost::mutex::scoped_lock lock( get_mutex() );
	#endif

	return SpatFldPtr<double>(new double, true, mtype, deviceIndex);
}

} // namespace SpatialOps

#endif
