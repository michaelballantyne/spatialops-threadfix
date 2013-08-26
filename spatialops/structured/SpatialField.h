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
 * ----------------------------------------------------------------------------
 * Available debugging flags:
 *
 * 	DEBUG_SF_ALL -- Enable all Spatial Field debugging flags.
 *
 */
//#define DEBUG_SF_ALL

#ifndef SpatialOps_SpatialField_h
#define SpatialOps_SpatialField_h

#include <iostream>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <map>
#include <algorithm>
#include <string.h> // for memcmp below...

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/structured/ExternalAllocators.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/GhostData.h>
#include <spatialops/structured/BoundaryCellInfo.h>
#include <spatialops/structured/MemoryPool.h>
#include <boost/static_assert.hpp>

#ifdef ENABLE_THREADS
#include <boost/thread/mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#endif

#ifdef ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace SpatialOps{
namespace structured{

  //Forward Declaration
  template <typename T>
    class Pool;

  enum StorageMode
  {
    InternalStorage,
    ExternalStorage
  };

  /**
   *  \class SpatialField
   *  \ingroup structured
   *  \ingroup fields
   *  \author James C. Sutherland
   *
   *  \brief Abstracts a field.
   *
   *  \tparam FieldLocation - type traits to describe the location of
   *    this field.  On staggered meshes, this will describe the mesh
   *    this field is associated with.  It also defines whether this
   *    field is on a volume or surface.
   *
   *  \tparam T - the underlying datatype (defaults to \c double)
   *
   *  \par Related classes:
   *   - \ref MemoryWindow
   *   - \ref SpatialFieldStore
   *
   *  \par Public Typedefs
   *   - \c field_type - this field's type
   *   - \c Location - the location type traits
   *   - \c value_type  - the type of underlying data being stored in this SpatialField
   *   - \c iterator, \c const_iterator - iterators to the elements in this field
   *   - \c interior_iterator, \c const_interior_iterator - iterators to the interior elements in this field (excludes ghost cells).
   */
  template< typename FieldLocation,
            typename GhostTraits,
            typename T=double >
  class SpatialField
  {
    typedef SpatialField<FieldLocation,GhostTraits,T> MyType;
    typedef std::map<unsigned short int, T*> ConsumerMap;

    MemoryWindow fieldWindow_;	        ///< Full representation of the window to the field ( includes ghost cells )
    MemoryWindow interiorFieldWindow_;  ///< Window representation sans ghost cells.

    const BoundaryCellInfo bcInfo_;     ///< information about this field's behavior on a boundary
    const GhostDataRT ghosts_;          ///< The total number of ghost cells on each face of this field.
    GhostDataRT validGhosts_;           ///< The number of valid ghost cells on each face of this field.

    T* fieldValues_;			///< Values associated with this field in the context of LOCAL_RAM
    T* fieldValuesExtDevice_;           ///< External field pointer ( This pointer will only be valid on the device it was created )
    const bool builtField_;		///< Indicates whether or not we created this field ( we could just be wrapping memory )

    MemoryType memType_; 		///< Indicates the type of device on which this field is allocated
    unsigned short deviceIndex_;        ///< Indicates which device is this field stored on

    bool readOnly_;                     ///< flag set to prevent write-access to the field
    bool disableInterior_;              ///< flag set to prevent interior iteration on this field

    //Note: Presently this is assumed to be a collection of GPU based < index, pointer pairs >;
    //      which is not as general is it likely should be, but GPUs are currently the only external
    //      device we're interested in supporting.
    ConsumerMap consumerFieldValues_;	///< Provides the ability to store and track copies of this field consumed on other devices.
    ConsumerMap myConsumerFieldValues_;	///< Provides the ability to correctly delete/release copies of this field that this field allocated

    unsigned long int allocatedBytes_;	///< Stores entire field size in bytes: sizeof(T) * glob.x * glob.y * glob.z

#   ifdef ENABLE_CUDA
    cudaStream_t cudaStream_;
#   endif

    inline void reset_values(const T* values);

#ifdef ENABLE_THREADS
    /**
     *  \class ExecMutex
     *  \brief Scoped lock.
     */
    class ExecMutex {
#   ifdef ENABLE_THREADS
      const boost::mutex::scoped_lock lock;
      inline boost::mutex& get_mutex() const {static boost::mutex m; return m;}

    public:
      ExecMutex() : lock( get_mutex() ) {}
      ~ExecMutex() {}
#   else
    public:
      ExecMutex() {
      }
      ~ExecMutex() {
      }
#   endif
    };
#endif

  public:

    typedef SpatialField<FieldLocation, GhostTraits, T> field_type;
    typedef GhostTraits Ghost;
    typedef FieldLocation Location;
    typedef T AtomicT;
    typedef T value_type;
    typedef MemoryWindow memory_window;
    typedef FieldIterator<field_type> iterator;
    typedef FieldIterator<field_type> interior_iterator;
    typedef ConstFieldIterator<field_type> const_iterator;
    typedef ConstFieldIterator<field_type> const_interior_iterator;

    /**
     *  \brief Construct a SpatialField
     *  \param window - the MemoryWindow that specifies this field
     *         including ghost cells.
     *  \param bc information on boundary treatment for this field
     *  \param ghosts information on ghost cells for this field
     *  \param fieldValues a pointer to memory to be wrapped by this field
     *  \param mode Storage options.  If InternalStorage then the
     *         fieldValues will be copied into an internal buffer.  If
     *         ExternalStorage then the fieldValues will be stored
     *         externally.  Efficiency suggests that ExternalStorage
     *         is best, since it will avoid excessive copies.  Safety
     *         suggests that InternalStorage is best, since it
     *         protects against memory corruption and inadvertent
     *         deletion of the field's underlying memory.
     *  \param consumerMemoryType describes where this field lives (e.g., CPU, GPU)
     *  \param devIdx the identifier for the GPU/accelerator if the field lives
     *         there. This allows for the case where multiple accelerators are
     *         on a given node.
     */
    SpatialField( const MemoryWindow& window,
                  const BoundaryCellInfo& bc,
                  const GhostDataRT& ghosts,
                  T* const fieldValues,
                  const StorageMode mode = InternalStorage,
                  const MemoryType consumerMemoryType = LOCAL_RAM,
                  const unsigned short int devIdx = 0 );

    /**
     *  \brief Shallow copy constructor.  This results in two fields
     *  that share the same underlying memory.
     */
    SpatialField(const SpatialField& other);

    /**
     *  \brief Shallow copy constructor with new window.
     */
    SpatialField(const MemoryWindow& window,
                 const SpatialField& other);

    virtual ~SpatialField();

    /**
     *  \brief Given the index for this field 0-based including
     *  ghosts, obtain a reference to the field value.
     *  WARNING: slow!
     *  NOTE: USEAGE IS DEPRECATED!! Not supported for external field types.
     */
    T& operator()(const size_t i, const size_t j, const size_t k);

    /**
     *  \brief Given the index for this field 0-based including
     *  ghosts, obtain a const reference to the field value.
     *  WARNING: slow!
     *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
     */
    const T& operator()(const size_t i, const size_t j, const size_t k) const;

    /**
     *  \brief Given the index for this field 0-based including
     *  ghosts, obtain a reference to the field value.
     *  WARNING: slow!
     *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
     */
    T& operator()(const IntVec& ijk);

    /**
     *  \brief Given the index for this field 0-based including
     *  ghosts, obtain a const reference to the field value.
     *  WARNING: slow!
     *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
     */
    const T& operator()(const IntVec& ijk) const;

    /**
     *  Index into this field (global index, 0-based in ghost cells).
     *  Note that if this field is windowed, this is still the global
     *  (not windowed) flat index.
     *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
     */
    inline T& operator[](const size_t i);
    inline const T& operator[](const size_t i) const;

    /**
     * \brief Iterator constructs for traversing memory windows.
     * Note: Iteration is not directly supported for external field types.
     * @return
     */
    inline const_iterator begin() const {
      if( memType_ != LOCAL_RAM && fieldValues_ == NULL ){
        std::ostringstream msg;
        msg << "Field type ( "
            << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
            << " does not support direct iteration, and has no local consumer field allocated.\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
      return const_iterator(fieldValues_, fieldWindow_);
    }

    inline iterator begin() {
      if( memType_ != LOCAL_RAM ){
        std::ostringstream msg;
        msg << "Field type ( "
            << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
            << " does not support direct iteration, and has no local consumer field allocated.\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
      if( readOnly_ ){
	std::ostringstream msg;
        msg << "Field type ( "
            << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
            << " is read-only, so it cannot support non-const iteration.\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
      return iterator(fieldValues_, fieldWindow_);
    }

    inline const_iterator end() const;
    inline iterator end();

    inline const_interior_iterator interior_begin() const {
      if( disableInterior_ ){
        std::ostringstream msg;
        msg << "Interior iterators cannot be obtained on resized fields" << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw( std::runtime_error(msg.str()) );
      }
      if( memType_ != LOCAL_RAM && fieldValues_ == NULL ){
        std::ostringstream msg;
        msg << "Field type ( "
            << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
            << " does not support direct iteration, and has no local consumer field allocated.\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
      return const_interior_iterator( fieldValues_, interiorFieldWindow_ );
    }

    inline interior_iterator interior_begin() {
      if( disableInterior_ ){
        std::ostringstream msg;
        msg << "Interior iterators cannot be obtained on resized fields" << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw( std::runtime_error(msg.str()) );
      }

      if (memType_ != LOCAL_RAM) {
        std::ostringstream msg;
        msg << "Field type ( "
            << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
            << " does not support non-const iteration.\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
      else if( readOnly_ ){
	std::ostringstream msg;
        msg << "Field type ( "
            << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
            << " has consumers, so it cannot support non-const iteration.\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
      return interior_iterator(fieldValues_, interiorFieldWindow_);
    }

    inline const_interior_iterator interior_end() const;
    inline interior_iterator interior_end();

    inline MyType& operator =(const MyType&);

    /**
     * @brief Comparison operators
     * WARNING: Slow in general and comparison with external fields will incur copy penalties.
     */
    bool operator!=(const MyType&) const;
    bool operator==(const MyType&) const;

    /**
     * @brief Make this field available on another device type, index pair. Adding consumer fields
     * 		  increases the memory held by the the spatial field by 'allocated_bytes' for each
     *		  unique device added.
     *
     *		  Note: This operation is guaranteed to be atomic
     *
     *		  Note: consumer fields are read-only and no functionality should ever depend on the
     *		  field being writable.
     *
     * @param mtype -- Device type where this field should be made available for consumption
     * @param deviceIndex -- Index to the proper device
     */
    void add_consumer(MemoryType consumerMemoryType, const unsigned short int consumerDeviceIndex);

    bool find_consumer(MemoryType consumerMemoryType, const unsigned short int consumerDeviceIndex) const;

    const BoundaryCellInfo& boundary_info() const{ return bcInfo_; }

    const MemoryWindow& window_without_ghost() const {
      return interiorFieldWindow_;
    }

    const MemoryWindow& window_with_ghost() const {
      return fieldWindow_;
    }

    /**
     *
     * @return Hardware device where this field is allocated
     */
    MemoryType memory_device_type() const {
      return memType_;
    }

    /**
     *
     * @return Index to the hardware device storing this field.
     */
    unsigned short int device_index() const {
      return deviceIndex_;
    }

#   ifdef ENABLE_CUDA
    void set_stream( const cudaStream_t& stream ) { cudaStream_ = stream; }

    cudaStream_t const & get_stream() const { return cudaStream_; }
#   endif

    /**
     * Field values will return a pointer to the field type, which is valid on the device and context supplied to
     * the function ( LOCAL_RAM, 0 ) by default.
     *
     * Note: If the desired field is intended to be a consumer field, then it must have been added previously via
     * a call to add_consumer(mtype, index)
     *
     * @param mtype -- Select the type of device we want a pointer to
     * @param deviceIndex -- Index of the device
     * @return
     */
    T* field_values(const MemoryType consumerMemoryType = LOCAL_RAM, const unsigned short int consumerDeviceIndex = 0);
    const T* field_values(const MemoryType consumerMemoryType = LOCAL_RAM, const unsigned short int consumerDeviceIndex = 0) const;

    /**
     * @brief obtain the ghost information for this field
     */
    const GhostDataRT& get_ghost_data() const{ return ghosts_; }

    /**
     * @brief Obtain the information about the valid number of ghosts for this
     *        field.  Manipulation of fields through nebo assignment operators
     *        can lead to modification of the number of valid ghost cells.  This
     *        information is recorded here.
     */
    const GhostDataRT& get_valid_ghost_data() const{ return validGhosts_; }

    unsigned int allocated_bytes() const {
    	return allocatedBytes_;
    }

    /**
     * @brief reset the active window given the provided number of valid ghosts
     * @param ghosts the number of ghosts, which cannot be larger than the number of valid ghosts on this field.
     *
     * This method should be used when a field assignment can only occur on a subset of the field due to invalidation of some of the ghost cells.
     */
    inline void reset_valid_ghosts( const GhostDataRT& ghosts ){
      const IntVec diff = validGhosts_.get_minus() - ghosts.get_minus();
      fieldWindow_ = MemoryWindow( fieldWindow_.glob_dim(),
                                   fieldWindow_.offset() + diff,
                                   fieldWindow_.extent() - diff - validGhosts_.get_plus() + ghosts.get_plus() );
      validGhosts_ = ghosts;
    }

    /**
     * @brief Obtain a child field that is reshaped.
     * @param extentModify the amount to modify the extent of the current field by
     * @param shift the number of grid points to shif the current field by
     * @return the reshaped child field
     *
     * The memory is the same as the parent field, but windowed differently.
     * Note that a reshaped field is considered read-only and you cannot obtain
     * interior iterators for these fields.
     */
    MyType
    inline reshape( const IntVec& extentModify,
                    const IntVec& shift ) const
    {
      MemoryWindow w( fieldWindow_.glob_dim(),
                      fieldWindow_.offset() + shift,
                      fieldWindow_.extent() + extentModify );
      return MyType( w, *this );
    }

    // jcs needs to be redone/removed
    template<typename NewGhost>
    inline field_type resize_ghost() const {
        typename GhostFromField<MyType>::result typedef OldGhost;

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * If one of the following Boost static asserts has failed,
         *  then the Nebo calculation you are trying cannot be
         *  performed because there are not enough ghost cells.
         *  Essentially, the chain of stencil operators you are
         *  trying requires more ghost cells than are present in
         *  the given fields.
         *
         * There are two possible solutions:
         * 1) Add more ghost cells to the fields.
         * 2) Change the calculation.
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
        BOOST_STATIC_ASSERT(int(NewGhost::nX) <= int(OldGhost::nX));
        BOOST_STATIC_ASSERT(int(NewGhost::pX) <= int(OldGhost::pX));
        BOOST_STATIC_ASSERT(int(NewGhost::bX) <= int(OldGhost::bX));
        BOOST_STATIC_ASSERT(int(NewGhost::nY) <= int(OldGhost::nY));
        BOOST_STATIC_ASSERT(int(NewGhost::pY) <= int(OldGhost::pY));
        BOOST_STATIC_ASSERT(int(NewGhost::bY) <= int(OldGhost::bY));
        BOOST_STATIC_ASSERT(int(NewGhost::nZ) <= int(OldGhost::nZ));
        BOOST_STATIC_ASSERT(int(NewGhost::pZ) <= int(OldGhost::pZ));
        BOOST_STATIC_ASSERT(int(NewGhost::bZ) <= int(OldGhost::bZ));

        // jcs should this use "ghosts_" or "validGhosts_" ????
        const GhostDataRT newGhosts( NewGhost::neg_int_vec(),
                                     NewGhost::pos_int_vec() );
        return MyType( window_with_ghost().resize_ghost(ghosts_,newGhosts), *this );
    }

    // jcs needs to be removed
    template<typename NewGhost>
    inline field_type resize_ghost_and_maintain_interior() const {
        typename GhostFromField<MyType>::result typedef OldGhost;
        typename MinimumGhostFromField<MyType>::result typedef Minimum;

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * If one of the following Boost static asserts has failed,
         *  then the Nebo calculation you are trying cannot be
         *  performed because there are not enough ghost cells.
         *  Essentially, the chain of stencil operators you are
         *  trying requires more ghost cells than are present in
         *  the given fields.
         *
         * There are two possible solutions:
         * 1) Add more ghost cells to the fields.
         * 2) Change the calculation.
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
        BOOST_STATIC_ASSERT(int(NewGhost::nX) <= int(OldGhost::nX));
        BOOST_STATIC_ASSERT(int(NewGhost::pX) <= int(OldGhost::pX));
        BOOST_STATIC_ASSERT(int(NewGhost::bX) <= int(OldGhost::bX));
        BOOST_STATIC_ASSERT(int(NewGhost::nY) <= int(OldGhost::nY));
        BOOST_STATIC_ASSERT(int(NewGhost::pY) <= int(OldGhost::pY));
        BOOST_STATIC_ASSERT(int(NewGhost::bY) <= int(OldGhost::bY));
        BOOST_STATIC_ASSERT(int(NewGhost::nZ) <= int(OldGhost::nZ));
        BOOST_STATIC_ASSERT(int(NewGhost::pZ) <= int(OldGhost::pZ));
        BOOST_STATIC_ASSERT(int(NewGhost::bZ) <= int(OldGhost::bZ));

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * If one of the following Boost static asserts has failed,
         *  then the Nebo calculation you are trying cannot be
         *  performed because in the presense of boundary conditions,
         *  the extra cells would not be populated with valid results.
         *
         * There are two possible solutions:
         * 1) Add more ghost cells to the fields.
         * 2) Change the calculation.
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
        BOOST_STATIC_ASSERT(int(NewGhost::nX) >= int(Minimum::nX));
        BOOST_STATIC_ASSERT(int(NewGhost::pX) >= int(Minimum::pX));
        BOOST_STATIC_ASSERT(int(NewGhost::bX) >= int(Minimum::bX));
        BOOST_STATIC_ASSERT(int(NewGhost::nY) >= int(Minimum::nY));
        BOOST_STATIC_ASSERT(int(NewGhost::pY) >= int(Minimum::pY));
        BOOST_STATIC_ASSERT(int(NewGhost::bY) >= int(Minimum::bY));
        BOOST_STATIC_ASSERT(int(NewGhost::nZ) >= int(Minimum::nZ));
        BOOST_STATIC_ASSERT(int(NewGhost::pZ) >= int(Minimum::pZ));
        BOOST_STATIC_ASSERT(int(NewGhost::bZ) >= int(Minimum::bZ));

        const GhostDataRT ghostNew( NewGhost::neg_int_vec(), NewGhost::pos_int_vec() );

//        IntVec plus = ghostNew.get_plus();
//        if( ghosts_.has_bc(0) ) plus[0] += MyType::Location::BCExtra::x_value();
//        if( ghosts_.has_bc(1) ) plus[1] += MyType::Location::BCExtra::y_value();
//        if( ghosts_.has_bc(2) ) plus[2] += MyType::Location::BCExtra::z_value();
//        ghostNew.set_plus( plus );

        // jcs should this use "ghosts_" or "validGhosts_" ????
        const MemoryWindow w = window_with_ghost().resize_ghost( ghosts_, ghostNew );
        return MyType( w, *this );
    }

    // jcs remove this
    template<typename Shift>
    inline field_type shift() const {
        typename GhostFromField<MyType>::result typedef OldGhost;

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * If one of the following Boost static asserts has failed,
         *  then the Nebo calculation you are trying cannot be
         *  performed because there are not enough ghost cells.
         *  Essentially, the chain of stencil operators you are
         *  trying requires more ghost cells than are present in
         *  the given fields.
         *
         * There are two possible solutions:
         * 1) Add more ghost cells to the fields.
         * 2) Change the calculation.
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
        BOOST_STATIC_ASSERT(Shift::X < 0 ? ((int)(Abs<Shift::X>::result) <= (int)(OldGhost::nX)) : true);
        BOOST_STATIC_ASSERT(Shift::X > 0 ? ((int)(Shift::X) <= (int)(OldGhost::pX)) : true);
        BOOST_STATIC_ASSERT(Shift::X > 0 ? ((int)(Shift::X) <= (int)(OldGhost::bX)) : true);
        BOOST_STATIC_ASSERT(Shift::Y < 0 ? ((int)(Abs<Shift::Y>::result) <= (int)(OldGhost::nY)) : true);
        BOOST_STATIC_ASSERT(Shift::Y > 0 ? ((int)(Shift::Y) <= (int)(OldGhost::pY)) : true);
        BOOST_STATIC_ASSERT(Shift::Y > 0 ? ((int)(Shift::Y) <= (int)(OldGhost::bY)) : true);
        BOOST_STATIC_ASSERT(Shift::Z < 0 ? ((int)(Abs<Shift::Z>::result) <= (int)(OldGhost::nZ)) : true);
        BOOST_STATIC_ASSERT(Shift::Z > 0 ? ((int)(Shift::Z) <= (int)(OldGhost::pZ)) : true);
        BOOST_STATIC_ASSERT(Shift::Z > 0 ? ((int)(Shift::Z) <= (int)(OldGhost::bZ)) : true);

        return MyType( window_with_ghost().shift(Shift::int_vec()), *this );
//        return MyType(window_with_ghost().template shift<Shift>(),
//		      *this);
    }

    // jcs remove this
    template<typename Shift>
    inline field_type shift_and_maintain_interior() const {

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * If one of the following Boost static asserts has failed,
         *  then somehow a shift is being performed on a field that
         *  is being written to.  That is, there is a shift on the
         *  lhs of an assignment.
         *
         * If one of these asserts has failed, you probably modified
         *  Nebo, and there is a side-effect that you did not
         *  anticipate.
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
        BOOST_STATIC_ASSERT(Shift::X == 0);
        BOOST_STATIC_ASSERT(Shift::Y == 0);
        BOOST_STATIC_ASSERT(Shift::Z == 0);

        return *this;
    }

    // jcs remove this
    template<typename NewGhost, typename Shift>
    inline field_type resize_ghost_and_shift() const {
        return resize_ghost<NewGhost>().template shift<Shift>();
    }

    // jcs remove this
    template<typename NewGhost, typename Shift>
    inline field_type resize_ghost_and_shift_and_maintain_interior() const {
        return resize_ghost_and_maintain_interior<NewGhost>().template shift_and_maintain_interior<Shift>();
    }
  };

//==================================================================
//
//                          Implementation
//
//==================================================================

template<typename Location, typename GhostTraits, typename T>
SpatialField<Location, GhostTraits, T>::
SpatialField( const MemoryWindow& window,
              const BoundaryCellInfo& bc,
              const GhostDataRT& ghost,
              T* const fieldValues,
              const StorageMode mode,
              const MemoryType mtype,
              const unsigned short int devIdx )
    : fieldWindow_(window),
      interiorFieldWindow_(window), // reset with correct info later
      bcInfo_( bc ),
      ghosts_     ( ghost ),
      validGhosts_( ghost ),
      fieldValues_( ( mtype == LOCAL_RAM ) ?
                    ( ( mode == ExternalStorage) ? fieldValues : (NULL) )
                    : ( NULL ) ),
      fieldValuesExtDevice_( (mtype == EXTERNAL_CUDA_GPU ) ?
    		  	  	  	  	 // Note: this assumes fieldValues is on the proper GPU....
                             ( ( mode == ExternalStorage ) ? fieldValues : (NULL) ) // reset gpu memory later
                             : ( NULL ) ),
      builtField_( mode == InternalStorage ),
      memType_( mtype ),
      deviceIndex_( devIdx ),
      readOnly_( false ),
      disableInterior_( false ),
      allocatedBytes_( 0 )
#     ifdef ENABLE_CUDA
      , cudaStream_( 0 )
#     endif
{ //InteriorStorage => we build a new field
  //Exterior storage => we wrap T*

# ifndef NDEBUG
  // ensure that we have a consistent BoundaryCellInfo object
  for( int i=0; i<3; ++i ){
    if( bcInfo_.has_bc(i) ){
      assert( bcInfo_.num_extra(i) == Location::BCExtra::int_vec()[i] );
    }
  }
# endif // NDEBUG

  // set the interior MemoryWindow
  IntVec ext = window.extent();
  IntVec ofs = window.offset();
  for (size_t i = 0; i < 3; ++i) {
    if (ext[i] > 1) {
      ext[i] -= ghosts_.get_minus(i) + ghosts_.get_plus(i);
      ofs[i] += ghosts_.get_minus(i);
    }
  }
  interiorFieldWindow_ = MemoryWindow( window.glob_dim(), ofs, ext );

  //Determine raw byte count -- this is sometimes required for external device allocation.
  allocatedBytes_ = sizeof(T) * ( window.glob_npts() );

  switch ( mtype ) {
      case LOCAL_RAM:
	if( mode == InternalStorage ){
	  try {
	    fieldValues_ = new T[window.glob_dim(0) * window.glob_dim(1) * window.glob_dim(2)];
	  }
	  catch(std::runtime_error& e){
	    std::cout << "Error occurred while allocating memory on LOCAL_RAM" << std::endl;
	    std::cout << e.what() << std::endl;
	    std::cout << __FILE__ << " : " << __LINE__ << std::endl;
	  }
	}
      break;
#ifdef ENABLE_CUDA
      case EXTERNAL_CUDA_GPU: {
        if( mode == InternalStorage ){
          // Allocate Memory, only if Storage Mode is INTERNAL.
          ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
          fieldValuesExtDevice_ = (T*)CDI.get_raw_pointer( allocatedBytes_, deviceIndex_ );
        }
        break;
      }
#endif
    default: {
      std::ostringstream msg;
      msg << "Unsupported attempt to create field of type ( "
          << DeviceTypeTools::get_memory_type_description(memType_)
          << " )\n" << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
  }

  if (mode == InternalStorage){
    reset_values(fieldValues);
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
SpatialField<Location, GhostTraits, T>::SpatialField( const SpatialField& other )
: fieldWindow_(other.fieldWindow_),
  interiorFieldWindow_(other.interiorFieldWindow_),
  bcInfo_( other.bcInfo_ ),
  ghosts_( other.ghosts_ ),
  validGhosts_( other.validGhosts_ ),
  fieldValues_(other.fieldValues_),
  fieldValuesExtDevice_(other.fieldValuesExtDevice_),
  builtField_(false),
  memType_(other.memType_),
  deviceIndex_(other.deviceIndex_),
  readOnly_( other.readOnly_ ),
  disableInterior_( other.disableInterior_ ),
  consumerFieldValues_(other.consumerFieldValues_),
  allocatedBytes_( other.allocatedBytes_ )
# ifdef ENABLE_CUDA
  , cudaStream_( other.cudaStream_ )
# endif
{}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
SpatialField<Location, GhostTraits, T>::
SpatialField( const MemoryWindow& window, const SpatialField& other )
: fieldWindow_(window),
  interiorFieldWindow_( other.interiorFieldWindow_ ), // This should not be used!
  bcInfo_( other.bcInfo_ ),
  ghosts_( other.ghosts_ ),
  validGhosts_( other.ghosts_ ),
  fieldValues_(other.fieldValues_),
  fieldValuesExtDevice_(other.fieldValuesExtDevice_),
  builtField_(false),
  memType_(other.memType_),
  deviceIndex_(other.deviceIndex_),
  readOnly_( other.readOnly_ ),
  disableInterior_( other.disableInterior_ ),
  consumerFieldValues_(other.consumerFieldValues_),
  allocatedBytes_( other.allocatedBytes_ )
# ifdef ENABLE_CUDA
    , cudaStream_( other.cudaStream_ )
# endif
{
  /*
   *  If the new window results in a view of the field that
   *  cuts into the interior, then disable interior iterators.
   */
  if( window.offset(0) > other.interiorFieldWindow_.offset(0) ||
      window.offset(1) > other.interiorFieldWindow_.offset(1) ||
      window.offset(2) > other.interiorFieldWindow_.offset(2) ||
      window.extent(0)+window.offset(0) < interiorFieldWindow_.offset(0)+interiorFieldWindow_.extent(0) ||
      window.extent(1)+window.offset(1) < interiorFieldWindow_.offset(1)+interiorFieldWindow_.extent(1) ||
      window.extent(2)+window.offset(2) < interiorFieldWindow_.offset(2)+interiorFieldWindow_.extent(2) )
  {
    // jcs I would like to be able to disable write access, but in multithreaded
    //     situations we must be able to write to subsets of the field.
    //readOnly_ = true;
    disableInterior_ = true;
  }

  // ensure that we are doing sane operations with the new window:
# ifndef NDEBUG
  assert( window.sanity_check() );

  const MemoryWindow& pWindow = other.window_with_ghost();
  for( size_t i=0; i<3; ++i ){
    assert( window.extent(i) + window.offset(i) <= pWindow.glob_dim(i) );
    assert( window.offset(i) < pWindow.glob_dim(i) );
  }
# endif
}

  //------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
SpatialField<Location, GhostTraits, T>::~SpatialField() {
#ifdef ENABLE_CUDA
  //Release any fields allocated for consumer use
  for( typename ConsumerMap::iterator i = myConsumerFieldValues_.begin(); i != myConsumerFieldValues_.end(); ++i ){
    Pool<T>::self().put( EXTERNAL_CUDA_GPU, i->second );
  }

  consumerFieldValues_.clear();
  myConsumerFieldValues_.clear();
#endif

  if ( builtField_ ) {
    switch ( memType_ ) {
    case LOCAL_RAM: {
      delete[] fieldValues_;
      fieldValues_ = NULL;
    }
    break;
#ifdef ENABLE_CUDA
    case EXTERNAL_CUDA_GPU: {

      //Deallocate local_ram consumer copy if it exists
      if( fieldValues_ != NULL ) {
        delete[] fieldValues_;
        fieldValues_ = NULL;
      }
      ema::cuda::CUDADeviceInterface::self().release( (void*)fieldValuesExtDevice_, deviceIndex_);
    }
    break;
#endif
    default:
      std::ostringstream msg;
      msg << "Attempt to release ( "
          << DeviceTypeTools::get_memory_type_description(memType_)
          << " ) field type, without supporting libraries -- this likely indicates a serious problem in field initialization\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw( std::runtime_error( msg.str() ) );
    }
  }
}

//------------------------------------------------------------------

template<typename FieldLocation, typename GhostTraits, typename T>
void SpatialField<FieldLocation, GhostTraits, T>::
reset_values( const T* values )
{
  switch ( memType_ ) {
  case LOCAL_RAM: {
    iterator ifld = begin();
    const iterator iflde = end();
    if ( values == NULL ) {
      for (; ifld != iflde; ++ifld)
        *ifld = 0.0;
    } else {
      for (; ifld != iflde; ++ifld, ++values)
        *ifld = *values;
    }
  }
  break;
#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();

    if( values == NULL ) {
      CDI.memset( fieldValuesExtDevice_, 0, allocatedBytes_, deviceIndex_ );
    } else {
      void* src = (void*)(values);
      CDI.memcpy_to( fieldValuesExtDevice_, src, allocatedBytes_, deviceIndex_ );
    }

  }
  break;
#endif

  default:
    std::ostringstream msg;
    msg << "Reset values called for unsupported field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " )"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
T* SpatialField<Location, GhostTraits, T>::
field_values( const MemoryType consumerMemoryType,
              const unsigned short int consumerDeviceIndex )
{
  if( readOnly_ ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support 'non-const' access to field values.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  switch( consumerMemoryType ){
  case LOCAL_RAM:{
    if( fieldValues_ == NULL ){
      std::ostringstream msg;
      msg << "Request for consumer field pointer on a device (Local RAM) for which it has not been allocated\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    return fieldValues_;
  }

#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    //Check local allocations first
    if( consumerMemoryType == memType_  &&  consumerDeviceIndex == deviceIndex_  ) {
      return fieldValuesExtDevice_;
    }

    typename ConsumerMap::const_iterator citer = consumerFieldValues_.find( consumerDeviceIndex );
    if( citer != consumerFieldValues_.end() ) {
      return ( citer->second );
    }

    std::ostringstream msg;
    msg << "Request for consumer field pointer on a device (GPU) for which it has not been allocated\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );

  }
#endif
  default:{
    std::ostringstream msg;
    msg << "Request for consumer field pointer to unknown or unsupported device\n";
    msg << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  }

}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
const T* SpatialField<Location, GhostTraits, T>::
field_values( const MemoryType consumerMemoryType,
    const unsigned short int consumerDeviceIndex ) const
{
  switch( consumerMemoryType ){
  case LOCAL_RAM:{
    if( fieldValues_ == NULL ){
      std::ostringstream msg;
      msg << "Request for consumer field pointer on a device (Local RAM) for which it has not been allocated\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    return fieldValues_;
  }

#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    //Check local allocations first
    if( consumerMemoryType == memType_  &&  consumerDeviceIndex == deviceIndex_  ) {
      return fieldValuesExtDevice_;
    }

    typename ConsumerMap::const_iterator citer = consumerFieldValues_.find( consumerDeviceIndex );
    if( citer != consumerFieldValues_.end() ) {
      return ( citer->second );
    }

    std::ostringstream msg;
    msg << "Request for consumer field pointer on a device (GPU) for which it has not been allocated\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );

  }
#endif
  default:{
    std::ostringstream msg;
    msg << "Request for consumer field pointer to unknown or unsupported device\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
void SpatialField<Location, GhostTraits, T>::
add_consumer( MemoryType consumerMemoryType,
              const unsigned short int consumerDeviceIndex )
{
  readOnly_ = true;  // once a consumer is added, we only allow read-only access
#ifdef DEBUG_SF_ALL
  std::cout << "Caught call to Spatial Field add_consumer for field : " << this->field_values() << "\n";
#endif
  //Check for local allocation
  if( consumerMemoryType == memType_ && consumerDeviceIndex == deviceIndex_ ) {
    return;
  }
#   ifdef ENABLE_THREADS
  //Make sure adding consumers is per-field atomic
  ExecMutex lock;
#   endif

  //Take action based on where the field must be available and where it currently is
  switch( consumerMemoryType ){
  case LOCAL_RAM: {
    switch( memType_ ) {
#ifdef ENABLE_CUDA
    case LOCAL_RAM: {
      // The only way we should get here is if for some reason a field was allocated as
      // LOCAL_RAM with a non-zero device index.
      // This shouldn't happen given how we define LOCAL_RAM at present, but I don't know if
      // we should consider it an error.
    }
    break;

    case EXTERNAL_CUDA_GPU: { // GPU field that needs to be available on the CPU
#ifdef DEBUG_SF_ALL
      std::cout << "EXTERNAL_CUDA_GPU field\n";
#endif

      if( fieldValues_ == NULL ) {  // Space is already allocated
#ifdef DEBUG_SF_ALL
        std::cout << "Consumer field does not exist, allocating...\n\n";
#endif
	fieldValues_ = Pool<T>::self().get( consumerMemoryType, ( allocatedBytes_/sizeof(T) ) );
#ifdef DEBUG_SF_ALL
        std::cout << "fieldValues_ == " << fieldValues_ << std::endl;
#endif
      }

#ifdef DEBUG_SF_ALL
      std::cout << "Calling memcpy with " << fieldValues_ << " " << fieldValuesExtDevice_
          << " " << allocatedBytes_ << " " << deviceIndex_ << std::endl;
#endif
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_from( fieldValues_, fieldValuesExtDevice_, allocatedBytes_, deviceIndex_ );
    }
    break;
#endif
    default:{
      std::ostringstream msg;
      msg << "Failed call to add_consumer on Spatial Field, unknown source device type\n"
          << "This error indicates a serious problem in how this field was originally created\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    }

  } // LOCAL_RAM
  break;

#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    switch( memType_ ) {
    case LOCAL_RAM: { //CPU Field needs to be available on a GPU
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();

      //Check to see if field memory exists
      if ( consumerFieldValues_.find( consumerDeviceIndex ) == consumerFieldValues_.end() ) {
        //Field doesn't exist, attempt to allocate it
        consumerFieldValues_[consumerDeviceIndex] = Pool<T>::self().get( consumerMemoryType, ( allocatedBytes_/sizeof(T) ) );
        myConsumerFieldValues_[consumerDeviceIndex] = consumerFieldValues_[consumerDeviceIndex];
      }

      CDI.memcpy_to( (void*)consumerFieldValues_[consumerDeviceIndex], fieldValues_, allocatedBytes_, consumerDeviceIndex );
    }
    break;

    case EXTERNAL_CUDA_GPU: {
      //GPU Field needs to be available on another GPU
      //Check to see if the field exists
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();

      if ( consumerFieldValues_.find( consumerDeviceIndex ) == consumerFieldValues_.end() ) {
        consumerFieldValues_[consumerDeviceIndex] = Pool<T>::self().get( consumerMemoryType, ( allocatedBytes_/sizeof(T) ) );
        myConsumerFieldValues_[consumerDeviceIndex] = consumerFieldValues_[consumerDeviceIndex];
      }

      CDI.memcpy_peer( (void*)consumerFieldValues_[consumerDeviceIndex],
          consumerDeviceIndex, fieldValuesExtDevice_, deviceIndex_, allocatedBytes_ );
    }
    break;

    default:{
      std::ostringstream msg;
      msg << "Failed call to add_consumer on Spatial Field, unknown source device type\n"
          << "This error indicates a serious problem in how this field was originally created\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    }
  } // EXTERNAL_CUDA_GPU
  break;
#endif

  default: {
    std::ostringstream msg;
    msg << "Failed call to add_consumer on Spatial Field, unknown destination device type\n"
        << "Ensure that you are compiling spatial ops with the proper end device support\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
bool SpatialField<Location, GhostTraits, T>::
find_consumer( MemoryType consumerMemoryType,
              const unsigned short int consumerDeviceIndex ) const
{
#ifdef DEBUG_SF_ALL
  std::cout << "Caught call to Spatial Field find_consumer for field : " << this->field_values() << "\n";
#endif
  //Check for local allocation
  if( consumerMemoryType == memType_ && consumerDeviceIndex == deviceIndex_ ) {
    return true;
  }
#   ifdef ENABLE_THREADS
  //Make sure adding consumers is per-field atomic
  ExecMutex lock;
#   endif

  //Take action based on where the field must be available and where it currently is
  switch( consumerMemoryType ){
  case LOCAL_RAM: {
    switch( memType_ ) {
#ifdef ENABLE_CUDA
    case EXTERNAL_CUDA_GPU: { // GPU field that needs to be available on the CPU
        return fieldValues_ != NULL;
    }
#endif
    default:{
      std::ostringstream msg;
      msg << "Failed call to find_consumer on Spatial Field, unknown source device type\n"
          << "This error indicates a serious problem in how this field was originally created\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    }

  } // LOCAL_RAM

#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
        return consumerFieldValues_.find( consumerDeviceIndex ) != consumerFieldValues_.end();
  } // EXTERNAL_CUDA_GPU
#endif

  default: {
    std::ostringstream msg;
    msg << "Failed call to find_consumer on Spatial Field, unknown destination device type\n"
        << "Ensure that you are compiling spatial ops with the proper end device support\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
typename SpatialField<Location, GhostTraits, T>::const_iterator
SpatialField<Location,GhostTraits, T>::end() const
{
  // We can allow constant iteration of the field even if its not local,
  // so long as it has a local consumer field allocated
  if( memType_ == LOCAL_RAM || fieldValues_ != NULL ){
    int extent = fieldWindow_.extent(0) * fieldWindow_.extent(1) * fieldWindow_.extent(2);
    const_iterator i(fieldValues_, fieldWindow_);
    return i + extent;
  } else {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__ << std::endl
        << "Unsupported request for const_iterator to field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " )"
        << "\t - No consumer allocated." << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
typename SpatialField<Location, GhostTraits, T>::iterator
SpatialField<Location,GhostTraits, T>::end()
{
  switch (memType_) {
  case LOCAL_RAM: {
    int extent = fieldWindow_.extent(0) * fieldWindow_.extent(1) * fieldWindow_.extent(2);
    iterator i(fieldValues_, fieldWindow_);
    return i + extent;
  }
  default:
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__ << std::endl
        << "Unsupported request for iterator to external field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " )"
        << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
typename SpatialField<Location, GhostTraits, T>::const_interior_iterator
SpatialField<Location,GhostTraits,T>::interior_end() const
{
  if( disableInterior_ ){
    std::ostringstream msg;
    msg << "Interior iterators cannot be obtained on resized fields" << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }
  if( memType_ == LOCAL_RAM || fieldValues_ != NULL ) {
    const size_t extent = interiorFieldWindow_.extent(0) * interiorFieldWindow_.extent(1) * interiorFieldWindow_.extent(2);
    const_interior_iterator i(fieldValues_, interiorFieldWindow_);
    return i + extent;
  } else {
    std::ostringstream msg;
    msg << "Unsupported request for iterator to external field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " )"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
typename SpatialField<Location, GhostTraits, T>::interior_iterator
SpatialField<Location,GhostTraits,T>::interior_end()
{
  if( disableInterior_ ){
    std::ostringstream msg;
    msg << "Interior iterators cannot be obtained on resized fields" << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }
  if( readOnly_ ){
    std::ostringstream msg;
    msg << "Cannot obtain a non-const iterator to a read-only field." << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }

  switch (memType_) {
  case LOCAL_RAM: {
    const size_t extent = interiorFieldWindow_.extent(0) * interiorFieldWindow_.extent(1) * interiorFieldWindow_.extent(2);
    interior_iterator i(fieldValues_, interiorFieldWindow_);
    return i + extent;
  }
  default:
    std::ostringstream msg;
    msg << "Unsupported request for iterator to external field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " )"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
T&
SpatialField<Location, GhostTraits, T>::
operator()( const size_t i, const size_t j, const size_t k )
{
  if( readOnly_ ){
    std::ostringstream msg;
    msg << "Cannot obtain a non-const iterator to a read-only field." << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }

  switch (memType_) {
  case LOCAL_RAM: {
    return (*this)(IntVec(i, j, k));
  }
  default:
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

template<typename Location, typename GhostTraits, typename T>
T&
SpatialField<Location, GhostTraits, T>::operator()(const IntVec& ijk)
{
  if( readOnly_ ){
    std::ostringstream msg;
    msg << "Cannot obtain a non-const iterator to a read-only field." << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }

  switch (memType_) {
  case LOCAL_RAM: {
#   ifndef NDEBUG
    assert(ijk[0] <  fieldWindow_.extent(0));
    assert(ijk[1] <  fieldWindow_.extent(1));
    assert(ijk[2] <  fieldWindow_.extent(2));
    assert(ijk[0] >= fieldWindow_.offset(0));
    assert(ijk[1] >= fieldWindow_.offset(1));
    assert(ijk[2] >= fieldWindow_.offset(2));
#   endif
    return fieldValues_[fieldWindow_.flat_index(ijk)];
  }
  default:
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
const T& SpatialField<Location, GhostTraits, T>::
operator()( const size_t i, const size_t j, const size_t k ) const
{
  if ( memType_ == LOCAL_RAM || fieldValues_ != NULL ) {
#   ifndef NDEBUG
    assert(i < fieldWindow_.extent(0));
    assert(j < fieldWindow_.extent(1));
    assert(k < fieldWindow_.extent(2));
    assert(i >= fieldWindow_.offset(0));
    assert(j >= fieldWindow_.offset(1));
    assert(k >= fieldWindow_.offset(2));
#   endif
    IntVec ijk(i,j,k);
    return fieldValues_[fieldWindow_.flat_index(ijk)];
  } else {
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support direct indexing, and has no consumer field allocated.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

template<typename Location, typename GhostTraits, typename T>
const T&
SpatialField<Location, GhostTraits, T>::
operator()( const IntVec& ijk ) const
{
  if( memType_ == LOCAL_RAM || fieldValues_ != NULL ){
#   ifndef NDEBUG
    assert( ijk[0] < fieldWindow_.extent(0) && ijk[0] >= fieldWindow_.offset(0) );
    assert( ijk[1] < fieldWindow_.extent(1) && ijk[1] >= fieldWindow_.offset(1) );
    assert( ijk[2] < fieldWindow_.extent(2) && ijk[2] >= fieldWindow_.offset(2) );
#   endif
    return fieldValues_[fieldWindow_.flat_index(ijk)];
  } else {
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support direct indexing, and has no consumer field allocated.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
T&
SpatialField<Location, GhostTraits, T>::operator[](const size_t i)
{
  if( readOnly_ ){
    std::ostringstream msg;
    msg << "Cannot obtain a non-const iterator to a read-only field." << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }

  switch (memType_) {
  case LOCAL_RAM: {
    return fieldValues_[i];
  }
  default:
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

//TODO: This function is not semantically consistent. What we want is
//		another [] overload which returns a constant reference and
//		functions in a consistent manner when a local consumer field
//		is allocated.
//		However, given the deprecated nature of the function, this may
//		not be an immediate issue.
template<typename Location, typename GhostTraits, typename T>
const T&
SpatialField<Location, GhostTraits, T>::operator[](const size_t i) const
{
  if( memType_ != LOCAL_RAM || fieldValues_ == NULL ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(memType_) << " ) ,"
        << " does not support direct indexing and has not allocated a consumer field.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  return fieldValues_[i];
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
SpatialField<Location, GhostTraits, T>&
SpatialField<Location, GhostTraits, T>::operator=(const MyType& other)
{
  if( readOnly_ ){
    std::ostringstream msg;
    msg << "Cannot obtain a non-const iterator to a read-only field." << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }

  if( allocatedBytes_ != other.allocatedBytes_ ) {
    std::ostringstream msg;
    msg << "Attempted assignment between fields of unequal size!\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  switch (memType_) {
  case LOCAL_RAM: {
    switch ( other.memory_device_type() ) {
    case LOCAL_RAM: { // LOCAL_RAM = LOCAL_RAM
      //Check for self assignment
      if( fieldValues_ == other.field_values() ){
        return *this;
      }
      std::copy( other.fieldValues_, other.fieldValues_+other.fieldWindow_.glob_npts(), fieldValues_ );
    }
    break;

#ifdef ENABLE_CUDA
    case EXTERNAL_CUDA_GPU: { //LOCAL_RAM = EXTERNAL_CUDA_GPU
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_from( fieldValues_, other.fieldValuesExtDevice_, allocatedBytes_, other.deviceIndex_ );
    }
    break;
#endif

    default:
      std::ostringstream msg;
      msg << "Attempted unsupported copy operation, at n\t"
          << __FILE__ << " : " << __LINE__ << std::endl
          << "\t - "
          << DeviceTypeTools::get_memory_type_description(memType_) << " = "
          << DeviceTypeTools::get_memory_type_description(
              other.memory_device_type()) << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    return *this;
  }
#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    switch( other.memory_device_type() ) {
    case LOCAL_RAM: {

      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_to( fieldValuesExtDevice_, other.fieldValues_, allocatedBytes_, deviceIndex_ );

    }
    break;

    case EXTERNAL_CUDA_GPU: {

      //Check for self assignment
      if( deviceIndex_ == other.deviceIndex_ && fieldValuesExtDevice_ == other.fieldValuesExtDevice_ ){
        return *this;
      }

      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_peer( fieldValuesExtDevice_, deviceIndex_, other.fieldValuesExtDevice_, other.deviceIndex_, allocatedBytes_ );

    }
    break;

    default: {
      std::ostringstream msg;
      msg << "Attempted unsupported copy operation, at " << std::endl
          << "\t" << __FILE__ << " : " << __LINE__ << std::endl
          << "\t - " << DeviceTypeTools::get_memory_type_description(memType_) << " = "
          << DeviceTypeTools::get_memory_type_description(other.memory_device_type()) << std::endl;
      throw( std::runtime_error ( msg.str() ));
    }

    } // end internal switch
    break;

  }
#endif
  default:
    std::ostringstream msg;
    msg << "Attempted unsupported copy operation, at \n\t"
        << __FILE__ << " : " << __LINE__ << std::endl
        << "\t - " << DeviceTypeTools::get_memory_type_description(memType_)
        << " = "
        << DeviceTypeTools::get_memory_type_description(
        other.memory_device_type()) << std::endl;
    throw(std::runtime_error(msg.str()));
  } // end outer switch
  return *this;
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
bool SpatialField<Location, GhostTraits, T>::operator!=(const MyType& other) const {
  return !(*this == other);
}

//------------------------------------------------------------------

template<typename Location, typename GhostTraits, typename T>
bool SpatialField<Location, GhostTraits, T>::operator==(const MyType& other) const
{
  switch (memType_) {
  case LOCAL_RAM: {
    switch (other.memory_device_type()) {
    case LOCAL_RAM: {
      const_iterator iother = other.begin();
      const_iterator iend = this->end();
      for (const_iterator ifld = this->begin(); ifld != iend; ++ifld, ++iother) {
        if (*ifld != *iother) return false;
      }
      return true;
    }
#ifdef ENABLE_CUDA
    case EXTERNAL_CUDA_GPU: {
      // Comparing LOCAL_RAM == EXTERNAL_CUDA_GPU
      // Note: This will incur a full copy penalty from the GPU and should not be used in a time sensitive context.
      if( allocatedBytes_ != other.allocatedBytes_ ) {
        throw( std::runtime_error( "Attempted comparison between fields of unequal size." ) );
      }

      void* temp = (void*)malloc(allocatedBytes_);
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_from( temp, other.fieldValuesExtDevice_, allocatedBytes_, other.deviceIndex_ );

      if( memcmp(temp, fieldValues_, allocatedBytes_) ) {
        free(temp);
        return false;
      }
      free(temp);
      return true;
    }
#endif
    default:{
      std::ostringstream msg;
      msg << "Attempted unsupported compare operation, at \n\t"
          << __FILE__ << " : " << __LINE__ << std::endl
          << "\t - "
          << DeviceTypeTools::get_memory_type_description(memType_) << " = "
          << DeviceTypeTools::get_memory_type_description(
              other.memory_device_type()) << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    } // switch( other.memory_device_type() )
  } // case LOCAL_RAM
#ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    switch( other.memory_device_type() ) {
    case LOCAL_RAM: {
      // Comparing EXTERNAL_CUDA_GPU == LOCAL_RAM
      // WARNING: This will incur a full copy penalty from the GPU and should not be used in a time sensitive context.
      if( allocatedBytes_ != other.allocatedBytes_ ) {
        throw( std::runtime_error( "Attempted comparison between fields of unequal size." ) );
      }

      void* temp = (void*)malloc(allocatedBytes_);
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_from( temp, fieldValuesExtDevice_, allocatedBytes_, deviceIndex_ );

      if( memcmp(temp, other.fieldValues_, allocatedBytes_) ) {
        free(temp);
        return false;
      }

      free(temp);
      return true;
    }

    case EXTERNAL_CUDA_GPU: {
      // Comparing EXTERNAL_CUDA_GPU == EXTERNAL_CUDA_GPU
      // WARNING: This will incur a full copy penalty from the GPU and should not be used in a time sensitive context.
      if( allocatedBytes_ != other.allocatedBytes_ ) {
        throw( std::runtime_error( "Attempted comparison between fields of unequal size." ) );
      }
      void* tempLHS = (void*)malloc(allocatedBytes_);
      void* tempRHS = (void*)malloc(allocatedBytes_);

      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_from( tempLHS, fieldValuesExtDevice_, allocatedBytes_, deviceIndex_ );
      CDI.memcpy_from( tempRHS, other.fieldValuesExtDevice_, allocatedBytes_, other.deviceIndex_ );

      free(tempLHS);
      free(tempRHS);

      if( memcmp(tempLHS, tempRHS, allocatedBytes_) ) {
        return false;
      }

      return true;
    }

    default: {
      std::ostringstream msg;
      msg << "Attempted unsupported compare operation, at \n\t"
          << __FILE__ << " : " << __LINE__ << std::endl
          << "\t - "
          << DeviceTypeTools::get_memory_type_description(memType_) << " = "
          << DeviceTypeTools::get_memory_type_description( other.memory_device_type() ) << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    }
  }
#endif
  default:
    std::ostringstream msg;
    msg << "Attempted unsupported compare operation, at \n\t"
        << __FILE__ << " : " << __LINE__ << std::endl
        << "\t - " << DeviceTypeTools::get_memory_type_description(memType_)
        << " = "
        << DeviceTypeTools::get_memory_type_description( other.memory_device_type()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_SpatialField_h
