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
#include <cmath>
#include <limits>

// Boost includes //
#include <boost/math/special_functions/next.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <string.h> // for memcmp below...

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/FVStaggeredLocationTypes.h>

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
  template <typename T> class Pool;

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
            typename T=double >
  class SpatialField
  {
  public:
    class SpatialFieldLoc {

    public:
      struct FieldInfo{
        T* field;
        bool isValid;
        bool _builtField;

        FieldInfo( T* field_, bool isValid_, bool builtField ){
          field = field_;
          isValid = isValid_;
          _builtField = builtField;
        }
        FieldInfo( T* field_, bool isValid_ )
        : field(field_), isValid(isValid_), _builtField(true)
        {}
        FieldInfo( T* field_ )
        : field(field_), isValid(true), _builtField(true)
        {}
        FieldInfo(){
          field = NULL;
          isValid = true;
          _builtField = true;
        }
      };

      typedef std::map<short int, FieldInfo> MultiFieldMap;
      MultiFieldMap multiFieldMap_;
      typedef typename MultiFieldMap::const_iterator cmapIter;
      typedef typename MultiFieldMap::iterator mapIter;

      MemoryWindow fieldWindow_;	        ///< Full representation of the window to the field ( includes ghost cells )
      MemoryWindow interiorFieldWindow_;  ///< Window representation sans ghost cells.

      const BoundaryCellInfo bcInfo_;     ///< information about this field's behavior on a boundary
      const GhostData ghosts_;          ///< The total number of ghost cells on each face of this field.
      GhostData validGhosts_;           ///< The number of valid ghost cells on each face of this field.

      T* fieldValues_;                    ///< Values associated with this field in the context of LOCAL_RAM
      T* fieldValuesExtDevice_;           ///< External field pointer ( This pointer will only be valid on the device it was created )

      const bool builtField_;             ///< Indicates whether or not we created this field ( we could just be wrapping memory )
      short int activeDeviceIndex_;
      short int actualDeviceIndex_;

      bool disableInterior_;              ///< flag set to prevent interior iteration on this field

      //Note: Presently this is assumed to be a collection of GPU based < index, pointer pairs >;
      //      which is not as general is it likely should be, but GPUs are currently the only external
      //      device we're interested in supporting.

      unsigned long int allocatedBytes_;	///< Stores entire field size in bytes: sizeof(T) * glob.x * glob.y * glob.z

  #   ifdef ENABLE_THREADS
      int partitionCount_; // Number of partitions Nebo uses in its thread-parallel backend when assigning to this field
  #   endif

  #   ifdef ENABLE_CUDA
      cudaStream_t cudaStream_;
  #   endif

      inline void reset_values(const T* values);

  #   ifdef ENABLE_THREADS
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
  #   endif

      typedef SpatialField<FieldLocation,T> field_type;
      typedef T value_type;
      typedef MemoryWindow memory_window;
      typedef FieldIterator<field_type> iterator;
      typedef FieldIterator<field_type> interior_iterator;
      typedef ConstFieldIterator<field_type> const_iterator;
      typedef ConstFieldIterator<field_type> const_interior_iterator;

      SpatialFieldLoc( const MemoryWindow& window,
                       const BoundaryCellInfo& bc,
                       const GhostData& ghosts,
                       T* const fieldValues,
                       const StorageMode mode = InternalStorage,
                       const short int devIdx = CPU_INDEX );

      SpatialFieldLoc(const SpatialFieldLoc& other);

      SpatialFieldLoc(const MemoryWindow& window,
                      const SpatialFieldLoc& other);

      ~SpatialFieldLoc();

      T& operator()(const size_t i, const size_t j, const size_t k);
      const T& operator()(const size_t i, const size_t j, const size_t k) const;
      T& operator()(const IntVec& ijk);
      const T& operator()(const IntVec& ijk) const;
      T& operator[](const size_t i);
      const T& operator[](const size_t i) const;

      inline const_iterator cbegin(MemoryWindow const window) const {
        cmapIter mapIter = multiFieldMap_.find( CPU_INDEX );
#       ifndef NDEBUG
        if( mapIter == multiFieldMap_.end() || mapIter->second.field == NULL ){
          std::ostringstream msg;
          msg << "Field type "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
              << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
        if( !mapIter->second.isValid ){
          std::ostringstream msg;
          msg << "Field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
              << " is not valid. const-iterators are not allowed on to a invalid field.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
#       endif
        return const_iterator(mapIter->second.field, window);
      }

      inline iterator begin( MemoryWindow const window )
      {
        mapIter mapIter = multiFieldMap_.find( CPU_INDEX );
#       ifndef NDEBUG
        if( mapIter == multiFieldMap_.end() ){
          std::ostringstream msg;
          msg << "Field type "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
              << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
        if( activeDeviceIndex_ != CPU_INDEX ) {
          std::ostringstream msg;
          msg << "Field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
              << " does not support direct iteration.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
#       endif
        return iterator( mapIter->second.field, window );
      }

      inline const_iterator cend( MemoryWindow const window ) const
      {
        cmapIter mapIter = multiFieldMap_.find( CPU_INDEX );
        if( mapIter == multiFieldMap_.end() || mapIter->second.field != NULL ){
          int extent = window.extent(0) * window.extent(1) * window.extent(2);
          const_iterator i(mapIter->second.field, window);
          return i + extent;
        }
        else {
          std::ostringstream msg;
          msg << "Unsupported request for const_iterator to field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " )"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw std::runtime_error( msg.str() );
        }
      }

      inline iterator end( MemoryWindow const window )
      {
        mapIter mapIter = multiFieldMap_.find( CPU_INDEX );
#        ifndef NDEBUG
        if( mapIter == multiFieldMap_.end() ){
          std::ostringstream msg;
          msg << "Field type "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
              << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
#       endif
        if( activeDeviceIndex_ == CPU_INDEX ){
          int extent = window.extent(0) * window.extent(1) * window.extent(2);
          iterator i(mapIter->second.field, window);
          return i + extent;
        }
        else{
          std::ostringstream msg;
          msg << __FILE__ << " : " << __LINE__ << std::endl
              << "Unsupported request for iterator to external field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " )"
              << std::endl;
          throw std::runtime_error( msg.str() );
        }
      }

      inline const_interior_iterator interior_begin() const
      {
        cmapIter const_mapIter = multiFieldMap_.find(CPU_INDEX);
#       ifndef NDEBUG
        if( disableInterior_ ){
          std::ostringstream msg;
          msg << "Interior iterators cannot be obtained on resized fields" << std::endl
              << __FILE__ << " : " << __LINE__ << std::endl;
          throw( std::runtime_error(msg.str()) );
        }
        if( const_mapIter == multiFieldMap_.end() || const_mapIter->second.field == NULL ){
          std::ostringstream msg;
          msg << "Field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
              << " doesn't exist in the map. Something wrong is with the memory.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
        if( !const_mapIter->second.isValid ){
          std::ostringstream msg;
          msg << "Field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
              << " is not valid. const_interior_iterator are not allowed on to a invalid field.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
#       endif
        return const_interior_iterator( const_mapIter->second.field, interiorFieldWindow_ );
      }

      inline interior_iterator interior_begin()
      {
        mapIter mapIter = multiFieldMap_.find( CPU_INDEX );
#       ifndef NDEBUG
        if( disableInterior_ ){
          std::ostringstream msg;
          msg << "Interior iterators cannot be obtained on resized fields" << std::endl
              << __FILE__ << " : " << __LINE__ << std::endl;
          throw( std::runtime_error(msg.str()) );
        }
        if( mapIter == multiFieldMap_.end() ){
          std::ostringstream msg;
          msg << "Field type "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
          << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
        if( activeDeviceIndex_ != CPU_INDEX ){
          std::ostringstream msg;
          msg << "Field type ( "
              << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
              << " does not support non-const iteration.\n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
#       endif
        return interior_iterator( mapIter->second.field, interiorFieldWindow_);
      }

      inline const_interior_iterator interior_end() const;
      inline interior_iterator interior_end();


  #   ifdef ENABLE_CUDA
      void add_field_loc(short int deviceLocation);

      void sync_location(short int deviceLocation);

      void remove_multiple_fields();
  #   endif

      void set_field_loc_active(short int deviceLocation);

      bool find_field_loc( const short int deviceLocation ) const;

      bool is_valid( const short int deviceLocation ) const;

      bool has_multiple_fields() const;

      const BoundaryCellInfo& boundary_info() const{ return bcInfo_; }

      const MemoryWindow& window_without_ghost() const {
        return interiorFieldWindow_;
      }

      const MemoryWindow& window_with_ghost() const {
        return fieldWindow_;
      }

      short int active_device_index() const {
        return activeDeviceIndex_;
      }

  #   ifdef ENABLE_THREADS
      void set_partition_count( const int count) { partitionCount_ = count; }

      int get_partition_count() { return partitionCount_; }
  #   endif

  #   ifdef ENABLE_CUDA
      void set_stream( const cudaStream_t& stream ) { cudaStream_ = stream; }

      cudaStream_t const & get_stream() const { return cudaStream_; }
  #   endif

      inline T* field_values(const short int deviceLocation = CPU_INDEX);
      inline const T* cfield_values(const short int deviceLocation = CPU_INDEX) const;

      const GhostData& get_ghost_data() const{ return ghosts_; }

      const GhostData& get_valid_ghost_data() const{ return validGhosts_; }

      unsigned int allocated_bytes() const {
          return allocatedBytes_;
      }

      inline void reset_valid_ghosts( const GhostData& ghosts ){
        const IntVec diff = validGhosts_.get_minus() - ghosts.get_minus();
        fieldWindow_ = MemoryWindow( fieldWindow_.glob_dim(),
                                     fieldWindow_.offset() + diff,
                                     fieldWindow_.extent() - diff - validGhosts_.get_plus() + ghosts.get_plus() );
        validGhosts_ = ghosts;
      }

      }; // SpatialFieldLocation



   typedef SpatialField<FieldLocation,T> field_type;
   typedef FieldLocation Location;
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

   inline SpatialField( const MemoryWindow& window,
                        const BoundaryCellInfo& bc,
                        const GhostData& ghosts,
                        T* const fieldValues,
                        const StorageMode mode = InternalStorage,
                        const short int devIdx = CPU_INDEX )
     : fieldWindow_(window),
       bcInfo_(bc.limit_by_extent(window.extent())),
       validGhosts_(ghosts.limit_by_extent(window.extent())),
       sfsharedPtr_(new SpatialFieldLoc( window, bc, ghosts, fieldValues, mode, devIdx ))
   {}

   /**
    *  \brief Shallow copy constructor.  This results in two fields
    *  that share the same underlying memory.
    */
   inline SpatialField(const SpatialField& other)
   : fieldWindow_(other.fieldWindow_),
     bcInfo_(other.bcInfo_),
     validGhosts_(other.validGhosts_),
     sfsharedPtr_(other.sfsharedPtr_)
   {}

   /**
    *  \brief Shallow copy constructor with new window.
    */
   inline SpatialField(const MemoryWindow& window,
                       const SpatialField& other)
   : fieldWindow_(window),
     bcInfo_(other.bcInfo_.limit_by_extent(window.extent())),
     validGhosts_( other.validGhosts_.limit_by_extent(window.extent())),
     sfsharedPtr_( other.sfsharedPtr_ )
   {
     //SpatialFieldLoc( window, *(other.sfsharedPtr_) );
   }

   virtual inline ~SpatialField();

   /**
    *  \brief Given the index for this field 0-based including
    *  ghosts, obtain a reference to the field value.
    *  WARNING: slow!
    *  NOTE: USEAGE IS DEPRECATED!! Not supported for external field types.
    */
   T& operator()(const size_t i, const size_t j, const size_t k){
     return sfsharedPtr_->operator()( i, j, k );
   }

   /**
    *  \brief Given the index for this field 0-based including
    *  ghosts, obtain a const reference to the field value.
    *  WARNING: slow!
    *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
    */
   const T& operator()(const size_t i, const size_t j, const size_t k) const{
     return sfsharedPtr_->operator()( i, j, k );
   }

   /**
    *  \brief Given the index for this field 0-based including
    *  ghosts, obtain a reference to the field value.
    *  WARNING: slow!
    *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
    */
   T& operator()(const IntVec& ijk){
     return sfsharedPtr_->operator()( ijk );
   }

   /**
    *  \brief Given the index for this field 0-based including
    *  ghosts, obtain a const reference to the field value.
    *  WARNING: slow!
    *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
    */
   const T& operator()(const IntVec& ijk) const{
     return sfsharedPtr_->operator()( ijk );
   }

   /**
    *  Index into this field (global index, 0-based in ghost cells).
    *  Note that if this field is windowed, this is still the global
    *  (not windowed) flat index.
    *  NOTE: USAGE IS DEPRECATED!! Not supported for external field types
    */
   T& operator[](const size_t i){
     return sfsharedPtr_->operator[]( i );
   }
   const T& operator[](const size_t i) const{
     return (const_cast<shared_ptrT&>(sfsharedPtr_))->operator[]( i );
   }

   /**
    * \brief Iterator constructs for traversing memory windows.
    * Note: Iteration is not directly supported for external field types.
    * @return
    */
   inline const_iterator begin() const {
     return sfsharedPtr_->cbegin(fieldWindow_);
   }
   inline iterator begin() {
     return sfsharedPtr_->begin(fieldWindow_);
   }

   inline const_iterator end() const{
     return sfsharedPtr_->cend(fieldWindow_);
   }
   inline iterator end(){
     return sfsharedPtr_->end(fieldWindow_);
   }

   inline const_interior_iterator interior_begin() const {
     return sfsharedPtr_->interior_begin();
   }

   inline interior_iterator interior_begin() {
     return sfsharedPtr_->interior_begin();
   }

   inline const_interior_iterator interior_end() const{
     return sfsharedPtr_->interior_end();
   }
   inline interior_iterator interior_end(){
     return sfsharedPtr_->interior_end();
   }

   inline field_type& operator =(const field_type&);

   /**
    * @brief Adds a location to a field. If the field location already exists,
    *        check if it is valid and make it as the active field Location.
    *        Else, use sync_location() and set the field as active.
    *
    *        increases the memory held by the the spatial field by 'allocated_bytes' for each
    *        unique device added.
    *        Note: This operation is guaranteed to be atomic
    *
    * @param deviceLocation -- Index to the proper device
    */

#   ifdef ENABLE_CUDA
   inline void add_field_loc(short int deviceLocation){
      sfsharedPtr_->add_field_loc( deviceLocation );
   };

  /**
   * @brief synchronizes all the invalid field location (deviceLocation)
   *        with the valid field location. Allocates memory and performs
   *        data-transfer whenever needed.
   *
   * @param deviceLocation -- Index to the proper device
   */
   inline void sync_location(short int deviceLocation){
     sfsharedPtr_->sync_location( deviceLocation );
   }

   /**
    * @brief Finds if the field has any other locations using has_multiple_fields().
    *        Frees all the multiple-field locations back to the relevant Memory Pool
    *        using remove_multiple_fields() expect the initial field.
    *
    */
   inline void remove_multiple_fields(){
     sfsharedPtr_->remove_multiple_fields();
   }
#   endif

  /**
   * @brief sets a field location as active. This method will update active device location so that
   *        writes can be performed.
   *
   * @param deviceLocation -- Device type where this field should be made as active
   */
   inline void set_field_loc_active(const short int deviceLocation){
     sfsharedPtr_->set_field_loc_active( deviceLocation );
   }

  /**
   * @brief queries for the availability of the field location which is mostly being used by Nebo
   *
   * @param deviceLocation -- Device type under query
   */
   inline bool find_field_loc( const short int deviceLocation ) const{
     return sfsharedPtr_->find_field_loc( deviceLocation );
   };

  /**
   * @brief checks if the field location is valid.
   *
   * @param deviceLocation -- Device type under query
   */
   bool is_valid( const short int deviceLocation ) const{
     return sfsharedPtr_->is_valid( deviceLocation );
   }

  /**
   * @brief reports if the spatial field has multiple field locations
   *
   */
   inline bool has_multiple_fields() const{
     return sfsharedPtr_->has_multiple_fields();
   }

   const BoundaryCellInfo& boundary_info() const{
     return sfsharedPtr_->bcInfo_;
   }

   const MemoryWindow& window_without_ghost() const {
     return sfsharedPtr_->interiorFieldWindow_;
   }

   const MemoryWindow& window_with_ghost() const {
     return fieldWindow_;
   }

  /**
   * @brief returns the current active device index.
   *
   */
   short int device_index() const {
     return sfsharedPtr_->active_device_index();
   }

#  ifdef ENABLE_THREADS
   /**
    * Sets number of partitions Nebo uses in its thread-parallel backend when assigning to this field
    *
    */
   void set_partition_count( const int count) { sfsharedPtr_->partitionCount_ = count; }

   /**
    * Returns number of partitions Nebo uses in its thread-parallel backend when assigning to this field
    *
    */
   int get_partition_count() { return sfsharedPtr_->partitionCount_; }
#   endif

#   ifdef ENABLE_CUDA
    /**
     * @brief sets the cuda stream to aid the kernel execution and async data transfer
     */
   void set_stream( const cudaStream_t& stream ) { sfsharedPtr_->cudaStream_ = stream; }

    /**
     * @brief returns the cuda stream for the SpatialField
     */
   cudaStream_t const & get_stream() const { return sfsharedPtr_->cudaStream_; }
#   endif

   /**
    * Field values will return a pointer to the field type, which is valid on the device and context supplied to
    * the function ( CPU_INDEX ) by default.
    *
    * Note: This method will invalidate all the other device locations apart from the deviceLocation.
    *
    * @param deviceLocation -- Index of the device
    * @return
    */
   inline T* field_values(const short int deviceLocation = CPU_INDEX){
     return sfsharedPtr_->field_values( deviceLocation );
   }

   /**
    * Field values will return a pointer to the field type, which is valid on the device and context supplied to
    * the function ( CPU_INDEX ) by default.
    *
    * Note: This method will perform a check if the deviceLocation is valid and not necessary to be active
    * field Location.
    *
    * @param deviceLocation -- Index of the device
    * @return
    */
   inline const T* field_values(const short int deviceLocation = CPU_INDEX) const{
     return sfsharedPtr_->cfield_values( deviceLocation );
   }

   /**
    * @brief obtain the ghost information for this field
    */
   const GhostData& get_ghost_data() const{ return sfsharedPtr_->ghosts_; }

   /**
    * @brief Obtain the information about the valid number of ghosts for this
    *        field.  Manipulation of fields through nebo assignment operators
    *        can lead to modification of the number of valid ghost cells.  This
    *        information is recorded here.
    */
   const GhostData& get_valid_ghost_data() const{ return validGhosts_; }

   unsigned int allocated_bytes() const{ return sfsharedPtr_->allocatedBytes_; }


   /**
    * @brief reset the active window given the provided number of valid ghosts
    * @param ghosts the number of ghosts, which cannot be larger than the number of valid ghosts on this field.
    *
    * This method should be used when a field assignment can only occur on a subset of the field due to invalidation of some of the ghost cells.
    */
   inline void reset_valid_ghosts( const GhostData& ghosts ){
     const IntVec diff = sfsharedPtr_->validGhosts_.get_minus() - ghosts.get_minus();
     fieldWindow_ = MemoryWindow( fieldWindow_.glob_dim(),
                                  fieldWindow_.offset() + diff,
                                  fieldWindow_.extent() - diff - sfsharedPtr_->validGhosts_.get_plus() + ghosts.get_plus() );
     sfsharedPtr_->reset_valid_ghosts( ghosts );
   }

   /**
    * @brief Obtain a child field that is reshaped.
    * @param extentModify the amount to modify the extent of the current field by
    * @param shift the number of grid points to shift the current field by
    * @return the reshaped child field
    *
    * The memory is the same as the parent field, but windowed differently.
    * Note that a reshaped field is considered read-only and you cannot obtain
    * interior iterators for these fields.
    */
   field_type
   inline reshape( const IntVec& extentModify,
                   const IntVec& shift ) const
   {
     MemoryWindow w( fieldWindow_.glob_dim(),
                     fieldWindow_.offset() + shift,
                     fieldWindow_.extent() + extentModify );
     return field_type( w, *this );
   }

  private:

   MemoryWindow fieldWindow_;
   BoundaryCellInfo bcInfo_;
   GhostData validGhosts_;
   typedef boost::shared_ptr<SpatialFieldLoc> shared_ptrT;
   shared_ptrT sfsharedPtr_;

  };  // SpatialField

//------------------------------------------------------------------

template<typename Location>
struct SingleValueCheck {
  static inline void check( MemoryWindow const & window,
                            GhostData const & ghosts ) {}
};

template<>
struct SingleValueCheck<SingleValue> {
  static inline void check( MemoryWindow const & window,
                            GhostData const & ghosts )
  {
#   ifndef NDEBUG
    if( (window.extent(0) > 1) &&
        (window.extent(1) > 1) &&
        (window.extent(2) > 1) ) {
      std::ostringstream msg;
      msg << "Single Value Field does not support window extents larger than 1\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }

    if( (ghosts.get_minus(0) > 0) &&
        (ghosts.get_minus(1) > 0) &&
        (ghosts.get_minus(2) > 0) &&
        (ghosts.get_plus(0) > 0) &&
        (ghosts.get_plus(1) > 0) &&
        (ghosts.get_plus(2) > 0) ) {
      std::ostringstream msg;
      msg << "Single Value Field does not support non-zero ghosts\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
#   endif
  }
};

//==================================================================
//
//                          Implementation
//
//==================================================================

template<typename Location, typename T>
SpatialField<Location,T>::SpatialFieldLoc::
SpatialFieldLoc( const MemoryWindow& window,
                 const BoundaryCellInfo& bc,
                 const GhostData& ghost,
                 T* const fieldValues,
                 const StorageMode mode,
                 const short int devIdx )
    : fieldWindow_(window),
      interiorFieldWindow_(window), // reset with correct info later
      bcInfo_( bc.limit_by_extent(window.extent()) ),
      ghosts_     ( ghost.limit_by_extent(window.extent()) ),
      validGhosts_( ghost.limit_by_extent(window.extent()) ),
      fieldValues_( (devIdx == CPU_INDEX) ?
                    ( ( mode == ExternalStorage) ? fieldValues : (NULL) )
                    : ( NULL ) ),
      fieldValuesExtDevice_( (IS_GPU_INDEX(devIdx)) ?
          // Note: this assumes fieldValues is on the proper GPU....
                             ( ( mode == ExternalStorage ) ? fieldValues : (NULL) ) // reset gpu memory later
                             : ( NULL ) ),
      builtField_( mode == InternalStorage ),
      activeDeviceIndex_( devIdx ),
      actualDeviceIndex_( devIdx ), // depricated with remove_multiple_fields()
      disableInterior_( false ),
      allocatedBytes_( 0 )
#     ifdef ENABLE_THREADS
      , partitionCount_( NTHREADS )
#     endif
#     ifdef ENABLE_CUDA
      , cudaStream_( 0 )
#     endif
{ //InteriorStorage => we build a new field
  //Exterior storage => we wrap T*
  // this error trapping is disabled currently because of the way that Wasatch is
  // hijacking the SpatialOps interface for flux limiters.  Once that gets folded
  // into a real nebo interface using operators we will be able to reinstate
  // these error trappings
//# ifndef NDEBUG
//  // ensure that we have a consistent BoundaryCellInfo object
//  for( int i=0; i<3; ++i ){
//    if( bcInfo_.has_bc(i) ){
//      assert( bcInfo_.num_extra(i) == Location::BCExtra::int_vec()[i] );
//    }
//  }
//# endif // NDEBUG
  SingleValueCheck<Location>::check(fieldWindow_, ghosts_);

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

  if( activeDeviceIndex_ == CPU_INDEX ){
    if( mode == InternalStorage ){
      fieldValues_ = Pool<T>::self().get( CPU_INDEX, window.glob_dim(0) * window.glob_dim(1) * window.glob_dim(2) );
    }
    multiFieldMap_[activeDeviceIndex_] = FieldInfo( fieldValues_, true, builtField_ );
  }
# ifdef ENABLE_CUDA
  else if( IS_GPU_INDEX(activeDeviceIndex_) ){
    if( mode == InternalStorage ){
      fieldValuesExtDevice_ = Pool<T>::self().get( activeDeviceIndex_, window.glob_dim(0) * window.glob_dim(1) * window.glob_dim(2) );
    }
    multiFieldMap_[activeDeviceIndex_] = FieldInfo( fieldValuesExtDevice_, true, builtField_ );
  }
# endif // ENABLE_CUDA
  else{
    std::ostringstream msg;
    msg << "Unsupported attempt to create field of type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
    << " )\n" << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  if (mode == InternalStorage ){
    reset_values(fieldValues);
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
SpatialField<Location,T>::SpatialFieldLoc::
SpatialFieldLoc( const SpatialFieldLoc& other )
: fieldWindow_(other.fieldWindow_),
  interiorFieldWindow_(other.interiorFieldWindow_),
  bcInfo_( other.bcInfo_ ),
  ghosts_( other.ghosts_ ),
  validGhosts_( other.validGhosts_ ),
  fieldValues_(other.fieldValues_),
  fieldValuesExtDevice_(other.fieldValuesExtDevice_),
  builtField_(false),
  activeDeviceIndex_(other.activeDeviceIndex_),
  actualDeviceIndex_(other.actualDeviceIndex_),
  disableInterior_( other.disableInterior_ ),
  allocatedBytes_( other.allocatedBytes_ )
# ifdef ENABLE_THREADS
, partitionCount_( other.partitionCount_ )
# endif
# ifdef ENABLE_CUDA
, cudaStream_( other.cudaStream_ )
# endif
{}

//------------------------------------------------------------------
// both the constructors can be taken out.
template<typename Location, typename T>
SpatialField<Location,T>::SpatialFieldLoc::
SpatialFieldLoc( const MemoryWindow& window, const SpatialFieldLoc& other )
: fieldWindow_(window),
  interiorFieldWindow_( other.interiorFieldWindow_ ), // This should not be used!
  bcInfo_( other.bcInfo_.limit_by_extent(window.extent()) ),
  ghosts_( other.ghosts_.limit_by_extent(window.extent()) ),
  validGhosts_( other.ghosts_.limit_by_extent(window.extent()) ),
  fieldValues_(other.fieldValues_),
  fieldValuesExtDevice_(other.fieldValuesExtDevice_),
  builtField_(false),
  activeDeviceIndex_(other.activeDeviceIndex_),
  actualDeviceIndex_(other.actualDeviceIndex_),
  disableInterior_( other.disableInterior_ ),
  allocatedBytes_( other.allocatedBytes_ )
# ifdef ENABLE_THREADS
, partitionCount_( other.partitionCount_ )
# endif
# ifdef ENABLE_CUDA
, cudaStream_( other.cudaStream_ )
# endif
{
  SingleValueCheck<Location>::check(fieldWindow_, ghosts_);

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

template<typename Location, typename T>
SpatialField<Location,T>::SpatialFieldLoc::
~SpatialFieldLoc()
{
  // abhi : When a field location with 'mode = ExternalStorage' has a additional
  //        field location from add_field_loc(), the mode will still be the same
  //        which prevents the memory from being released but the memory has to be released
  //        as the FieldInfo reports the _builtField is true.
  for(typename MultiFieldMap::iterator mapIter = multiFieldMap_.begin(); mapIter != multiFieldMap_.end(); mapIter++ ){
    if( !IS_VALID_INDEX(mapIter->first) ){
      std::ostringstream msg;
      msg << "Attempt to release field ( "
          << DeviceTypeTools::get_memory_type_description(mapIter->first)
          << " ) with invalid device Location. \n"
          <<     "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw( std::runtime_error( msg.str() ) );
    }
    if( builtField_ || mapIter->second._builtField ){
      //std::cout << "Putting memory back to the pool ~SpatialFieldLoc() : " << mapIter->second.field << std::endl;
      Pool<T>::self().put( mapIter->first, mapIter->second.field );
      mapIter->second.field = NULL;
      multiFieldMap_.erase(mapIter->first);
    }
  }
}
//------------------------------------------------------------------

template<typename Location, typename T>
SpatialField<Location,T>::~SpatialField()
{}

//------------------------------------------------------------------

template<typename FieldLocation, typename T>
void SpatialField<FieldLocation,T>::SpatialFieldLoc::
reset_values( const T* values )
{
  if( activeDeviceIndex_ == CPU_INDEX ){
    iterator ifld = begin(fieldWindow_);
    const iterator iflde = end(fieldWindow_);
    if ( values == NULL ) {
      for (; ifld != iflde; ++ifld)           *ifld = 0.0;
    } else {
      for (; ifld != iflde; ++ifld, ++values) *ifld = *values;
    }
  }
# ifdef ENABLE_CUDA
  else if( IS_GPU_INDEX(activeDeviceIndex_) ) {
    ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
    mapIter mapIter = multiFieldMap_.find(activeDeviceIndex_);
    if( values == NULL ) {
      CDI.memset( mapIter->second.field, 0, allocatedBytes_, activeDeviceIndex_ );
    } else {
      void* src = (void*)(values);
      CDI.memcpy_to( mapIter->second.field, src, allocatedBytes_, activeDeviceIndex_, cudaStream_ );
    }
  }
# endif // ENABLE_CUDA
  else{
    std::ostringstream msg;
    msg << "Reset values called for unsupported field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " )"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
T* SpatialField<Location,T>::SpatialFieldLoc::
field_values( const short int deviceLocation )
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to non-const SpatialField::field_values() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif

  // NOTE : check if the field being accessed is "ACTIVE"
  typename MultiFieldMap::iterator mapIter = multiFieldMap_.find( deviceLocation );

  if( mapIter == multiFieldMap_.end() ){
    if( !IS_VALID_INDEX( deviceLocation ) ){
      std::ostringstream msg;
      msg << "Request for field pointer on a device Location for which it has not been allocated\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    else{
      std::ostringstream msg;
      msg << "Request for field pointer on an unknown or unsupported device Location. \n"
          << "Please check the arguments passed into the function."
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
  }
  else{
    // Check if we are accessing an "active" field Location.
    if( deviceLocation != activeDeviceIndex_ ){
      std::ostringstream msg;
      msg << "The current field location is not set as active while accessing non-const version of field_values(). \n"
          << "Use set_field_loc_active() to set the device location.\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    else{
      // Make all the other device locations except active device Location as "INVALID"
      for(typename MultiFieldMap::iterator mapIter2 = multiFieldMap_.begin(); mapIter2 != multiFieldMap_.end(); mapIter2++ ){
        if( !(mapIter2->first == activeDeviceIndex_) ) mapIter2->second.isValid = false;
      }

      if( mapIter->second.field != NULL ) return mapIter->second.field;
      else{
        std::ostringstream msg;
        msg << "The field values from the non-const field_values() is NULL!\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
    }
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
const T* SpatialField<Location,T>::SpatialFieldLoc::
cfield_values( const short int deviceLocation ) const
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to const SpatialField::field_values() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif
  // Note : Make sure the field that is being accessed is not "Valid"
  typename MultiFieldMap::const_iterator mapIter = multiFieldMap_.find( deviceLocation );

  if( mapIter == multiFieldMap_.end()){
    if( !IS_VALID_INDEX( deviceLocation ) ){
      std::ostringstream msg;
      msg << "Request for field pointer on a device Location for which it has not been allocated\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    else{
      std::ostringstream msg;
      msg << "Request for field pointer on an unknown or unsupported device Location. \n"
          << "Please check the arguments passed into the function."
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
  }
  else{
    // Found the field pointer, check if it is "VALID"
    // cannot call non-const operation operator[] from a const member function.
    if( mapIter->second.isValid && mapIter->second.field != NULL ) return mapIter->second.field;
    else{
      std::ostringstream msg;
      msg << "Requested field pointer on a device Location is not validated! \n"
          << "Please check the synchronization options. \n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
  }
}

//------------------------------------------------------------------

#ifdef ENABLE_CUDA

template<typename Location, typename T>
void SpatialField<Location,T>::SpatialFieldLoc::
    add_field_loc( const short int deviceLocation )
{
#ifdef DEBUG_SF_ALL
  std::cout << "Call to SpatialField::add_field_loc() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif

  typename MultiFieldMap::iterator mapiter = multiFieldMap_.find( deviceLocation );
# ifdef ENABLE_THREADS
  //Make sure adding consumers is per-field atomic
  ExecMutex lock;
# endif

  if( mapiter != multiFieldMap_.end() ){
    if( mapiter->second.isValid ) return;
    else{
      sync_location( deviceLocation );
      return;
    }
  }
  else{
    multiFieldMap_[deviceLocation] = FieldInfo( Pool<T>::self().get(deviceLocation, (allocatedBytes_/sizeof(T))), false );
    sync_location(deviceLocation);
  }
}

//------------------------------------------------------------------

// Currently the sync_location() is called from add_field_loc() only
template<typename Location, typename T>
void SpatialField<Location, T>::SpatialFieldLoc::
     sync_location( const short int deviceLocation )
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to SpatialField::sync_location() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif

  // validate the deviceLocation from a valid field location
  if( !has_multiple_fields() ){
    if( !multiFieldMap_[deviceLocation].isValid ){
      std::ostringstream msg;
      msg << "sync_location() called on the field"
          << DeviceTypeTools::get_memory_type_description(deviceLocation)
          << " that didn't have any copies and it is invalid. \n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    return;
  }
  else{
    short int validIndex;
    T* validfieldValues_ = NULL;
    typename std::map<short int, FieldInfo>::iterator validIter = multiFieldMap_.find(activeDeviceIndex_);
    if( validIter != multiFieldMap_.end() ){
      // Convention : active Field Locations are always valid.
        validIndex  = validIter->first;
        validfieldValues_  = validIter->second.field;
    } else{
      std::ostringstream msg;
      msg << "Error : sync_location() didn't find a valid field entry in the map. \n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }

# ifndef NDEBUG
    if(validfieldValues_ == NULL){
      std::ostringstream msg;
      msg << "Error : valid fieldValues is NULL in sync_location() \n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }
# endif

    ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
    for( typename std::map<short int, FieldInfo>::iterator iter =  multiFieldMap_.begin(); iter != multiFieldMap_.end(); ++iter){
      if( !iter->second.isValid ){
        if( validIndex == CPU_INDEX && IS_GPU_INDEX(deviceLocation) ){
# ifdef DEBUG_SF_ALL
          std::cout << "data transfer from CPU to GPU (0/1)" << std::endl;
# endif
          CDI.memcpy_to( (void*)multiFieldMap_[deviceLocation].field, validfieldValues_, allocatedBytes_, deviceLocation, cudaStream_);
          iter->second.isValid = true;
          return;
        }
        else if( IS_GPU_INDEX(validIndex) && deviceLocation == CPU_INDEX ){
# ifdef DEBUG_SF_ALL
          std::cout << "data transfer from GPU (0/1) to CPU" << std::endl;
# endif
          CDI.memcpy_from( (void*)multiFieldMap_[deviceLocation].field, validfieldValues_, allocatedBytes_, validIndex, cudaStream_);
          iter->second.isValid = true;
          return;
        }
        else{
          std::ostringstream msg;
          msg << "Error : sync_location() called on the field"
              << DeviceTypeTools::get_memory_type_description(deviceLocation)
              << "for a GPU-GPU peer copy that is not yet supported. \n"
              << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
          throw(std::runtime_error(msg.str()));
        }
      }
    } // for loop
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
void SpatialField<Location,T>::SpatialFieldLoc::
  remove_multiple_fields()
{
  // Note this method is only accessed by fields that are locked.
# ifdef DEBUG_SF_ALL
  std::cout << "Caught call to SpatialField::remove_multiple_fields() for location : "
            << DeviceTypeTools::get_memory_type_description(actualDeviceIndex_) << std::endl;
# endif

  //Release any fields allocated for "multiple_fields" use
  if( !has_multiple_fields() ) return;

  for(typename MultiFieldMap::iterator mapIter = multiFieldMap_.begin(); mapIter != multiFieldMap_.end(); mapIter++ ){
    if( mapIter->first != actualDeviceIndex_ ){
      Pool<T>::self().put( mapIter->first, mapIter->second.field );
      mapIter->second.field = NULL;
      multiFieldMap_.erase(mapIter->first);
    }
  }
}

#endif // ENABLE_CUDA

//------------------------------------------------------------------

template<typename Location, typename T>
void SpatialField<Location,T>::SpatialFieldLoc::
set_field_loc_active( const short int deviceLocation )
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to SpatialField::set_field_loc_active() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif

# ifndef NDEBUG
  if( multiFieldMap_.find( deviceLocation ) == multiFieldMap_.end() ){
    std::ostringstream msg;
    msg << "Error : Requesting to set a field location as active that doesn't exist\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  // check if the field location that is active is indeed "VALID"
  if( !multiFieldMap_[deviceLocation].isValid ){
    std::ostringstream msg;
    msg << "Error : FieldLocation " << DeviceTypeTools::get_memory_type_description(deviceLocation)
        << " trying to set active, is not valid.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif
  activeDeviceIndex_ = deviceLocation;
}

//------------------------------------------------------------------

template<typename Location, typename T>
bool SpatialField<Location,T>::SpatialFieldLoc::
find_field_loc( const short int deviceLocation ) const
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to SpatialField::find_field_loc() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif

  // this logic is only used for the testing purpose only. Should be ripped out.
# ifndef NDEBUG
  if( multiFieldMap_.size() == 0 ){
    std::ostringstream msg;
    msg << "Error : Couldn't find an entry of the field in the map,"
        << DeviceTypeTools::get_memory_type_description(deviceLocation) << " field address, " << this << std::endl
        << ". This is a serious problem places to look for is the SpatialField constructor and add_fields_loc() method. \n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif

  typename MultiFieldMap::const_iterator mapIter = multiFieldMap_.find(deviceLocation);
  // Check if it is the only copy.
  if( multiFieldMap_.size() == 1 ){
    if( deviceLocation != activeDeviceIndex_ ){
      std::ostringstream msg;
      msg << "Error : Only one copy of a field exists, the active field location,\n"
          << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " is different \n"
          << "to that of the requested field, " << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }else{
      // we mostly use this method in the Nebo RHS (read). Hence, check to ensure
      // reading a NULL memory is not done.
      if( mapIter->second.isValid ){
        return mapIter->second.field != NULL;
      }else{
        std::ostringstream msg;
        msg << "Error : Requested field memory location,"
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << " is not valid. \n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
    }
  }
  else{
    if( mapIter == multiFieldMap_.end() ){
      std::ostringstream msg;
      msg << "Error : Requested field location," << DeviceTypeTools::get_memory_type_description(deviceLocation) << " doesn't exist.\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
    }else{
      if( mapIter->second.isValid ){
        return mapIter->second.field != NULL;
      }else{
        std::ostringstream msg;
        msg << "Error : Requested field memory location,"
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << " is not valid. \n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
      }
    }
  }

}

//------------------------------------------------------------------

template<typename Location, typename T>
bool SpatialField<Location,T>::SpatialFieldLoc::
is_valid( const short int deviceLocation ) const
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to SpatialField::is_valid() for field Location : "
            << DeviceTypeTools::get_memory_type_description(deviceLocation) << std::endl;
# endif

  cmapIter MapIter = multiFieldMap_.find( deviceLocation );
  if( MapIter == multiFieldMap_.end() ) return false;

  // check if it is valid
  if( !MapIter->second.isValid ) return false;
  else                           return true;
}

//------------------------------------------------------------------

template<typename Location, typename T>
bool SpatialField<Location,T>::SpatialFieldLoc::
  has_multiple_fields() const
{
# ifdef DEBUG_SF_ALL
  std::cout << "Call to SpatialField::has_multiple_fields() for field Location : "
            << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << std::endl;
# endif

  if( multiFieldMap_.size() > 1 ) return true;
  else                            return false;
}

//------------------------------------------------------------------

template<typename Location, typename T>
typename SpatialField<Location,T>::const_interior_iterator
SpatialField<Location,T>::SpatialFieldLoc::interior_end() const
{
  cmapIter const_mapIter = multiFieldMap_.find(CPU_INDEX);
# ifndef NDEBUG
  if( disableInterior_ ){
    std::ostringstream msg;
    msg << "Interior iterators cannot be obtained on resized fields" << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }
  if( const_mapIter == multiFieldMap_.end() || const_mapIter->second.field == NULL ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " doesn't exist in the map. Something wrong is with the memory. \n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif
  if( const_mapIter->second.isValid ) {
    const size_t extent = interiorFieldWindow_.extent(0) * interiorFieldWindow_.extent(1) * interiorFieldWindow_.extent(2);
    const_interior_iterator i(const_mapIter->second.field, interiorFieldWindow_);
    return i + extent;
  } else {
    std::ostringstream msg;
    msg << "Unsupported request for iterator to external field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " )"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
typename SpatialField<Location,T>::interior_iterator
SpatialField<Location,T>::SpatialFieldLoc::interior_end()
{
  mapIter mapIter = multiFieldMap_.find( CPU_INDEX );
# ifndef NDEBUG
  if( disableInterior_ ){
    std::ostringstream msg;
    msg << "Interior iterators cannot be obtained on resized fields" << std::endl
        << __FILE__ << " : " << __LINE__ << std::endl;
    throw( std::runtime_error(msg.str()) );
  }
  if( mapIter == multiFieldMap_.end() ){
    std::ostringstream msg;
    msg << "Field type "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
        << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif
  if( activeDeviceIndex_ == CPU_INDEX ){
    const size_t extent = interiorFieldWindow_.extent(0) * interiorFieldWindow_.extent(1) * interiorFieldWindow_.extent(2);
    interior_iterator i(mapIter->second.field, interiorFieldWindow_);
    return i + extent;
  }
  else{
    std::ostringstream msg;
    msg << "Unsupported request for iterator to external field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " )"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
T&
SpatialField<Location,T>::SpatialFieldLoc::
operator()( const size_t i, const size_t j, const size_t k )
{
  if( activeDeviceIndex_ == CPU_INDEX ){
    return (*this)(IntVec(i, j, k));
  }
  else{
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

template<typename Location, typename T>
T&
SpatialField<Location,T>::SpatialFieldLoc::operator()(const IntVec& ijk)
{
  mapIter mapIter = multiFieldMap_.find( CPU_INDEX );
# ifndef NDEBUG
  if( mapIter == multiFieldMap_.end() ){
    std::ostringstream msg;
    msg << "Field  type "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
        << "doesn't exist in the map. Something is wrong with the memory allocation.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif

  if( activeDeviceIndex_ == CPU_INDEX ){
#   ifndef NDEBUG
    assert( (size_t)ijk[0] <  fieldWindow_.extent(0) );
    assert( (size_t)ijk[1] <  fieldWindow_.extent(1) );
    assert( (size_t)ijk[2] <  fieldWindow_.extent(2) );
    assert( (size_t)ijk[0] >= fieldWindow_.offset(0) );
    assert( (size_t)ijk[1] >= fieldWindow_.offset(1) );
    assert( (size_t)ijk[2] >= fieldWindow_.offset(2) );
#   endif
    return mapIter->second.field[fieldWindow_.flat_index(ijk)];
  }
  else{
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
const T& SpatialField<Location,T>::SpatialFieldLoc::
operator()( const size_t i, const size_t j, const size_t k ) const
{
  cmapIter const_mapIter = multiFieldMap_.find(CPU_INDEX);
  if( (const_mapIter!=multiFieldMap_.end() || const_mapIter->second.field != NULL) && const_mapIter->second.isValid ){
#   ifndef NDEBUG
    assert( i <  fieldWindow_.extent(0) );
    assert( j <  fieldWindow_.extent(1) );
    assert( k <  fieldWindow_.extent(2) );
    assert( i >= fieldWindow_.offset(0) );
    assert( j >= fieldWindow_.offset(1) );
    assert( k >= fieldWindow_.offset(2) );
#   endif
    IntVec ijk(i,j,k);
    return const_mapIter->second.field[fieldWindow_.flat_index(ijk)];
  }
# ifndef NDEBUG
  else if( const_mapIter == multiFieldMap_.end() ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  else if( !const_mapIter->second.isValid ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " is invalid.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif
  else {
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

template<typename Location, typename T>
const T&
SpatialField<Location,T>::SpatialFieldLoc::
operator()( const IntVec& ijk ) const
{
  cmapIter const_mapIter = multiFieldMap_.find(CPU_INDEX);
  if( (const_mapIter == multiFieldMap_.end() || const_mapIter->second.field != NULL) && const_mapIter->second.isValid){
#   ifndef NDEBUG
    assert( (size_t)ijk[0] < fieldWindow_.extent(0) && (size_t)ijk[0] >= fieldWindow_.offset(0) );
    assert( (size_t)ijk[1] < fieldWindow_.extent(1) && (size_t)ijk[1] >= fieldWindow_.offset(1) );
    assert( (size_t)ijk[2] < fieldWindow_.extent(2) && (size_t)ijk[2] >= fieldWindow_.offset(2) );
#   endif
    return const_mapIter->second.field[fieldWindow_.flat_index(ijk)];
  }
# ifndef NDEBUG
  else if( const_mapIter == multiFieldMap_.end() ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  else if( !const_mapIter->second.isValid ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " is invalid.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif
  else {
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " does not support direct indexing.\n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
}

//------------------------------------------------------------------

template<typename Location, typename T>
T&
SpatialField<Location,T>::SpatialFieldLoc::operator[](const size_t i)
{
  mapIter mapIter = multiFieldMap_.find( CPU_INDEX );
# ifndef NDEBUG
  if( mapIter == multiFieldMap_.end() ){
    std::ostringstream msg;
    msg << "Field type "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_)
        << " doesn't exist in the map. Something is wrong with the memory allocation.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif

  if( activeDeviceIndex_ == CPU_INDEX ){
    return mapIter->second.field[i];
    //return fieldValues_[i];
  }
  else{
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
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
template<typename Location, typename T>
const T&
SpatialField<Location,T>::SpatialFieldLoc::operator[](const size_t i) const
{
  cmapIter mapIter = multiFieldMap_.find(CPU_INDEX);
# ifndef NDEBUG
  if( mapIter == multiFieldMap_.end() || mapIter->second.field == NULL ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " does not support direct indexing \n"
        << "Note: this function is DEPRECATED and is not recommended for future use.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  if( !mapIter->second.isValid ){
    std::ostringstream msg;
    msg << "Field type ( "
        << DeviceTypeTools::get_memory_type_description(activeDeviceIndex_) << " ) ,"
        << " is not valid. const-indexing is not allowed to a invalid field.\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }
# endif
  return mapIter->second.field[i];
}

//------------------------------------------------------------------

template<typename Location, typename T>
SpatialField<Location,T>&
SpatialField<Location,T>::operator=(const SpatialField& other)
{
  if( sfsharedPtr_->fieldWindow_ != other.sfsharedPtr_->fieldWindow_ ) {
    std::ostringstream msg;
    msg << "Error : Attempted assignment between fields of unequal size!\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  // check if the other is a valid field
  if( !other.is_valid( other.device_index() ) ){
    std::ostringstream msg;
    msg << "Error : Attempted to assign from field that is invalid !\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  short int currentIdx = other.device_index();

  if( IS_CPU_INDEX(sfsharedPtr_->activeDeviceIndex_) ) {
    if( IS_CPU_INDEX( currentIdx ) ){ // local_ram
      //Check for self assignment
      if( sfsharedPtr_->fieldValues_ == other.field_values() ){ return *this; }
      std::copy( other.sfsharedPtr_->fieldValues_, other.sfsharedPtr_->fieldValues_+other.sfsharedPtr_->fieldWindow_.glob_npts(), sfsharedPtr_->fieldValues_ );
    }

#   ifdef ENABLE_CUDA
    else if( IS_GPU_INDEX( currentIdx )){
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_from( sfsharedPtr_->fieldValues_, other.sfsharedPtr_->fieldValuesExtDevice_, sfsharedPtr_->allocatedBytes_, other.sfsharedPtr_->activeDeviceIndex_, sfsharedPtr_->cudaStream_ );
    }
#   endif

    else{
      std::ostringstream msg;
      msg << "Attempted unsupported copy operation, at n\t"
          << __FILE__ << " : " << __LINE__ << std::endl
          << "\t - "
          << DeviceTypeTools::get_memory_type_description(sfsharedPtr_->activeDeviceIndex_) << " = "
          << DeviceTypeTools::get_memory_type_description(
              other.device_index()) << std::endl;
      throw(std::runtime_error(msg.str()));
    }
    return *this;
  }
# ifdef ENABLE_CUDA
  else if( IS_GPU_INDEX(sfsharedPtr_->activeDeviceIndex_) ) {

    if( IS_CPU_INDEX( currentIdx ) ){ // local_ram
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_to( sfsharedPtr_->fieldValuesExtDevice_, other.sfsharedPtr_->fieldValues_, sfsharedPtr_->allocatedBytes_, sfsharedPtr_->activeDeviceIndex_, sfsharedPtr_->cudaStream_ );
    }

    else if( IS_GPU_INDEX( other.device_index() ) ){
      //Check for self assignment
      if( sfsharedPtr_->activeDeviceIndex_ == other.sfsharedPtr_->activeDeviceIndex_ && sfsharedPtr_->fieldValuesExtDevice_ == other.sfsharedPtr_->fieldValuesExtDevice_ ){
        return *this;
      }

      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      CDI.memcpy_peer( sfsharedPtr_->fieldValuesExtDevice_, sfsharedPtr_->activeDeviceIndex_, other.sfsharedPtr_->fieldValuesExtDevice_, other.sfsharedPtr_->activeDeviceIndex_, sfsharedPtr_->allocatedBytes_ );
    }

    else {
      std::ostringstream msg;
      msg << "Attempted unsupported copy operation, at " << std::endl
          << "\t" << __FILE__ << " : " << __LINE__ << std::endl
          << "\t - " << DeviceTypeTools::get_memory_type_description(sfsharedPtr_->activeDeviceIndex_) << " = "
          << DeviceTypeTools::get_memory_type_description(other.device_index()) << std::endl;
      throw( std::runtime_error ( msg.str() ));
    }

    } // else if
# endif // ENABLE_CUDA
  else{
    std::ostringstream msg;
    msg << "Attempted unsupported copy operation, at \n\t"
        << __FILE__ << " : " << __LINE__ << std::endl
        << "\t - " << DeviceTypeTools::get_memory_type_description(sfsharedPtr_->activeDeviceIndex_)
        << " = "
        << DeviceTypeTools::get_memory_type_description(
        other.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  } // else
  return *this;
}

//------------------------------------------------------------------

  /**
   *  \fn BoundaryCellInfo create_new_boundary_cell_info<FieldType>( const PrototypeType )
   *
   *  \brief create a boundary cell info for a field of type FieldType from a field of type PrototypeType
   *
   *  \param prototype the prototype field
   *
   *  \return the new boundary cell info
   */
  template<typename FieldType, typename PrototypeType>
  inline BoundaryCellInfo create_new_boundary_cell_info(const PrototypeType & prototype) {
    return BoundaryCellInfo::build<FieldType>(prototype.boundary_info().has_bc());
  }

//------------------------------------------------------------------

  /**
   *  \fn MemoryWindow create_new_memory_window<FieldType>( const PrototypeType )
   *
   *  \brief create a memory window for a field of type FieldType from a field of type PrototypeType
   *
   *  \param prototype the prototype field
   *
   *  \return the new memory window (with correct boundary conditions)
   */
  template<typename FieldType, typename PrototypeType>
  inline MemoryWindow create_new_memory_window(const PrototypeType & prototype) {
    const BoundaryCellInfo prototypeBC = prototype.boundary_info();
    const BoundaryCellInfo newBC = create_new_boundary_cell_info<FieldType, PrototypeType>(prototype);

    const MemoryWindow prototypeWindow = prototype.window_with_ghost();

    return MemoryWindow(prototypeWindow.glob_dim() - prototypeBC.has_extra() + newBC.has_extra(),
                        prototypeWindow.offset(),
                        prototypeWindow.extent() - prototypeBC.has_extra() + newBC.has_extra());
  }

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_SpatialField_h
