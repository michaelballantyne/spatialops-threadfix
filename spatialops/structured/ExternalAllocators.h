/*
 * ExternalMemoryAllocators.h
 *
 *  Created on: Nov 21, 2011
 *      Author: Devin Robison
 */

#ifndef EXTERNALALLOCATORS_H_
#define EXTERNALALLOCATORS_H_
#include <spatialops/SpatialOpsConfigure.h>

namespace ema {
#ifdef ENABLE_CUDA
namespace cuda { //ema::cuda

class CUDADeviceInterface;  // Intermediate structure to talk to CUDA capable devices
class CUDASharedPointer;

/** \brief Structure for storing memory statistics **/
struct CUDAMemStats{
    size_t f;
    size_t t;
    size_t inuse;
};

/** \brief Interface to the CUDADeviceManager
 *  This is used to obtain device information and memory allocation
 *  without forcing a recompile of libraries compiled with nvcc and
 *  prevents recompiling existing code that needs to pass GPU Device
 *  pointers.
 */

class CUDADeviceInterface {
    friend class CUDASharedPointer;

    CUDADeviceInterface();
  public:
    ~CUDADeviceInterface();

    /** \brief Returns a pointer to the global CU device interface object **/
    static CUDADeviceInterface& self();

    /** \brief attempts to allocate N bytes on device K */
    CUDASharedPointer get(unsigned long int N, int K = 0);

    /** \brief attempts to free an allocated shared pointer object */
    void release(CUDASharedPointer& x);

    /** \brief copy a data block into a cuda pointer **/
    void copy_to(CUDASharedPointer& dest, void* src, size_t sz);

    /** \brief copy a data block from a cuda pointer **/
    void copy_from(void* dest, CUDASharedPointer& src, size_t sz);

    /** \brief Returns the number of available CUDA capable compute devices */
    int get_device_count() const;

    /** \brief Returns the memory structure associted with device K */
    const CUDAMemStats& get_memory_statistics(int K = 0) const;

    /** \brief Updates the 'device_stats' structs with the most current memory usage statistics
     * Please note that it is possible memory can be allocated from other sources, this
     * is simply to provide a metric for checking relative memory availability.
     * */
    void update_memory_statistics();

    /** \brief output a list of all available CUDA hardware and capabilities */
    void print_device_info() const;

  private:

};

/** Not sure if were going to use this anymore as it duplicates
 *SpatialFieldPtr functionality, but it may be handy to have for now.
/* \brief Wrapper structure for a ref-counted GPU memory pointer */
class CUDASharedPointer {
    friend class CUDADeviceInterface;

    CUDASharedPointer(void*, int);

  public:
    ~CUDASharedPointer();
    CUDASharedPointer();
    CUDASharedPointer(const CUDASharedPointer& x);

    /** \brief test for equality between two shared pointers **/
    bool operator==(const CUDASharedPointer& x) const;

    /** \brief assign this pointer to another CUDASharedPointer **/
    CUDASharedPointer& operator=( const CUDASharedPointer& x );

    /** \brief assign this pointer to another CUDASharedPointer
     *
     * IMPORTANT: This will default pointer device to the current thread's device context
     * If this is not what you want, then you will need to manually assign the pointer deviceID.
     *  **/
    CUDASharedPointer& operator=( void* x );

    /** \brief detaches the pointer from what it references, setting it to NULL, returns
     * the NULL pointer
     * **/
    CUDASharedPointer& detatch();

    /** \brief dereference operators, return the location pointed to by location_ **/
    void* operator->();
    const void* operator->() const;

    /** \brief does this point to NULL **/
    bool isnull() const;

    /** \brief return deviceID **/
    int get_deviceID() const;

    /** \brief return remaining references **/
    int get_refcount() const;

  private:
    int* deviceID_;  ///< Device ID this pointer is valid on
    int* refCount_;  ///< Number of shared pointers referencing this position

    void* ptr_; ///< GPU memory pointe (valid ONLY on deviceID_)
};

} // End Namespace 'ema::cuda'
#endif // ENABLE_CUDA
} // End Namespace 'ema'

#endif /* EXTERNALMEMORYALLOCATORS_H_ */
