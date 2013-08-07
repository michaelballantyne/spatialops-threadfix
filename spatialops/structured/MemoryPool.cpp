/**
 *  \file   MemoryPool.cpp
 *  \date   Jul 17, 2013
 *  \author "James C. Sutherland"
 *
 *
 * The MIT License
 *
 * Copyright (c) 2013 The University of Utah
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
 *
 */

#include <spatialops/structured/MemoryPool.h>
#include <spatialops/structured/SpatialField.h>

#include <spatialops/structured/CudaMemoryAllocator.h>

#include <boost/type_traits.hpp>

namespace SpatialOps{
namespace structured{


template< typename T >
Pool<T>::Pool() : deviceIndex_(0)
{
  destroyed_ = false;
  pad_ = 0;
  cpuhighWater_ = 0;
# ifdef ENABLE_CUDA
  gpuhighWater_ = 0;
# endif
}

template< typename T >
Pool<T>::~Pool()
{
  destroyed_ = true;

  for( typename FQSizeMap::iterator i=cpufqm_.begin(); i!=cpufqm_.end(); ++i ){
    FieldQueue& fq = i->second;
    while( !fq.empty() ){
      delete [] fq.top();
      fq.pop();
    }
  }

# ifdef ENABLE_CUDA
  for( typename FQSizeMap::iterator i=gpufqm_.begin(); i!=gpufqm_.end(); ++i ){
    FieldQueue& fq = i->second;
    while( !fq.empty() ){
            ema::cuda::CUDADeviceInterface::self().release( fq.top(), deviceIndex_ );
            fq.pop();
    }
  }
# endif
}

template< typename T >
Pool<T>&
Pool<T>::self()
{
  //ensure CUDA driver is loaded before pool is initialized
  static cudaError_t junk = cudaGetDeviceCount(0);
  // see Modern C++ (Alexandrescu) chapter 6 for an excellent discussion on singleton implementation
  static Pool<T> p;
  return p;
}

template< typename T >
T*
Pool<T>::get( const MemoryType mtype, const size_t _n )
{
  assert( !destroyed_ );

  if( pad_==0 ) pad_ = _n/10;
  size_t n = _n+pad_;

  switch(mtype) {
  case LOCAL_RAM: {
    T* field = NULL;
    typename FQSizeMap::iterator ifq = cpufqm_.lower_bound( n );
    if( ifq == cpufqm_.end() ){
      ifq = cpufqm_.insert( ifq, make_pair(n,FieldQueue()) );
    }
    else{
      n = ifq->first;
    }
    FieldQueue& fq = ifq->second;
    if( fq.empty() ){
      ++cpuhighWater_;
      try{
        field = new T[n];
      }
      catch(std::runtime_error& e){
        std::cout << "Error occurred while allocating memory on LOCAL_RAM" << std::endl
                  << e.what() << std::endl
                  << __FILE__ << " : " << __LINE__ << std::endl;
      }
      fsm_[field] = n;
    }
    else{
      field = fq.top(); fq.pop();
    }
    return field;
  }
# ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    T* field = NULL;
    typename FQSizeMap::iterator ifq = gpufqm_.lower_bound( n );
    if( ifq == gpufqm_.end() ){
      ifq = gpufqm_.insert( ifq, make_pair(n,FieldQueue()) );
    }
    else{
      n = ifq->first;
    }
    FieldQueue& fq = ifq->second;
    if( fq.empty() ) {
      ++gpuhighWater_;
      ema::cuda::CUDADeviceInterface& CDI = ema::cuda::CUDADeviceInterface::self();
      field = (T*)CDI.get_raw_pointer( n * sizeof(T), deviceIndex_ );
      fsm_[field] = n;
    }
    else{
      field = fq.top(); fq.pop();
    }
    return field;
  }
# endif
  default: {
    std::ostringstream msg;
      msg << "Attempt to get unsupported memory pool ( "
          << DeviceTypeTools::get_memory_type_description(mtype)
          << " ) \n";
      msg << "\t " << __FILE__ << " : " << __LINE__;
      throw(std::runtime_error(msg.str()));
  }
  } //switch
}

template< typename T >
void
Pool<T>::put( const MemoryType mtype, T* t )
{
  // in some cases (notably in the LBMS code), singleton destruction order
  // causes the pool to be prematurely deleted.  Then subsequent calls here
  // would result in undefined behavior.  If the singleton has been destroyed,
  // then we will just ignore calls to return resources to the pool.  This will
  // leak memory on shut-down in those cases.
  if( destroyed_ ) return;

  switch(mtype) {
  case LOCAL_RAM: {
    const size_t n = fsm_[t];
    const typename FQSizeMap::iterator ifq = cpufqm_.lower_bound( n );
    assert( ifq != cpufqm_.end() );
    ifq->second.push(t);
    break;
  }
# ifdef ENABLE_CUDA
  case EXTERNAL_CUDA_GPU: {
    const size_t n = fsm_[t];
    const typename FQSizeMap::iterator ifq = gpufqm_.lower_bound( n );
    assert( ifq != gpufqm_.end() );
    ifq->second.push(t);
    break;
  }
# endif
  default: {
    std::ostringstream msg;
        msg << "Error occurred while restoring memory back to memory pool ( "
            << DeviceTypeTools::get_memory_type_description(mtype)
            << " ) \n";
        msg << "\t " << __FILE__ << " : " << __LINE__;
        throw(std::runtime_error(msg.str()));
  }
  }//switch
}

template<typename T>
size_t
Pool<T>::active() const{
  size_t n=0;
  for( typename FQSizeMap::const_iterator ifq=cpufqm_.begin(); ifq!=cpufqm_.end(); ++ifq ){
    n += ifq->second.size();
  }
  return cpuhighWater_-n;
}


// explicit instantiation
template class Pool<double>;
template class Pool<float>;

}
}
