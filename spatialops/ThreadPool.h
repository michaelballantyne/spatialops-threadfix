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

#ifndef Field_Expr_ThreadPool_h
#define Field_Expr_ThreadPool_h

#include <stdio.h>
#include <map>
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/threadpool/threadpool.hpp>

#include <boost/thread/mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

#ifdef  FIELD_EXPRESSION_THREADS
#define USE_THREADS_FLAG
#endif

#ifdef  TILL_THREADS
#define USE_THREADS_FLAG
#endif

#ifdef  STENCIL_THREADS
#define USE_THREADS_FLAG
#endif

#ifndef USE_THREADS_FLAG
#  error Must define FIELD_EXPRESSION_THREADS, TILL_THREADS, or STENCIL_THREADS flag.
#endif

#ifndef NTHREADS
#  error Must define NTHREADS flag.
#endif

namespace SpatialOps{

  class ThreadPoolResourceManager {
    typedef std::map<void*,int>::iterator ResourceIter;
    ThreadPoolResourceManager(){};
    ~ThreadPoolResourceManager(){};

  public:
    static ThreadPoolResourceManager& self(){
      static ThreadPoolResourceManager tprm;
      return tprm;
    }

    template<class VoidType>
    static const bool insert(VoidType& rID, int threads){
      ExecutionMutex lock;
      if( threads < 1 ) { threads = 1; }
      //Make sure we don't have the threadpool
      ResourceIter rit = resourceMap_.find(&rID);
      if ( rit == resourceMap_.end() ){
        resourceMap_.insert(std::make_pair(&rID, threads));
      } else {
        printf("Warning: attempting to insert a ThreadPool that already exists!\n");
        return false;
      }

      return true;
    }

    template<class VoidType>
    static const bool remove(VoidType& rID, int threads){
      ExecutionMutex lock;
      ResourceIter rit = resourceMap_.find(&rID);
      if ( rit != resourceMap_.end() ){
        resourceMap_.erase(rit);
      } else {
        printf("Warning: attempting to remove ThreadPool that does not exist!\n");
        return false;
      }
      return true;
    }

    template<class VoidType>
    static const int resize(VoidType& rID, int threads){
      ExecutionMutex lock;
      //Make sure we have the threadpool
      ResourceIter rit = resourceMap_.find(&rID);
      if( rit == resourceMap_.end() ) {
        fprintf(stderr, "Error: ThreadPool does not exist!\n");
        return -1;
      }

      //Fast exit
      if( rit->second == threads ) { return threads; }

      //Connect the right resource interface
      VoidType* resource = (VoidType*)rit->first;

      if( threads < 1 ) { threads = 1; }
      rit->second = threads;
      resource->size_controller().resize(threads);
      return threads;
    }

    template<class VoidType>
    static const int resize_active(VoidType& rID, int threads){
      ExecutionMutex lock;
      //Make sure we have the threadpool
      ResourceIter rit = resourceMap_.find(&rID);
      if( rit == resourceMap_.end() ) {
        fprintf(stderr, "Error: ThreadPool does not exist!\n");
        return -1;
      }
      //Connect the right resource interface
      VoidType* resource = (VoidType*)rit->first;
      resource->size_controller().resize_active(std::max(1,threads));
      return threads;
    }

    template<class VoidType>
    static const int get_worker_count(VoidType& rID) {
      ExecutionMutex lock;
      ResourceIter rit = resourceMap_.find(&rID);
      if( rit == resourceMap_.end() ) {
        fprintf(stderr, "Error: Threadpool does not exist!\n");
        return -1;
      }
      VoidType* resource = (VoidType*)rit->first;
      return resource->size();
    }

    template<class VoidType>
    static const int get_max_active_worker_count(VoidType& rID) {
      ExecutionMutex lock;
      //Make sure we have the threadpool
      ResourceIter rit = resourceMap_.find(&rID);
      if( rit == resourceMap_.end() ) {
        fprintf(stderr, "Error: ThreadPool does not exist!\n");
        return -1;
      }

      //Connect the right resource interface
      VoidType* resource = (VoidType*)rit->first;
      return resource->max_active();
    }

  private:
    static std::map<void*, int> resourceMap_;

    class ExecutionMutex{
      const boost::mutex::scoped_lock lock;

      inline boost::mutex& get_mutex() const { static boost::mutex m; return m; }

    public:
      ExecutionMutex() : lock( get_mutex() ){}
      ~ExecutionMutex(){}
    };
  };

  /**
   * \brief Wrapper for a priority thread pool
   */
  class ThreadPool : public boost::threadpool::prio_pool{
    ThreadPool( const int nthreads ) : boost::threadpool::prio_pool( nthreads ){}
    ~ThreadPool(){}

  public:
    static ThreadPool& self() {
      static int threads = NTHREADS < 1 ? 1 : NTHREADS;
      static ThreadPool tp(threads);
      ThreadPoolResourceManager& tprm = ThreadPoolResourceManager::self();
      if( init == false ){
        tprm.insert<boost::threadpool::prio_pool>(tp, threads);
        init = true;
      }
      return tp;
    }

  private:
    static bool init;
  };

  /**
   * \brief Wrapper for a FIFO thread pool
   */
  class ThreadPoolFIFO : public boost::threadpool::fifo_pool{
    ThreadPoolFIFO( const int nthreads ) : boost::threadpool::fifo_pool( nthreads ){}
    ~ThreadPoolFIFO(){}

  public:

    static ThreadPoolFIFO& self(){
      static int threads = NTHREADS < 1 ? 1 : NTHREADS;
      static ThreadPoolFIFO tp(threads);
      ThreadPoolResourceManager& tprm = ThreadPoolResourceManager::self();
      if( init == false ){
        tprm.insert<boost::threadpool::fifo_pool>(tp, threads);
        init = true;
      }
      return tp;
    }

    static bool is_nebo_thread_parallel() { return nebo_parallelism; };

    static bool set_nebo_thread_parallel(bool newBool) {
      nebo_parallelism = newBool;
      return nebo_parallelism;
    };

  private:
    static bool init;
    static bool nebo_parallelism;
  };
} // namespace SpatialOps

#endif // Field_Expr_ThreadPool_h
