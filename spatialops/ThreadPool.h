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
         static bool insert(VoidType& rID, int threads){
            ResourceIter rit;
            ExecutionMutex lock;

            //printf("Inserting ThreadPool<%s><0x%x>\n", typeid(rID).name(), &rID);
            //Make sure we don't have the threadpool
            rit = resourceMap_.find(&rID);
            if ( rit == resourceMap_.end() ){
               resourceMap_.insert(std::make_pair(&rID, std::min(maxWorkerThreads_, threads))); 
            } else {
               printf("Warning: attempting to insert a ThreadPool that already exists!\n");
               return false;
            }

            return true;
         }
         
         // Can the boost threadpool's smart pointer ever reach zero when ThreadPool is static?
         template<class VoidType>
         static bool remove(VoidType& rID, int threads){
            ResourceIter rit;
            ExecutionMutex lock;

            rit = resourceMap_.find(&rID);
            if ( rit != resourceMap_.end() ){
               resourceMap_.erase(rit);
            } else {
               printf("Warning: attempting to remove ThreadPool that does not exist!\n");
               return false;
            }

            return true;
         }

         template<class VoidType>
         static int resize(VoidType& rID, int threads){
            VoidType* resource;
            ResourceIter rit;
            ExecutionMutex lock;

            //Make sure we have the threadpool
            rit = resourceMap_.find(&rID);
            if( rit == resourceMap_.end() ) { 
               fprintf(stderr, "Error: ThreadPool does not exist!\n"); 
               return -1;
            }

            //Fast exit 
            if( rit->second == threads ) { return threads; }
            
            //Connect the right resource interface
            resource = (VoidType*)rit->first;
            
            //Clamp
            threads = std::min(maxWorkerThreads_, threads);
            rit->second = threads;
            resource->size_controller().resize(threads);

            return threads;
         }

         static const int getMaxThreadingResources(){ return maxWorkerThreads_; }
         static const int getThreadingResources(){ return availableThreadingResources_; }
      private:
         static const int maxWorkerThreads_ = NTHREADS;
         static int availableThreadingResources_;
         static std::map<void*, int> resourceMap_;

         class ExecutionMutex{
            const boost::mutex::scoped_lock lock;
            
            inline boost::mutex& get_mutex() const { static boost::mutex m; return m; }

            public:
               ExecutionMutex() : lock( get_mutex() ){}
               ~ExecutionMutex(){}
         };
   };

  class ThreadPool : public boost::threadpool::prio_pool{
      ThreadPool( const int nthreads ) : boost::threadpool::prio_pool( nthreads ){}
      ~ThreadPool(){}

      public:
         static ThreadPool& self(){
            static ThreadPool tp(1);
            ThreadPoolResourceManager& tprm = ThreadPoolResourceManager::self();
            if( init == false ){ 
               tprm.insert<boost::threadpool::prio_pool>(tp, NTHREADS);
               init = true; 
            }
            return tp;
         }

      private:
         static bool init;
  };

  class ThreadPoolFIFO : public boost::threadpool::fifo_pool{
      ThreadPoolFIFO( const int nthreads ) : boost::threadpool::fifo_pool( nthreads ){}
      ~ThreadPoolFIFO(){}
      
     public:

         static ThreadPoolFIFO& self(){
            static ThreadPoolFIFO tp(NTHREADS);
            ThreadPoolResourceManager& tprm = ThreadPoolResourceManager::self();
            if( init == false ){ 
               tprm.insert<boost::threadpool::fifo_pool>(tp, NTHREADS); 
               init = true; 
            }
            return tp;
         }

      private:
         static bool init;
  };
} // namespace SpatialOps

#endif // Field_Expr_ThreadPool_h
