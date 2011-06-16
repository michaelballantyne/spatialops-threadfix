#include <spatialops/ThreadPool.h>

namespace SpatialOps {

   const int ThreadPoolResourceManager::maxWorkerThreads_;
   int ThreadPoolResourceManager::availableThreadingResources_ = NTHREADS;
   std::map<void*, int> ThreadPoolResourceManager::resourceMap_;
   bool ThreadPool::init = false;
   bool ThreadPoolFIFO::init = false;

}
