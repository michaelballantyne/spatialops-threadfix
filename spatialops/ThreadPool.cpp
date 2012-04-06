#include <spatialops/ThreadPool.h>

namespace SpatialOps {
   std::map<void*, int> ThreadPoolResourceManager::resourceMap_;
   bool ThreadPool::init = false;
   bool ThreadPoolFIFO::init = false;
   bool ThreadPoolFIFO::nebo_parallelism = NTHREADS > 0;
}
