#include <spatialops/ThreadPool.h>

namespace SpatialOps {
   std::map<void*, int> ThreadPoolResourceManager::resourceMap_;
   bool ThreadPool::init = false;
   bool ThreadPoolFIFO::init = false;
}
