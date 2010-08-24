#ifndef SpatialOpsTools_h
#define SpatialOpsTools_h

#include <spatialops/SpatialOpsConfigure.h>

namespace SpatialOps{

  template<typename T1, typename T2>
  struct IsSameType{ enum{ result=0 }; };

  template< typename T1 >
  struct IsSameType<T1,T1>{ enum{ result=1 }; };


  template<typename FieldT>
  inline size_t nghost(){ return 2*FieldT::Ghost::NGHOST; }

  template<>
  inline size_t nghost<double>(){ return 0; }

}

#endif

