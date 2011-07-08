#ifndef PointFieldTypes_h
#define PointFieldTypes_h

#include <spatialops/structured/SpatialField.h>

namespace SpatialOps{
namespace Point{

  struct LinAlg{
    typedef int VecType;
    VecType& setup_vector( const int, double* ){ static VecType vt=0; return vt; }
  };
 
  struct PointFieldGhostTraits{ enum{ NGHOST=0 }; };
  
  struct PointFieldTraits{};
  
  typedef structured::SpatialField< Point::LinAlg,
                                    Point::PointFieldTraits,
                                    Point::PointFieldGhostTraits >  PointField;
  
}  // namespace Point
} // namespace SpatialOps

#endif // PointFieldTypes_h
