#ifndef PointFieldTypes_h
#define PointFieldTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/structured/SpatialField.h>

namespace SpatialOps{
namespace Point{

  struct LinAlg{
    typedef int VecType;
    VecType& setup_vector( const int, double* ){ static VecType vt=0; return vt; }
  };
 
  struct PointFieldGhostTraits{ enum{ NGHOST=0 }; };
  
  struct PointFieldTraits{ typedef NODIR FaceDir; typedef NODIR StagLoc; };

  /**
   *  \brief The PointField type is intended for use in extracting and
   *         working with individual points from within another field type.
   *
   *  This field type is not compatible with operations such as
   *  interpolants, gradients, etc.  Operators are provided to extract
   *  points from a parent field and return them back to a parent
   *  field.
   */  
  typedef structured::SpatialField< Point::LinAlg,
                                    Point::PointFieldTraits,
                                    Point::PointFieldGhostTraits >  PointField;
  
}  // namespace Point
} // namespace SpatialOps

#endif // PointFieldTypes_h
