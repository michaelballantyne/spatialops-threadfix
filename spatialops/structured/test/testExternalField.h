#ifndef TESTEXTERNALFIELD_H_
#define TESTEXTERNALFIELD_H_

#define DEBUG_SPATIAL_FIELD
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/structured/MemoryWindow.h>

namespace SpatialOps {
namespace Point {

struct PointFieldGhostTraits {
    typedef SpatialOps::structured::IndexTriplet<0,0,0> NGhostMinus;
    typedef SpatialOps::structured::IndexTriplet<0,0,0> NGhostPlus;
};

struct PointFieldTraits { typedef NODIR FaceDir; typedef NODIR StagLoc; };

/**
 *  \brief The PointField type is intended for use in extracting and
 *         working with individual points from within another field type.
 *
 *  This field type is not compatible with operations such as
 *  interpolants, gradients, etc.  Operators are provided to extract
 *  points from a parent field and return them back to a parent
 *  field.
 */
typedef structured::SpatialField<Point::PointFieldTraits,
    Point::PointFieldGhostTraits> PointField;

} // namespace Point
} // namespace SpatialOps

#endif /* TESTEXTERNALFIELD_H_ */
