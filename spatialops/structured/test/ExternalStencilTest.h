/*
 * ExternalStencilTest.h
 *
 *  Created on: Jan 5, 2012
 *      Author: Devin Robison
 */

#ifndef EXTERNALSTENCILTEST_H_
#define EXTERNALSTENCILTEST_H_
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps {
	namespace Point {

	struct PointFieldGhostTraits {
		typedef structured::IndexTriplet<0,0,0> NGhostMinus;
		typedef structured::IndexTriplet<0,0,0> NGhostPlus;
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
		Point::PointFieldGhostTraits, float> PointFloatField;

	} // namespace Point
} // namespace SpatialOps



#endif /* EXTERNALSTENCILTEST_H_ */
