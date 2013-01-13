#include <iostream>
#include <vector>

//--- SpatialOps includes ---//
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/FieldExpressions.h>

#include <test/FieldHelper.h>
#include <spatialops/structured/GPUtest.h>

using namespace SpatialOps;

void addsin(Field & result,
            Field const & src1,
            Field const & src2) {
    result |= src1 + sin(src2);
};
