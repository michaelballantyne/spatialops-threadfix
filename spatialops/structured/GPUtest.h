#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/Nebo.h>

typedef SpatialOps::structured::SVolField Field;

#ifdef ENABLE_CUDA
void addsin(Field & result,
            Field const & src1,
            Field const & src2);
#else
void addsin(Field & result,
            Field const & src1,
            Field const & src2) {
    using namespace SpatialOps;
    result <<= src1 + sin(src2);
};
#endif
void test(int nx, int ny, int nz);
