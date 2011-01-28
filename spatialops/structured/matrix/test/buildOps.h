#ifndef buildOps_h
#define buildOps_h

#include <vector>
namespace SpatialOps{ class OperatorDatabase; } // forward

void build_ops( const SpatialOps::structured::IntVec& dim,
                const std::vector<double>& spacing,
                const std::vector<bool>& bcFlag,
                SpatialOps::OperatorDatabase& opDB );

#endif
