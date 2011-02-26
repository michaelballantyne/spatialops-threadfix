#ifndef FVStaggered_h
#define FVStaggered_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/OperatorDatabase.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVStaggeredOperatorTypes.h>

#include <spatialops/structured/matrix/FVStaggeredBCOp.h>

#include <spatialops/structured/matrix/FVStaggeredInterpolant.h>
#include <spatialops/structured/matrix/FVStaggeredGradient.h>
#include <spatialops/structured/matrix/FVStaggeredDivergence.h>
#include <spatialops/structured/matrix/FVStaggeredScratch.h>
#include <spatialops/structured/matrix/FVTopHatFilter.h>
#include <spatialops/structured/matrix/FVRestrictOp.h>

#include <spatialops/structured/SpatialFieldStore.h>

#include <spatialops/FieldExpressionsExtended.h>
#include <spatialops/FieldReductions.h>

namespace SpatialOps{
namespace structured{

/**
 *  \file FVStaggered.h
 *  \defgroup fields		Field definitions and tools
 *  \defgroup operators		Operator definitions and tools
 *  \defgroup structured	Tools and types for structured meshes
 *
 *
 *  \todo Need to build & test scratch operator assemblers.
 *  \todo Need to fix linear system interface for spatial operators.
 *
 *
 *
 *  \page structured Structured Mesh Tools
 *
 *  SpatialOps has been initially targeted at providing facilities for
 *  structured meshes.  Although data structures are simpler for
 *  structured meshes, in the case of a staggered mesh, there are many
 *  more field types and operators for structured meshes than for
 *  unstructured ones.
 *
 *  Most codes intending to use SpatialOps for structured meshes should include the following:
 *  \code
 *  #include <spatialops/structured/FVStaggered.h>
 *  \endcode
 *  This defines fields and operators for use on a staggered mesh.
 *
 *
 *  \section structured_key_concepts Key Concepts
 *
 *  \section structured_key_classes Key Classes  
 *    - MemoryWindow
 *    - SpatialField
 *
 *
 *  \section Examples
 *
 *  \code
 *  #include <spatialops/FVStaggered.h>
 *  using namespace SpatialOps;
 *  using namespace structured;
 *
 *  MemoryWindow vw( IntVec(10,1,1) );  // create a 1-D window for a field with 10 points
 *  SVolField v1( vw, NULL );           // create a field with 10 points in it.
 *  SVolField v2( vw, NULL );           // create a field with 10 points in it.
 *
 *  // populate the fields
 *  v1 <<= 5.0*exp(2.0);
 *  v2 << v1 / ( 2.0 + sin(v1) );
 *  \endcode
 */

} // namespace structured
} // namespace SpatialOps

#endif
