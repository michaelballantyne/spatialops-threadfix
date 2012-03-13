#ifndef ParticleFieldTypes_h
#define ParticleFieldTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps{
namespace Particle{


  /**
   *  @file ParticleFieldTypes.h
   *
   *  Particle fields are dimensioned by the number of particles.  They
   *  are distinctly different types than fields on the underlying mesh,
   *  and we must define operators to move between particle fields and
   *  mesh fields.
   */

  struct ParticleGhostTraits{ enum{ NGHOST=0 }; };

  struct ParticleFieldTraits{
    typedef NODIR FaceDir;
    typedef structured::IndexTriplet< 0, 0, 0> Offset;
    typedef structured::IndexTriplet<0,0,0>  BCExtra;
  };

  typedef structured::SpatialField< ParticleFieldTraits, ParticleGhostTraits  > ParticleField;


} // namespace Particle
} // namespace SpatialOps

#endif // ParticleFieldTypes_h
