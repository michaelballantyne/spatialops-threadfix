#ifndef ParticleFieldTypes_h
#define ParticleFieldTypes_h

#include <spatialops/SpatialField.h>

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

  struct ParticleFieldTraits{};

  class ParticleLinAlg
  {
  public:
    typedef double* VecType;  // required for compatibility
    ParticleLinAlg();
    VecType& setup_vector(const size_t,double*const);
    void destroy_vector();
    void print_vec( std::ostream& s ) const;
    ~ParticleLinAlg();
  private:
    VecType vec;
    size_t length;
  };

  typedef SpatialOps::SpatialField<
    ParticleLinAlg,
    ParticleFieldTraits,
    ParticleGhostTraits  > ParticleField;


} // namespace Particle
} // namespace SpatialOps

#endif // ParticleFieldTypes_h
