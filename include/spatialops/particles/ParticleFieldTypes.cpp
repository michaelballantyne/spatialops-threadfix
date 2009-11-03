#include <spatialops/particles/ParticleFieldTypes.h>

namespace Particle{

  ParticleLinAlg::ParticleLinAlg()
    : vec(NULL), length(0)
  {}

  ParticleLinAlg::VecType&
  ParticleLinAlg::setup_vector( const size_t len,
                                double*const f )
  {
    vec = f;
    length = len;
    return vec;
  }

  void
  ParticleLinAlg::destroy_vector()
  {}

  ParticleLinAlg::~ParticleLinAlg()
  {
    destroy_vector();
  }

  void
  ParticleLinAlg::print_vec( std::ostream& s ) const
  {
    for( size_t i=0; i<length; ++i ){
      s << vec[i] << std::endl;
    }
  }

} // namespace Particle
