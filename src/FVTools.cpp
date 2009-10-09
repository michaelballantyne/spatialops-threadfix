#include <spatialops/FVTools.h>

namespace SpatialOps{
namespace FVStaggered{

  void _ghost_set_( const int ngm, const int ngp,
                    const int nxt, const int nyt, const int nzt,
                    const std::vector<int>& dim,
                    const bool hasPlusXSideFaces,
                    const bool hasPlusYSideFaces,
                    const bool hasPlusZSideFaces,
                    int& ix,
                    std::set<int>& ghostSet )
  {
    const int ngxm = dim[0]>1 ? ngm : 0;
    const int ngxp = dim[0]>1 ? ngp : 0;
    const int ngym = dim[1]>1 ? ngm : 0;
    const int ngyp = dim[1]>1 ? ngp : 0;
    const int ngzm = dim[2]>1 ? ngm : 0;
    const int ngzp = dim[2]>1 ? ngp : 0;

    // -z side ghost layer
    if( dim[2]>1 ){
      for( int kg=0; kg<ngzm; ++kg )
        for( int j=0; j<nyt; ++j )
          for( int i=0; i<nxt; ++i )
            ghostSet.insert(ix++);
    }

    // z interior
    for( int k=ngzm; k<nzt-ngzp; ++k ){

      // -y side ghost layer
      if( dim[1]>1 ){
        for( int i=0; i<nxt; ++i )
          for( int jg=0; jg<ngym; ++jg )
            ghostSet.insert(ix++);
      }

      // y interior
      for( int j=ngym; j<nyt-ngyp; ++j ){
        // -x side ghost layer
        if( dim[0]>1 ) for( int ig=0; ig<ngxm; ++ig ) ghostSet.insert(ix++);
        // x interior
        ix+=nxt-ngxm-ngxp;
        // +x side ghost layer
        if( dim[0]>1 ) for( int ig=0; ig<ngxp; ++ig ) ghostSet.insert(ix++);
      }

      // +y side ghost layer
      if( dim[1]>1 ){
        for( int i=0; i<nxt; ++i )
          for( int jg=0; jg<ngyp; ++jg )
            ghostSet.insert(ix++);
      }
    }

    // +z side ghost layer
    if( dim[2]>1 ){
      for( int kg=0; kg<ngzp; ++kg )
        for( int i=0; i<nxt; ++i )
          for( int j=0; j<nyt; ++j )
            ghostSet.insert(ix++);
    }
  }

} // namespace FVStaggered
} // namespace SpatialOps
