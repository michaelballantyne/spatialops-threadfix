#include <vector>
#include <set>
#include <iostream>
using std::cout;
using std::endl;


#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>

#include <spatialops/particles/ParticleFieldTypes.h>
#include <spatialops/particles/ParticleOperators.h>

#include "TestHelper.h"

typedef SpatialOps::structured::SVolField CellField;

int main()
{
  size_t np=4;
  double dx = 1.0;
  std::vector<int> dim(3,1);
  dim[0] = 10;

  //
  // build the fields
  //
  CellField cCoord( SpatialOps::structured::npts<CellField>(SpatialOps::XDIR::value,dim,true),
                    SpatialOps::structured::get_ghost_set<CellField>( dim, true, true, true ),
                    NULL );

  CellField ctmp( SpatialOps::structured::npts<CellField>(SpatialOps::XDIR::value,dim,true),
                  SpatialOps::structured::get_ghost_set<CellField>( dim, true, true, true ),
                  NULL );

  SpatialOps::Particle::ParticleField pCoord( np, std::set<int>(), NULL );
  SpatialOps::Particle::ParticleField ptmp  ( np, std::set<int>(), NULL );
                  
  //
  // set the cCoord coordinates.  These go from -0.5 to 10.5
  //
  for( size_t i=0; i<cCoord.get_ntotal(); ++i ){
    cCoord[i] = i*dx - dx/2.0;
    ctmp[i] = -double(i)*dx;
  }

  //
  // set the particle coordinates
  //
  pCoord[0] = 1;
  pCoord[1] = 2.25;
  pCoord[2] = 4.1;
  pCoord[3] = 7.5;

  //
  // build the operator
  //
  typedef SpatialOps::Particle::CellToParticle<CellField> C2P;
  const C2P* const c2p = new C2P( pCoord, cCoord );

  //
  // interpolate to particles
  //
  c2p->apply_to_field( ctmp, ptmp );

  TestHelper status;

  status( ptmp[0] == -1.50 );
  status( ptmp[1] == -2.75 );
  status( ptmp[2] == -4.60 );
  status( ptmp[3] == -8.00 );

  if( status.isfailed() ){
    cCoord.write_matlab("x");
    ctmp.write_matlab("ctmp");
    ptmp.write_matlab("pmesh");
    pCoord.write_matlab("px");
    return -1;
  }

  return 0;
}
