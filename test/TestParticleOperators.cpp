#include <vector>
#include <set>
#include <iostream>
using std::cout;
using std::endl;


#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/MemoryWindow.h>

#include <spatialops/particles/ParticleFieldTypes.h>
#include <spatialops/particles/ParticleOperators.h>

#include <spatialops/WriteMatlab.h>

#include "TestHelper.h"

typedef SpatialOps::structured::SVolField CellField;
using SpatialOps::structured::IntVec;
using SpatialOps::write_matlab;

int main()
{
  const size_t np=4;
  const double dx = 1.0;
  std::vector<int> dim(3,1);
  dim[0] = 10;

  IntVec totDim(&dim[0]);
  for( size_t i=0; i<3; ++i ) 
    if( dim[i]>1 ) totDim[i] += 2*CellField::Ghost::NGHOST;
  //
  // build the fields
  //
  CellField cCoord( totDim, NULL );
  CellField   ctmp( totDim, NULL );

  SpatialOps::Particle::ParticleField pCoord( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField ptmp  ( IntVec(np,1,1), NULL );
                  
  //
  // set the cCoord coordinates.  These go from -0.5 to 10.5
  //
  CellField::iterator icoord=cCoord.begin(), ictmp=ctmp.begin();
  CellField::iterator icoordEnd=cCoord.end();
  int i=0;
  for( ; icoord!=icoordEnd; ++icoord, ++ictmp, ++i ){
    *icoord = i*dx - dx/2.0;
    *ictmp = -double(i)*dx;
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
  const C2P* const c2p = new C2P( cCoord );

  //
  // interpolate to particles
  //
  c2p->apply_to_field( pCoord, ctmp, ptmp );

  TestHelper status;

  status( ptmp[0] == -1.50 );
  status( ptmp[1] == -2.75 );
  status( ptmp[2] == -4.60 );
  status( ptmp[3] == -8.00 );

  if( status.isfailed() ){
    write_matlab(cCoord,"x");
    write_matlab(ctmp,"ctmp");
    write_matlab(ptmp,"pmesh");
    write_matlab(pCoord,"px");
    return -1;
  }

  return 0;
}
