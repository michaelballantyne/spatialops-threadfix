#include <vector>
#include <set>
#include <iostream>
using std::cout;
using std::endl;

#include <spatialops/structured/FVStaggeredFieldTypes.h>

#include <spatialops/particles/ParticleFieldTypes.h>
#include <spatialops/particles/ParticleOperators.h>

#include <spatialops/WriteMatlab.h>

#include <test/TestHelper.h>

typedef SpatialOps::SVolField CellField;
using namespace SpatialOps;

int main()
{
  const size_t np=1;
  //const double dx = 1.0;

  IntVec dim(10,1,1);

  const GhostData cg(1);
  const BoundaryCellInfo cbc = BoundaryCellInfo::build<CellField>();
  const MemoryWindow mw = get_window_with_ghost( dim, cg, BoundaryCellInfo::build<CellField>(false,false,false) );

  //
  // build the fields
  //
  CellField       cellField( mw, cbc, cg, NULL );
  CellField cellFieldvalues( mw, cbc, cg, NULL );
  CellField            ctmp( mw, cbc, cg, NULL );

  const GhostData pg(0);
  const BoundaryCellInfo pbc = BoundaryCellInfo::build<SpatialOps::Particle::ParticleField>();
  const MemoryWindow pmw( IntVec(np,1,1) );

  SpatialOps::Particle::ParticleField pCoord( pmw, pbc, pg, NULL );
  SpatialOps::Particle::ParticleField pSize ( pmw, pbc, pg, NULL );
  SpatialOps::Particle::ParticleField pfield( pmw, pbc, pg, NULL );
  SpatialOps::Particle::ParticleField ptmp  ( pmw, pbc, pg, NULL );

  //
  // set the cCoord coordinates.  These go from -0.5 to 10.5
  //

  size_t i=0;
  CellField::iterator fielditer=cellFieldvalues.begin();
  for( CellField::iterator icoord=cellField.begin(); icoord!=cellField.end(); ++icoord, ++i, ++fielditer ){
    *icoord = i - 0.5 ;
    *fielditer = *icoord * 10 ;
  }

  //
  // set the particle coordinates
  //
  pCoord[0] = 3.25;/*
  pCoord[1] = 7.0;
  pCoord[2] = 1;
  pCoord[3] = 7.5;*/

  pSize[0] = 5;/*
  pSize[1] = 2;
  pSize[2] = 0;
  pSize[3] = 0;*/

  pfield[0] = 20;/*
  pfield[1] = 20;
  pfield[2] = 40;
  pfield[3] = 70;*/

  for( size_t i=0; i<np; ++i )
    std::cout << "  particle coord : " << pCoord[i] << "  particle field : " << pfield[i] << std::endl;

  //
  // build the operators
  //
  typedef SpatialOps::Particle::CellToParticle<CellField> C2P;
  C2P c2p( cellField[1]-cellField[0], cellField[0] );
  c2p.set_coordinate_information( &pCoord, NULL, NULL, &pSize );

  typedef SpatialOps::Particle::ParticleToCell<CellField> P2C;
  P2C p2c( cellField[1] - cellField[0], cellField[0] );
  p2c.set_coordinate_information( &pCoord, NULL, NULL, &pSize );

  //
  // interpolate to particles
  //
  c2p.apply_to_field( cellFieldvalues, ptmp );
  p2c.apply_to_field( pfield, ctmp );

  for( size_t i=0; i<np; ++i )
    std::cout << "  Interpolated gas value to particle field : " << ptmp[i] << std::endl;

  CellField::const_iterator ix=cellField.begin();
  for( CellField::const_iterator i=ctmp.begin(); i!=ctmp.end(); ++i, ++ix )
    std::cout << "  Interpolated particle value to gas field at x=" << *ix << " = "<< *i << std::endl;

  TestHelper status;
  for( size_t i=0; i<np; ++i )
    status( ptmp[i] == 32.5, "c2p" );

  status( ctmp[0] == 0, "p2c [0]" );
  status( ctmp[1] == 1, "p2c [1]" );
  status( ctmp[2] == 4, "p2c [2]" );
  status( ctmp[3] == 4, "p2c [3]" );
  status( ctmp[4] == 4, "p2c [4]" );
  status( ctmp[5] == 4, "p2c [5]" );
  status( ctmp[6] == 3, "p2c [6]" );
  status( ctmp[7] == 0, "p2c [7]" );
  status( ctmp[8] == 0, "p2c [8]" );

  if( status.ok() ) return 0;
  return -1;
}
