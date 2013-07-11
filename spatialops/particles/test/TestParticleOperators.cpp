#include <vector>
#include <set>
#include <iostream>
using std::cout;
using std::endl;


#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/MemoryWindow.h>

#include <spatialops/particles/ParticleFieldTypes.h>
#include <spatialops/particles/ParticleOperators.h>

#include <spatialops/WriteMatlab.h>

#include <test/TestHelper.h>

typedef SpatialOps::structured::SVolField CellField;
using SpatialOps::structured::IntVec;
using SpatialOps::write_matlab;
using SpatialOps::structured::MemoryWindow;
namespace SS=SpatialOps::structured;

int main()
{
  const size_t np=1;
  //const double dx = 1.0;

  IntVec totDim(10,1,1);
  for( size_t i=0; i<3; ++i )
    if( totDim[i]>1 )
      totDim[i] += CellField::Ghost::NGhostMinus::int_vec()[i] +
                   CellField::Ghost::NGhostPlus::int_vec()[i];

  const SS::GhostDataRT cg(1);
  const SS::BoundaryCellInfo cbc = SS::BoundaryCellInfo::build<CellField>();
  const MemoryWindow mw( totDim );

  //
  // build the fields
  //
  CellField       cellField( mw, cbc, cg, NULL );
  CellField cellFieldvalues( mw, cbc, cg, NULL );
  CellField            ctmp( mw, cbc, cg, NULL );

  const SS::GhostDataRT pg(0);
  const SS::BoundaryCellInfo pbc = SS::BoundaryCellInfo::build<SpatialOps::Particle::ParticleField>();
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
    *fielditer = i * 10 ;
  }

  //
  // set the particle coordinates
  //
  pCoord[0] = 3;/*
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
    std::cout<<"  particle coord : "<<pCoord[i]<<"  particle field : "<<pfield[i]<<std::endl;

  //
  // build the operators
  //
  typedef SpatialOps::Particle::CellToParticle<CellField> C2P;
  const C2P c2p( cellField );
  typedef SpatialOps::Particle::ParticleToCell<CellField> P2C;
  const P2C p2c( cellField );

  //
  // interpolate to particles
  //
  c2p.apply_to_field( pCoord, pSize, cellFieldvalues, ptmp );
  p2c.apply_to_field( pCoord, pSize, pfield, ctmp );

  for( size_t i=0; i<np; ++i )
    std::cout<<"  Interpolated particle field : "<<ptmp[i]<<std::endl;

  CellField::const_iterator ix=cellField.begin();
  for( CellField::const_iterator i=ctmp.begin(); i!=ctmp.end(); ++i, ++ix )
    std::cout<<"  Interpolated Gas field at x=" << *ix << " = "<< *i <<std::endl;

  TestHelper status;
  for( size_t i=0; i<np; ++i )
    status( ptmp[i] == 35 );

  status( ctmp[0] == 0, "gas interp [0]" );
  status( ctmp[1] == 2, "gas interp [1]" );
  status( ctmp[2] ==  4, "gas interp [2]" );
  status( ctmp[3] ==  4, "gas interp [3]" );
  status( ctmp[4] == 4, "gas interp [4]" );
  status( ctmp[5] ==  4, "gas interp [5]" );
  status( ctmp[6] ==  2, "gas interp [6]" );

  if( status.ok() ) return 0;
  return -1;
}
