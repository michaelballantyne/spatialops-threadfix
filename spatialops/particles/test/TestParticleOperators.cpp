#include <vector>
#include <set>
#include <iostream>
using std::cout;
using std::endl;


#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

#include <spatialops/particles/ParticleFieldTypes.h>
#include <spatialops/particles/ParticleOperators.h>

#include <spatialops/WriteMatlab.h>

#include <test/TestHelper.h>

typedef SpatialOps::structured::SVolField CellField;
using SpatialOps::structured::IntVec;
using SpatialOps::write_matlab;

int main()
{
  const size_t np=4;
  const double dx = 1.0;

  IntVec totDim(10,1,1);
  for( size_t i=0; i<3; ++i ) 
    if( totDim[i]>1 ) totDim[i] += 2*CellField::Ghost::NGHOST;

  //
  // build the fields
  //
  CellField cellField( totDim, NULL );
  CellField   ctmp( totDim, NULL );

  SpatialOps::Particle::ParticleField pCoord( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField pSize ( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField pfield( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField ptmp  ( IntVec(np,1,1), NULL );
   
  //
  // set the cCoord coordinates.  These go from -0.5 to 10.5
  //
  
  size_t i=0;
  for( CellField::iterator icoord=cellField.begin(); icoord!=cellField.end(); ++icoord, ++i ){
    *icoord = i*2 ;
  }

  //
  // set the particle coordinates
  //
  pCoord[0] = 1;
  pCoord[1] = 2.25;
  pCoord[2] = 1;
  pCoord[3] = 7.5;

  pSize[0] = 0;
  pSize[1] = 0;
  pSize[2] = 0;
  pSize[3] = 0;

  pfield[0] = 10;
  pfield[1] = 20;
  pfield[2] = 40;
  pfield[3] = 70;

  for( size_t i=0; i<np; ++i )
    std::cout<<"  particle coord : "<<pCoord[i]<<"  particle field : "<<pfield[i]<<std::endl;

  //
  // build the operator
  //
  typedef SpatialOps::Particle::CellToParticle<CellField> C2P;
  const C2P* const c2p = new C2P( cellField );
  typedef SpatialOps::Particle::ParticleToCell<CellField> P2C;
  const P2C* const p2c = new P2C( cellField );

  //
  // interpolate to particles
  //
  c2p->apply_to_field( pCoord, cellField, ptmp );
  p2c->apply_to_field( pCoord, pSize, pfield, ctmp ); 
 
  for( size_t i=0; i<np; ++i )
    std::cout<<"  Interpolated particle field : "<<ptmp[i]<<std::endl;

  CellField::const_iterator ix=cellField.begin();
  for( CellField::const_iterator i=ctmp.begin(); i!=ctmp.end(); ++i, ++ix )
    std::cout<<"  Interpolated Gas field at x=" << *ix << " = "<< *i <<std::endl;

  TestHelper status;
  for( size_t i=0; i<np; ++i )
    status( ptmp[i] == pCoord[i] );

  status( ctmp[0] == 50, "gas interp [0]" );
  status( ctmp[1] == 20, "gas interp [1]" );
  status( ctmp[2] ==  0, "gas interp [2]" );
  status( ctmp[3] ==  0, "gas interp [3]" );
  status( ctmp[4] == 70, "gas interp [4]" );
  status( ctmp[5] ==  0, "gas interp [5]" );
  status( ctmp[6] ==  0, "gas interp [6]" );

  if( status.ok() ) return 0;
  return -1;
}
