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
  std::vector<int> dim(3,1);
  dim[0] = 10;
  std::vector<double> coordvec;

  IntVec totDim(&dim[0]);
  for( size_t i=0; i<3; ++i ) 
    if( dim[i]>1 ) totDim[i] += 2*CellField::Ghost::NGHOST;
  //
  // build the fields
  //
  CellField cellfiled( totDim, NULL );
  CellField   ctmp( totDim, NULL );

  SpatialOps::Particle::ParticleField pCoord( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField pSize  ( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField pfield( IntVec(np,1,1), NULL );
  SpatialOps::Particle::ParticleField ptmp  ( IntVec(np,1,1), NULL );
   

  //
  // set the cCoord coordinates.  These go from -0.5 to 10.5
  //
  CellField::iterator icoord = cellfiled.begin();
  CellField::iterator icoordEnd = cellfiled.end();
  for( size_t i=0; icoord != icoordEnd; ++icoord, ++i ){
    //cellfiled[i] = i*2;
    *icoord = i*2 ;
    coordvec.push_back(i*dx - dx/2.0);    
    std::cout<<"in for"<<std::endl;
  }
  /*
  for( size_t i=0; i<cellfiled.get_ntotal(); ++i )
    std::cout<<" Gas coord : " <<coordvec[i]<< "  cell field value : "<<cellfiled[i]<<std::endl;
*/

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
  const C2P* const c2p = new C2P( coordvec );
  typedef SpatialOps::Particle::ParticleToCell<CellField> P2C;
  const P2C* const p2c = new P2C( coordvec );

 //
  // interpolate to particles
  //
  c2p->apply_to_field( pCoord, cellfiled, ptmp );
  p2c->apply_to_field( pCoord,pSize, pfield, ctmp ); 
 

  for( size_t i=0; i<np; ++i )
    std::cout<<"  Interpolated particle field : "<<ptmp[i]<<std::endl;

  icoord = ctmp.begin();
  icoordEnd = ctmp.end();
  for( size_t i=0; icoord != icoordEnd; ++icoord, ++i )
    std::cout<<"  Interpolated Gas field : "<< *icoord <<std::endl;

  TestHelper status;
  status( ptmp[0] == 3 );
  status( ptmp[1] == 5.5 );
  status( ptmp[2] == 3);
  status( ptmp[3] == 16 );


  return 0;
}
