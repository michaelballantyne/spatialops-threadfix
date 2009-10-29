#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>

#include "TestHelper.h"


int main()
{
  using namespace SpatialOps;
  using namespace structured;

  std::vector<int> n(3,10);

  bool hasPlusFace[] = {true, true, true};

  TestHelper status(true);

  const size_t xdir = SpatialOps::XDIR::value;
  const size_t ydir = SpatialOps::YDIR::value;
  const size_t zdir = SpatialOps::ZDIR::value;

  status( npts<SVolField  >( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "SVol nx" );
  status( npts<SSurfXField>( xdir, n, hasPlusFace[xdir] ) == n[0]+3, "SSX  nx" );
  status( npts<SSurfYField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "SSY  nx" );
  status( npts<SSurfZField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "SSZ  nx" );

  status( npts<SVolField  >( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "SVol ny" );
  status( npts<SSurfXField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "SSX  ny" );
  status( npts<SSurfYField>( ydir, n, hasPlusFace[ydir] ) == n[1]+3, "SSY  ny" );
  status( npts<SSurfZField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "SSZ  ny" );

  status( npts<SVolField  >( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "SVol nz" );
  status( npts<SSurfXField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "SSX  nz" );
  status( npts<SSurfYField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "SSY  nz" );
  status( npts<SSurfZField>( zdir, n, hasPlusFace[zdir] ) == n[2]+3, "SSZ  nz" );


  status( npts<XVolField  >( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "XVol nx" );
  status( npts<XSurfXField>( xdir, n, hasPlusFace[xdir] ) == n[0]+3, "XSX  nx" );
  status( npts<XSurfYField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "XSY  nx" );
  status( npts<XSurfZField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "XSZ  nx" );

  status( npts<XVolField  >( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "XVol ny" );
  status( npts<XSurfXField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "XSX  ny" );
  status( npts<XSurfYField>( ydir, n, hasPlusFace[ydir] ) == n[1]+3, "XSY  ny" );
  status( npts<XSurfZField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "XSZ  ny" );

  status( npts<XVolField  >( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "XVol nz" );
  status( npts<XSurfXField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "XSX  nz" );
  status( npts<XSurfYField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "XSY  nz" );
  status( npts<XSurfZField>( zdir, n, hasPlusFace[zdir] ) == n[2]+3, "XSZ  nz" );


  status( npts<YVolField  >( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "YVol nx" );
  status( npts<YSurfXField>( xdir, n, hasPlusFace[xdir] ) == n[0]+3, "YSX  nx" );
  status( npts<YSurfYField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "YSY  nx" );
  status( npts<YSurfZField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "YSZ  nx" );

  status( npts<YVolField  >( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "YVol ny" );
  status( npts<YSurfXField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "YSX  ny" );
  status( npts<YSurfYField>( ydir, n, hasPlusFace[ydir] ) == n[1]+3, "YSY  ny" );
  status( npts<YSurfZField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "YSZ  ny" );

  status( npts<YVolField  >( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "YVol nz" );
  status( npts<YSurfXField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "YSX  nz" );
  status( npts<YSurfYField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "YSY  nz" );
  status( npts<YSurfZField>( zdir, n, hasPlusFace[zdir] ) == n[2]+3, "YSZ  nz" );


  status( npts<ZVolField  >( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "ZVol nx" );
  status( npts<ZSurfXField>( xdir, n, hasPlusFace[xdir] ) == n[0]+3, "ZSX  nx" );
  status( npts<ZSurfYField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "ZSY  nx" );
  status( npts<ZSurfZField>( xdir, n, hasPlusFace[xdir] ) == n[0]+2, "ZSZ  nx" );

  status( npts<ZVolField  >( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "ZVol ny" );
  status( npts<ZSurfXField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "ZSX  ny" );
  status( npts<ZSurfYField>( ydir, n, hasPlusFace[ydir] ) == n[1]+3, "ZSY  ny" );
  status( npts<ZSurfZField>( ydir, n, hasPlusFace[ydir] ) == n[1]+2, "ZSZ  ny" );

  status( npts<ZVolField  >( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "ZVol nz" );
  status( npts<ZSurfXField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "ZSX  nz" );
  status( npts<ZSurfYField>( zdir, n, hasPlusFace[zdir] ) == n[2]+2, "ZSY  nz" );
  status( npts<ZSurfZField>( zdir, n, hasPlusFace[zdir] ) == n[2]+3, "ZSZ  nz" );

  if( status.ok() ) return 0;
  return -1;
}
