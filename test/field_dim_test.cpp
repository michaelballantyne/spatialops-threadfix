#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

#include "TestHelper.h"


int main()
{
  using namespace SpatialOps;
  using namespace structured;

  IntVec n(10,10,10);

  bool hasPlusFace[] = {true, true, true};

  TestHelper status(true);

  status( get_nx_with_ghost<SVolField  >( n[0], hasPlusFace[0] ) == n[0]+2, "SVol nx" );
  status( get_nx_with_ghost<SSurfXField>( n[0], hasPlusFace[0] ) == n[0]+3, "SSX  nx" );
  status( get_nx_with_ghost<SSurfYField>( n[0], hasPlusFace[0] ) == n[0]+2, "SSY  nx" );
  status( get_nx_with_ghost<SSurfZField>( n[0], hasPlusFace[0] ) == n[0]+2, "SSZ  nx" );

  status( get_ny_with_ghost<SVolField  >( n[1], hasPlusFace[1] ) == n[1]+2, "SVol ny" );
  status( get_ny_with_ghost<SSurfXField>( n[1], hasPlusFace[1] ) == n[1]+2, "SSX  ny" );
  status( get_ny_with_ghost<SSurfYField>( n[1], hasPlusFace[1] ) == n[1]+3, "SSY  ny" );
  status( get_ny_with_ghost<SSurfZField>( n[1], hasPlusFace[1] ) == n[1]+2, "SSZ  ny" );

  status( get_nz_with_ghost<SVolField  >( n[2], hasPlusFace[2] ) == n[2]+2, "SVol nz" );
  status( get_nz_with_ghost<SSurfXField>( n[2], hasPlusFace[2] ) == n[2]+2, "SSX  nz" );
  status( get_nz_with_ghost<SSurfYField>( n[2], hasPlusFace[2] ) == n[2]+2, "SSY  nz" );
  status( get_nz_with_ghost<SSurfZField>( n[2], hasPlusFace[2] ) == n[2]+3, "SSZ  nz" );


  status( get_nx_with_ghost<XVolField  >( n[0], hasPlusFace[0] ) == n[0]+2, "XVol nx" );
  status( get_nx_with_ghost<XSurfXField>( n[0], hasPlusFace[0] ) == n[0]+3, "XSX  nx" );
  status( get_nx_with_ghost<XSurfYField>( n[0], hasPlusFace[0] ) == n[0]+2, "XSY  nx" );
  status( get_nx_with_ghost<XSurfZField>( n[0], hasPlusFace[0] ) == n[0]+2, "XSZ  nx" );

  status( get_ny_with_ghost<XVolField  >( n[1], hasPlusFace[1] ) == n[1]+2, "XVol ny" );
  status( get_ny_with_ghost<XSurfXField>( n[1], hasPlusFace[1] ) == n[1]+2, "XSX  ny" );
  status( get_ny_with_ghost<XSurfYField>( n[1], hasPlusFace[1] ) == n[1]+3, "XSY  ny" );
  status( get_ny_with_ghost<XSurfZField>( n[1], hasPlusFace[1] ) == n[1]+2, "XSZ  ny" );

  status( get_nz_with_ghost<XVolField  >( n[2], hasPlusFace[2] ) == n[2]+2, "XVol nz" );
  status( get_nz_with_ghost<XSurfXField>( n[2], hasPlusFace[2] ) == n[2]+2, "XSX  nz" );
  status( get_nz_with_ghost<XSurfYField>( n[2], hasPlusFace[2] ) == n[2]+2, "XSY  nz" );
  status( get_nz_with_ghost<XSurfZField>( n[2], hasPlusFace[2] ) == n[2]+3, "XSZ  nz" );


  status( get_nx_with_ghost<YVolField  >( n[0], hasPlusFace[0] ) == n[0]+2, "YVol nx" );
  status( get_nx_with_ghost<YSurfXField>( n[0], hasPlusFace[0] ) == n[0]+3, "YSX  nx" );
  status( get_nx_with_ghost<YSurfYField>( n[0], hasPlusFace[0] ) == n[0]+2, "YSY  nx" );
  status( get_nx_with_ghost<YSurfZField>( n[0], hasPlusFace[0] ) == n[0]+2, "YSZ  nx" );

  status( get_ny_with_ghost<YVolField  >( n[1], hasPlusFace[1] ) == n[1]+2, "YVol ny" );
  status( get_ny_with_ghost<YSurfXField>( n[1], hasPlusFace[1] ) == n[1]+2, "YSX  ny" );
  status( get_ny_with_ghost<YSurfYField>( n[1], hasPlusFace[1] ) == n[1]+3, "YSY  ny" );
  status( get_ny_with_ghost<YSurfZField>( n[1], hasPlusFace[1] ) == n[1]+2, "YSZ  ny" );

  status( get_nz_with_ghost<YVolField  >( n[2], hasPlusFace[2] ) == n[2]+2, "YVol nz" );
  status( get_nz_with_ghost<YSurfXField>( n[2], hasPlusFace[2] ) == n[2]+2, "YSX  nz" );
  status( get_nz_with_ghost<YSurfYField>( n[2], hasPlusFace[2] ) == n[2]+2, "YSY  nz" );
  status( get_nz_with_ghost<YSurfZField>( n[2], hasPlusFace[2] ) == n[2]+3, "YSZ  nz" );


  status( get_nx_with_ghost<ZVolField  >( n[0], hasPlusFace[0] ) == n[0]+2, "ZVol nx" );
  status( get_nx_with_ghost<ZSurfXField>( n[0], hasPlusFace[0] ) == n[0]+3, "ZSX  nx" );
  status( get_nx_with_ghost<ZSurfYField>( n[0], hasPlusFace[0] ) == n[0]+2, "ZSY  nx" );
  status( get_nx_with_ghost<ZSurfZField>( n[0], hasPlusFace[0] ) == n[0]+2, "ZSZ  nx" );

  status( get_ny_with_ghost<ZVolField  >( n[1], hasPlusFace[1] ) == n[1]+2, "ZVol ny" );
  status( get_ny_with_ghost<ZSurfXField>( n[1], hasPlusFace[1] ) == n[1]+2, "ZSX  ny" );
  status( get_ny_with_ghost<ZSurfYField>( n[1], hasPlusFace[1] ) == n[1]+3, "ZSY  ny" );
  status( get_ny_with_ghost<ZSurfZField>( n[1], hasPlusFace[1] ) == n[1]+2, "ZSZ  ny" );

  status( get_nz_with_ghost<ZVolField  >( n[2], hasPlusFace[2] ) == n[2]+2, "ZVol nz" );
  status( get_nz_with_ghost<ZSurfXField>( n[2], hasPlusFace[2] ) == n[2]+2, "ZSX  nz" );
  status( get_nz_with_ghost<ZSurfYField>( n[2], hasPlusFace[2] ) == n[2]+2, "ZSY  nz" );
  status( get_nz_with_ghost<ZSurfZField>( n[2], hasPlusFace[2] ) == n[2]+3, "ZSZ  nz" );

  if( status.ok() ) return 0;
  return -1;
}
