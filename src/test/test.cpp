#include <vector>

#include <FV2ndOrderTypes.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>

using namespace std;

#define PI 3.14159265358979

void setup_geom( const std::vector<int> & dim,
		 std::vector<double> & spacing,
		 std::vector<double> & area,
		 double & volume )
{
  for( int i=0; i<3; ++i )
    spacing[i] = (dim[i]==1) ? 1.0 : PI/(dim[i]-1);

  area[0] = spacing[1]*spacing[2];
  area[1] = spacing[0]*spacing[2];
  area[2] = spacing[0]*spacing[1];

  volume = spacing[0]*spacing[1]*spacing[2];
}

//--------------------------------------------------------------------

void build_ops( const std::vector<int> & dim )
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  std::vector<double> spacing(3), area(3);
  double volume;
  setup_geom( dim, spacing, area, volume );

  //
  // build the assemblers
  //
  InterpXC2F::Assembler rxcfasmbl( dim );
  InterpYC2F::Assembler rycfasmbl( dim );
  InterpZC2F::Assembler rzcfasmbl( dim );

  InterpXF2C::Assembler rxfcasmbl( dim );
  InterpYF2C::Assembler ryfcasmbl( dim );
  InterpZF2C::Assembler rzfcasmbl( dim );

  GradXC2F::Assembler gxcfasmbl( spacing, dim );
  GradYC2F::Assembler gycfasmbl( spacing, dim );
  GradZC2F::Assembler gzcfasmbl( spacing, dim );

  GradXF2C::Assembler gxfcasmbl( spacing, dim );
  GradYF2C::Assembler gyfcasmbl( spacing, dim );
  GradZF2C::Assembler gzfcasmbl( spacing, dim );

  DivXF2C::Assembler dxfcasmbl( dim, area, volume );
  DivYF2C::Assembler dyfcasmbl( dim, area, volume );
  DivZF2C::Assembler dzfcasmbl( dim, area, volume );

  InterpX_YF2ZE::Assembler rx_yfze_a( dim );
  InterpX_ZF2YE::Assembler rx_zfye_a( dim );
  InterpY_XF2ZE::Assembler ry_xfze_a( dim );
  InterpY_ZF2XE::Assembler ry_zfxe_a( dim );
  InterpZ_XF2YE::Assembler rz_xfye_a( dim );
  InterpZ_YF2XE::Assembler rz_yfxe_a( dim );

  GradX_YF2ZE::Assembler gx_yfze_a( spacing, dim );
  GradX_ZF2YE::Assembler gx_zfye_a( spacing, dim );
  GradY_XF2ZE::Assembler gy_xfze_a( spacing, dim );
  GradY_ZF2XE::Assembler gy_zfxe_a( spacing, dim );
  GradZ_XF2YE::Assembler gz_xfye_a( spacing, dim );
  GradZ_YF2XE::Assembler gz_yfxe_a( spacing, dim );

  DivX_YE2ZF::Assembler dx_yfze_a( dim, area, volume );
  DivX_ZE2YF::Assembler dx_zfye_a( dim, area, volume );
  DivY_XE2ZF::Assembler dy_xfze_a( dim, area, volume );
  DivY_ZE2XF::Assembler dy_zfxe_a( dim, area, volume );
  DivZ_XE2YF::Assembler dz_xfye_a( dim, area, volume );
  DivZ_YE2XF::Assembler dz_yfxe_a( dim, area, volume );

  SxCell::Assembler  sxcasmbl( dim );
  SyCell::Assembler  sycasmbl( dim );
  SzCell::Assembler  szcasmbl( dim );

  SxCellSide::Assembler  sxcfasmbl( dim );
  SyCellSide::Assembler  sycfasmbl( dim );
  SzCellSide::Assembler  szcfasmbl( dim );

  // 
  // Inerpolants
  //
  SpatialOpDatabase<InterpXC2F>::self().register_new_operator( new InterpXC2F( rxcfasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<InterpYC2F>::self().register_new_operator( new InterpYC2F( rycfasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<InterpZC2F>::self().register_new_operator( new InterpZC2F( rzcfasmbl ) );

  SpatialOpDatabase<InterpXF2C>::self().register_new_operator( new InterpXF2C( rxfcasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<InterpYF2C>::self().register_new_operator( new InterpYF2C( ryfcasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<InterpZF2C>::self().register_new_operator( new InterpZF2C( rzfcasmbl ) );

  SpatialOpDatabase<InterpX_YF2ZE>::self().register_new_operator( new InterpX_YF2ZE( rx_yfze_a ) );
  SpatialOpDatabase<InterpX_ZF2YE>::self().register_new_operator( new InterpX_ZF2YE( rx_zfye_a ) );
  if( dim[1]>1 ){
   SpatialOpDatabase<InterpY_XF2ZE>::self().register_new_operator( new InterpY_XF2ZE( ry_xfze_a ) );
   SpatialOpDatabase<InterpY_ZF2XE>::self().register_new_operator( new InterpY_ZF2XE( ry_zfxe_a ) );
  }
  if( dim[2]>1 ){
    SpatialOpDatabase<InterpZ_XF2YE>::self().register_new_operator( new InterpZ_XF2YE( rz_xfye_a ) );
    SpatialOpDatabase<InterpZ_YF2XE>::self().register_new_operator( new InterpZ_YF2XE( rz_yfxe_a ) );
  }

  //
  // Gradient
  //
  SpatialOpDatabase<GradXC2F>::self().register_new_operator( new GradXC2F( gxcfasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<GradYC2F>::self().register_new_operator( new GradYC2F( gycfasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<GradZC2F>::self().register_new_operator( new GradZC2F( gzcfasmbl ) );

  SpatialOpDatabase<GradXF2C>::self().register_new_operator( new GradXF2C( gxfcasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<GradYF2C>::self().register_new_operator( new GradYF2C( gyfcasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<GradZF2C>::self().register_new_operator( new GradZF2C( gzfcasmbl ) );

  SpatialOpDatabase<GradX_YF2ZE>::self().register_new_operator( new GradX_YF2ZE( gx_yfze_a ) );
  SpatialOpDatabase<GradX_ZF2YE>::self().register_new_operator( new GradX_ZF2YE( gx_zfye_a ) );
  if( dim[1]>1 ){
    SpatialOpDatabase<GradY_XF2ZE>::self().register_new_operator( new GradY_XF2ZE( gy_xfze_a ) );
    SpatialOpDatabase<GradY_ZF2XE>::self().register_new_operator( new GradY_ZF2XE( gy_zfxe_a ) );
  }
  if( dim[2]>1 ){
    SpatialOpDatabase<GradZ_XF2YE>::self().register_new_operator( new GradZ_XF2YE( gz_xfye_a ) );
    SpatialOpDatabase<GradZ_YF2XE>::self().register_new_operator( new GradZ_YF2XE( gz_yfxe_a ) );
  }

  //
  // Divergence
  //
  SpatialOpDatabase<DivXF2C>::self().register_new_operator( new DivXF2C( dxfcasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<DivYF2C>::self().register_new_operator( new DivYF2C( dyfcasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<DivZF2C>::self().register_new_operator( new DivZF2C( dzfcasmbl ) );

  SpatialOpDatabase<DivX_YE2ZF>::self().register_new_operator( new DivX_YE2ZF( dx_yfze_a ) );
  SpatialOpDatabase<DivX_ZE2YF>::self().register_new_operator( new DivX_ZE2YF( dx_zfye_a ) );
  if( dim[1]>1 ){
    SpatialOpDatabase<DivY_XE2ZF>::self().register_new_operator( new DivY_XE2ZF( dy_xfze_a ) );
    SpatialOpDatabase<DivY_ZE2XF>::self().register_new_operator( new DivY_ZE2XF( dy_zfxe_a ) );
  }
  if( dim[2]>1 ){
    SpatialOpDatabase<DivZ_XE2YF>::self().register_new_operator( new DivZ_XE2YF( dz_xfye_a ) );
    SpatialOpDatabase<DivZ_YE2XF>::self().register_new_operator( new DivZ_YE2XF( dz_yfxe_a ) );
  }

  //
  // Scratch
  //
  SpatialOpDatabase<SxCellSide>::self().register_new_operator( new SxCellSide( sxcfasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<SyCellSide>::self().register_new_operator( new SyCellSide( sycfasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<SzCellSide>::self().register_new_operator( new SzCellSide( szcfasmbl ) );

  SpatialOpDatabase<SxCell>::self().register_new_operator( new SxCell( sxcasmbl ) );
  if( dim[1]>1 ) SpatialOpDatabase<SyCell>::self().register_new_operator( new SyCell( sycasmbl ) );
  if( dim[2]>1 ) SpatialOpDatabase<SzCell>::self().register_new_operator( new SzCell( szcasmbl ) );

}



bool test( const std::vector<int> & dim )
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  std::vector<double> spacing(3), area(3);
  double volume;
  setup_geom( dim, spacing, area, volume );

  // get the operators
  InterpXC2F & rxcf = *SpatialOpDatabase<InterpXC2F>::self().retrieve_operator( dim );
  InterpXF2C & rxfc = *SpatialOpDatabase<InterpXF2C>::self().retrieve_operator( dim );
  GradXC2F   & gxcf = *SpatialOpDatabase<GradXC2F  >::self().retrieve_operator( dim );
  DivXF2C    & dxfc = *SpatialOpDatabase<DivXF2C   >::self().retrieve_operator( dim );
  SxCellSide & sxcf = *SpatialOpDatabase<SxCellSide>::self().retrieve_operator( dim );
  SxCell     & sxc  = *SpatialOpDatabase<SxCell    >::self().retrieve_operator( dim );

  /*
  EpetraExt::RowMatrixToMatrixMarketFile( "Dx.mm", dxfc.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Gx.mm", gxcf.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Rx.mm", rxcf.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Rxfc.mm", rxfc.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Sx.mm", sxc .get_linalg_mat(), "", "" );
  */

  // build the spatial fields
  CellField      x( dim, NULL, InternalStorage );
  CellField      f( dim, NULL, InternalStorage );
  XSideField  xint( dim, NULL, InternalStorage );
  XSideField  fint( dim, NULL, InternalStorage );
  XSideField gradF( dim, NULL, InternalStorage );
  CellField    d2f( dim, NULL, InternalStorage );
  CellField  fint2( dim, NULL, InternalStorage );


  const double dx = spacing[0]; //1.0/double(dim[0]-1);
  const int ngxlo = CellFieldTraits::GhostTraits::get<XDIR,SideMinus>();
  const int khi = dim[2]==1 ? 1 : dim[2]+CellFieldTraits::GhostTraits::get<ZDIR,SideMinus>() + CellFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
  const int jhi = dim[1]==1 ? 1 : dim[1]+CellFieldTraits::GhostTraits::get<YDIR,SideMinus>() + CellFieldTraits::GhostTraits::get<YDIR,SidePlus>();
  const int ihi = dim[0]+CellFieldTraits::GhostTraits::get<XDIR,SideMinus>() + CellFieldTraits::GhostTraits::get<XDIR,SidePlus>();
  const int ilo=0, jlo=0, klo=0;
  int ix = 0;
  for( int k=klo; k<khi ; ++k ){
    for( int j=jlo; j<jhi; ++j ){
      for( int i=ilo; i<ihi; ++i ){
	x[ix] = (i-ngxlo)*dx;
	f[ix] = std::sin(x[ix]) + 2.0;
	++ix;
      }
    }
  }

  // apply the operators to the fields
  gxcf.apply_to_field( f, gradF );
  rxcf.apply_to_field( x, xint  );
  rxcf.apply_to_field( f, fint  );

  rxfc.apply_to_field( fint, fint2 );

  dxfc.apply_to_op( gxcf, sxc );
  //  EpetraExt::RowMatrixToMatrixMarketFile( "Lx.mm", sxc.get_linalg_mat(), "", "" );

  sxc.apply_to_field( f, d2f );

  //gxfc.apply_to_field( gradF, d2f );


  sxcf  = gxcf;
  sxcf += gxcf;
  sxcf -= gxcf;
  sxcf += fint;


  /*
  EpetraExt::VectorToMatrixMarketFile( "x.mm",   x.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "fx.mm",   f.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "xx.mm", xint.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "fintx.mm", fint.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "dfdx.mm", gradF.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "d2fdx2.mm", d2f.get_linalg_vec(), "", "" );
  */

  // build a field without any ghost information in it.
  CellFieldNoGhost tmpField( dim, NULL, InternalStorage );
  RHS tmp(dim);
  tmp.reset();
  tmp.add_field_contribution( x );
  tmpField = tmp;

  // this is really a hack, but it is a way to eliminate ghost values from the output...
  tmp.reset(); tmp.add_field_contribution(x    ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "x.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(f    ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(xint ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "xx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(fint ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fintx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(gradF); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "dfdx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(d2f  ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "d2fdx2.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(fint2); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fintx2.mm", tmpField.get_linalg_vec(), "", "" );

  tmpField += tmp;  // test += operator
  tmpField -= tmp;  // test -= operator
  CellFieldNoGhost tmp2( dim,NULL,InternalStorage );
  tmp2 = tmp;
  cout << endl << "SpatialField operators -=, +=, ==, = ... ";
  if( tmp2!=tmpField ) cout << "FAIL!" << std::endl;
  else cout << "PASS" << std::endl;

  if( dim[1] > 1 ){

    // get the operators
    InterpYC2F & rycf = *SpatialOpDatabase<InterpYC2F>::self().retrieve_operator( dim );
    InterpYF2C & ryfc = *SpatialOpDatabase<InterpYF2C>::self().retrieve_operator( dim );
    GradYC2F   & gycf = *SpatialOpDatabase<GradYC2F  >::self().retrieve_operator( dim );
    DivYF2C    & dyfc = *SpatialOpDatabase<DivYF2C   >::self().retrieve_operator( dim );
    SyCell     & syc  = *SpatialOpDatabase<SyCell    >::self().retrieve_operator( dim );

    /*
    EpetraExt::RowMatrixToMatrixMarketFile( "Dy.mm", dyfc.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Gy.mm", gycf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Ry.mm", rycf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Ryfc.mm", ryfc.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Sy.mm", syc.get_linalg_mat(), "", "" );
    */

    CellField     y( dim, NULL, InternalStorage );
    YSideField yint( dim, NULL, InternalStorage );
    YSideField   fy( dim, NULL, InternalStorage );
    YSideField dfdy( dim, NULL, InternalStorage );
    CellField   fy2( dim, NULL, InternalStorage );

    const double dy = spacing[1];
    const int ngylo = CellFieldTraits::GhostTraits::get<YDIR,SideMinus>();
    int ix=0;
    for( int k=klo; k<khi ; ++k ){
      for( int j=jlo; j<jhi; ++j ){
	for( int i=ilo; i<ihi; ++i ){
	  y[ix] = (j-ngylo)*dy;
	  f[ix] = std::sin(y[ix]) + 2.0;
	  ++ix;
	}
      }
    }

    rycf.apply_to_field( f, fy );
    rycf.apply_to_field( y, yint );
    ryfc.apply_to_field( fy, fy2 );
    gycf.apply_to_field( f, dfdy );
    dyfc.apply_to_field( dfdy, d2f );

    dyfc.apply_to_op( gycf, syc );

    syc.apply_to_field( f, d2f );
    
    //  EpetraExt::RowMatrixToMatrixMarketFile( "Ly.mm", syc.get_linalg_mat(), "", "" );

   // this is really a hack, but it is a way to eliminate ghost values from the output...
    tmp.reset(); tmp.add_field_contribution(y    ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "y.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(f    ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fy.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(yint ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "yy.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fy   ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "finty.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fy2  ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "finty2.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(dfdy ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "dfdy.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(d2f  ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "d2fdy2.mm", tmpField.get_linalg_vec(), "", "" );
 
  }


  if( dim[2] > 1 ){

    // get the operators
    InterpZC2F & rzcf = *SpatialOpDatabase<InterpZC2F>::self().retrieve_operator( dim );
    InterpZF2C & rzfc = *SpatialOpDatabase<InterpZF2C>::self().retrieve_operator( dim );
    GradZC2F   & gzcf = *SpatialOpDatabase<GradZC2F  >::self().retrieve_operator( dim );
    DivZF2C    & dzfc = *SpatialOpDatabase<DivZF2C   >::self().retrieve_operator( dim );
    SzCell     & szc  = *SpatialOpDatabase<SzCell    >::self().retrieve_operator( dim );

    /*
    EpetraExt::RowMatrixToMatrixMarketFile( "Dz.mm", dzfc.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Gz.mm", gzcf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Rz.mm", rzcf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Sz.mm", szc .get_linalg_mat(), "", "" );
    */

    CellField     z( dim, NULL, InternalStorage );
    ZSideField zint( dim, NULL, InternalStorage );
    ZSideField   fz( dim, NULL, InternalStorage );
    ZSideField dfdz( dim, NULL, InternalStorage );
    CellField   fz2( dim, NULL, InternalStorage );
    CellField  d2fa( dim, NULL, InternalStorage );

    const double dz = spacing[2];
    const int ngzlo = CellFieldTraits::GhostTraits::get<ZDIR,SideMinus>();
    int ix=0;
    for( int k=klo; k<khi ; ++k ){
      for( int j=jlo; j<jhi; ++j ){
	for( int i=ilo; i<ihi; ++i ){
	  z[ix] = (k-ngzlo)*dz;
	  f[ix] = std::sin(z[ix]) + 2.0;
	  ++ix;
	}
      }
    }

    rzcf.apply_to_field( f, fz );
    rzfc.apply_to_field( fz, fz2 );
    rzcf.apply_to_field( z, zint );
    gzcf.apply_to_field( f, dfdz );
    dzfc.apply_to_field( dfdz, d2fa );

    dzfc.apply_to_op( gzcf, szc );

    szc.apply_to_field( f, d2f );

    cout << endl << "Test on grad/div consistency ... ";
    ix=0;
    bool trouble = false;
    for( int k=klo; k<khi ; ++k ){
      for( int j=jlo; j<jhi; ++j ){
	for( int i=ilo; i<ihi; ++i ){
	  const double err = std::abs( d2fa[ix] - d2f[ix] );
	  if( err > 1.0e-11 ){
	    std::cout << "("<<i<<","<<j<<","<<k<<") err: " << err << std::endl;
	    trouble = true;
	  }
	  assert( err < 1.0e-10 );
	  ++ix;
	}
      }
    }
    if( trouble ) cout << "FAIL" << std::endl;
    else cout << "PASS" << std::endl;
    
    //  EpetraExt::RowMatrixToMatrixMarketFile( "Lz.mm", szc.get_linalg_mat(), "", "" );

   // this is really a hack, but it is a way to eliminate ghost values from the output...
    tmp.reset(); tmp.add_field_contribution(z    ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "z.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(f    ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(zint ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "zz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fz   ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fintz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fz2  ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "fintz2.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(dfdz ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "dfdz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(d2f  ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "d2fdz2.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(d2fa ); tmpField=tmp;  //EpetraExt::VectorToMatrixMarketFile( "d2fdz2a.mm", tmpField.get_linalg_vec(), "", "" );
 
  } // z-dir

  return true;
}

//====================================================================

bool test_linsys( const std::vector<int> & dim )
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  std::vector<double> spacing(3), area(3);
  double volume;
  setup_geom( dim, spacing, area, volume );

  LinearSystem & linSys = LinSysFactory::self().get_linsys( LinSysInfo(dim) );
  LHS & A = linSys.get_lhs();
  A.reset();


  GradXC2F   & gradX = *SpatialOpDatabase<GradXC2F  >::self().retrieve_operator( dim );
  DivXF2C    & divX = *SpatialOpDatabase<DivXF2C   >::self().retrieve_operator( dim );
  SxCell     & sxc  = *SpatialOpDatabase<SxCell    >::self().retrieve_operator( dim );
  divX.apply_to_op( gradX, sxc );
  A.add_op_contribution( sxc );

  if( dim[1]>1 ){
    GradYC2F   & gradY = *SpatialOpDatabase<GradYC2F  >::self().retrieve_operator( dim );
    DivYF2C    & divY = *SpatialOpDatabase<DivYF2C   >::self().retrieve_operator( dim );
    SyCell     & syc  = *SpatialOpDatabase<SyCell    >::self().retrieve_operator( dim );
    divY.apply_to_op( gradY, syc );
    A.add_op_contribution( syc );
  }
  if( dim[2]>1 ){
    GradZC2F   & gradZ = *SpatialOpDatabase<GradZC2F  >::self().retrieve_operator( dim );
    DivZF2C    & divZ = *SpatialOpDatabase<DivZF2C   >::self().retrieve_operator( dim );
    SzCell     & szc  = *SpatialOpDatabase<SzCell    >::self().retrieve_operator( dim );
    divZ.apply_to_op( gradZ, szc );
    A.add_op_contribution( szc );
  }

  CellField d(dim,NULL,InternalStorage);
  d = 1.0;
  A.add_field_contribution( d );

  RHS & b = linSys.get_rhs();
  b.reset( 1.0 );

  linSys.solve();

  CellFieldNoGhost tmpField( dim, NULL, InternalStorage );
  tmpField = linSys.get_soln_field();

  //EpetraExt::VectorToMatrixMarketFile( "soln.mm", tmpField.get_linalg_vec(), "", "" );

  tmpField = b;
  //EpetraExt::VectorToMatrixMarketFile( "rhs.mm", tmpField.get_linalg_vec(), "", "" );

  //EpetraExt::RowMatrixToMatrixMarketFile( "A.mm", A.epetra_mat(), "", "" );

  return true;
}

//====================================================================

void test_daixt( const vector<int>& dim )
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  cout << endl << "SpatialField operators (+, -, /, *, =) involving SpatialField objects (daixtrose)  ... ";
  bool isOkay = true;

  vector<double> d1, d2;
  const int n = dim[0]*dim[1]*dim[2];
  for( int i=0; i<n; ++i){
    d1.push_back(i);
    d2.push_back(-i+1.234);
  }

  CellFieldNoGhost f1( dim, &d1[0], ExternalStorage );
  CellFieldNoGhost f2( dim, &d2[0], ExternalStorage );

  const SpatFldPtr<CellFieldNoGhost> fp1(f1), fp2(f2);
  SpatFldPtr<CellFieldNoGhost> fp3;

  SpatFldPtr<CellFieldNoGhost> tmp = f1+f2;
  fp3 = fp1+fp2;
  for( int i=0; i<n; ++i ){
    if( std::abs((*fp3)[i] - 1.234)>1.0e-10  || (*fp3)[i] != (*tmp)[i] ){
      isOkay = false;
    }
  }

  fp3 = fp3-(f1+f2);
  for( int i=0; i<n; ++i ){
    if( std::abs((*fp3)[i])>1.0e-10 ){
      isOkay = false;
    }
  }


  fp3 = (f1+f2)*fp1;
  for( int i=0; i<n; ++i ){
    const double ans = f1[i]*1.234;
    const double abserr = std::abs((*fp3)[i] - ans);
    const double relerr = abserr/std::abs(ans);
    if( abserr>1.0e-10 && relerr>1.0e-8 ){
      isOkay = false;
    }
  }

  fp3 = f1*f2+f1/f2*fp3+f2/f1+f2*f1*f2;

  if( isOkay )  cout << "PASS" << endl;
  else          cout << "FAIL!" << endl;
}

//====================================================================

template<typename FieldT>
bool test_op_rhs_assign( const vector<int>& dim )
{
  RHS rhs(dim);
  FieldT field( dim, NULL, SpatialOps::InternalStorage );

  bool err = false;

  int ix=0;
  for( typename FieldT::iterator ifld=field.begin(); ifld!=field.end(); ++ifld, ++ix ){
    *ifld = ix;
  }
  rhs.reset();
  rhs.add_field_contribution( field );
  field = rhs;
  ix=0;
  for( typename FieldT::iterator ifld=field.begin(); ifld!=field.end(); ++ifld, ++ix ){
    if( *ifld != ix )  err = true;
  }
  return err;
}

//====================================================================

void test_rhs_ops()
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  vector<int> dim(3,1);
  dim[0] = 10;
  dim[1] = 7;
  dim[2] = 13;

  cout << endl
       << "SpatialField operators involving RHS objects: "
       << endl;

  cout << "  CellField  ... ";
  if( test_op_rhs_assign<CellField>( dim ) )    cout << "FAIL!" << endl;
  else cout << "PASS" << endl;

  cout << "  XSideField ... ";
  if( test_op_rhs_assign<XSideField>( dim ) )   cout << "FAIL!" << endl;
  else cout << "PASS" << endl;

  cout << "  YSideField ... ";
  if( test_op_rhs_assign<XSideField>( dim ) )   cout << "FAIL!" << endl;
  else cout << "PASS" << endl;

  cout << "  ZSideField ... ";
  if( test_op_rhs_assign<XSideField>( dim ) )   cout << "FAIL!" << endl;
  else cout << "PASS" << endl;

  return;
}

//====================================================================

#include <FVStaggeredBCTools.h>

template<typename OpT,typename FieldT,typename DestFieldT,typename Side>
bool test_bc_helper( const vector<int>&dim,
		     const int ii,
		     const int jj,
		     const int kk,
		     const double bcVal )
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  const OpT& op = *SpatialOpDatabase<OpT>::self().retrieve_operator(dim);

  FieldT      f( dim, NULL );
  DestFieldT df( dim, NULL );

  int icnt=0;
  for( typename FieldT::iterator ifld=f.begin(); ifld!=f.end(); ++ifld,++icnt )
    *ifld = icnt;

  // assign the BC.
  assign_bc_point<OpT,FieldT,Side>( op, ii, jj, kk, bcVal, f );

  // calculate the dest field
  op.apply_to_field( f, df );

  // verify that the BC was set properly - this is a bit of a hack.
  int ix=-1;
  ijk_interior_to_flatix<typename DestFieldT::Ghost>( ii,jj,kk,dim,ix );
  BCToolsLocal::RowIndexShifter<typename DestFieldT::Ghost,typename OpT::DirType,Side>::apply(dim,ix);

  const double abserr = abs(df[ix]-bcVal);
  const double relerr = abserr/abs(bcVal);
//   if( abserr>1.0e-10 && relerr>1.0e-8 ){
//     cout << "(" << ii << "," << jj << "," << kk << ") : " << df[ix] << " ; " << bcVal << endl;
//   }
  return ( abserr<1.0e-10 && relerr<1.0e-8 );
}
//--------------------------------------------------------------------
void test_bc()
{
  using namespace SpatialOps;
  using namespace FVStaggeredUniform;

  vector<int> dim;
  dim.push_back(4);
  dim.push_back(4);
  dim.push_back(4);

  // be sure we have operators built for this mesh.
  build_ops( dim );

  cout << endl << "Testing BC setting stuff:" << endl;

  bool isOkay = true;;

  // X BCs - Left side
  cout << "  X Dir, (-) side ... ";
  int i=0;
  for( int j=0; j<dim[1]; ++j ){
    for( int k=0; k<dim[2]; ++k ){
      const bool result1 = test_bc_helper<GradXC2F,  CellField,XSideField,SideMinus>( dim, i,j,k, 1.2345 );
      const bool result2 = test_bc_helper<InterpXC2F,CellField,XSideField,SideMinus>( dim, i,j,k, 123.456 );
      if( !result1 || !result2 ) isOkay=false;
    }
  }
  if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;

  // X BCs - Right side
  cout << "  X Dir, (+) side ... ";
  isOkay = true;
  i=dim[0]-1;
  for( int j=0; j<dim[1]; ++j ){
    for( int k=0; k<dim[2]; ++k ){
      const bool result1 = test_bc_helper<GradXC2F,  CellField,XSideField,SidePlus >( dim, i,j,k, 5.4321 );
      const bool result2 = test_bc_helper<InterpXC2F,CellField,XSideField,SidePlus >( dim, i,j,k, 123.456 );
      if( !result1 || !result2 ) isOkay=false;
    }
  }
  if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;


  // Y BCs - Left side
  cout << "  Y Dir, (-) side ... ";
  int j=0;
  for( int i=0; i<dim[0]; ++i ){
    for( int k=0; k<dim[2]; ++k ){
      const bool result1 = test_bc_helper<GradYC2F,  CellField,YSideField,SideMinus>( dim, i,j,k, 1.2345 );
      const bool result2 = test_bc_helper<InterpYC2F,CellField,YSideField,SideMinus>( dim, i,j,k, 123.456 );
      if( !result1 || !result2 ) isOkay=false;
    }
  }
  if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;

  // Y BCs - Right side
  cout << "  Y Dir, (+) side ... ";
  isOkay = true;
  j=dim[1]-1;
  for( int i=0; i<dim[0]; ++i ){
    for( int k=0; k<dim[2]; ++k ){
      const bool result1 = test_bc_helper<GradYC2F,  CellField,YSideField,SidePlus >( dim, i,j,k, 5.4321 );
      const bool result2 = test_bc_helper<InterpYC2F,CellField,YSideField,SidePlus >( dim, i,j,k, 123.456 );
      if( !result1 || !result2 ) isOkay=false;
    }
  }
  if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;



  // Z BCs - Left side
  cout << "  Z Dir, (-) side ... ";
  int k=0;
  for( int i=0; i<dim[0]; ++i ){
    for( int j=0; j<dim[1]; ++j ){
      const bool result1 = test_bc_helper<GradZC2F,  CellField,ZSideField,SideMinus>( dim, i,j,k, 1.2345 );
      const bool result2 = test_bc_helper<InterpZC2F,CellField,ZSideField,SideMinus>( dim, i,j,k, 123.456 );
      if( !result1 || !result2 ) isOkay=false;
    }
  }
  if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;

  // Z BCs - Right side
  cout << "  Z Dir, (+) side ... ";
  isOkay = true;
  k=dim[2]-1;
  for( int i=0; i<dim[0]; ++i ){
    for( int j=0; j<dim[1]; ++j ){
      const bool result1 = test_bc_helper<GradZC2F,  CellField,ZSideField,SidePlus >( dim, i,j,k, 5.4321 );
      const bool result2 = test_bc_helper<InterpZC2F,CellField,ZSideField,SidePlus >( dim, i,j,k, 123.456 );
      if( !result1 || !result2 ) isOkay=false;
    }
  }
  if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;

  cout << endl;
}

//====================================================================

int main()
{
  vector<int> dim(3,1);
  dim[0]=  4;
  dim[1]=  23;
  dim[2]=  18;

  try{

    build_ops( dim );

    test_daixt( dim );

    test( dim );
    test_linsys( dim );

    test_rhs_ops();

    test_bc();
  }
  catch( std::runtime_error& e ){
    cout << endl << e.what() << endl;
  }
  return 0;
}
