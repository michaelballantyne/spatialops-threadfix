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

  // build the assemblers
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

  SxCell::Assembler  sxcasmbl( dim );
  SyCell::Assembler  sycasmbl( dim );
  SzCell::Assembler  szcasmbl( dim );

  SxCellSide::Assembler  sxcfasmbl( dim );
  SyCellSide::Assembler  sycfasmbl( dim );
  SzCellSide::Assembler  szcfasmbl( dim );


  // build the operators
  InterpXC2F *rxcf = new InterpXC2F( rxcfasmbl );
  InterpYC2F *rycf = new InterpYC2F( rycfasmbl );
  InterpZC2F *rzcf = new InterpZC2F( rzcfasmbl );

  InterpXF2C *rxfc = new InterpXF2C( rxfcasmbl );
  InterpYF2C *ryfc = new InterpYF2C( ryfcasmbl );
  InterpZF2C *rzfc = new InterpZF2C( rzfcasmbl );

  GradXC2F   *gxcf = new GradXC2F( gxcfasmbl );
  GradYC2F   *gycf = new GradYC2F( gycfasmbl );
  GradZC2F   *gzcf = new GradZC2F( gzcfasmbl );

  GradXF2C   *gxfc = new GradXF2C( gxfcasmbl );
  GradYF2C   *gyfc = new GradYF2C( gyfcasmbl );
  GradZF2C   *gzfc = new GradZF2C( gzfcasmbl );

  DivXF2C    *dxfc = new DivXF2C( dxfcasmbl );
  DivYF2C    *dyfc = new DivYF2C( dyfcasmbl );
  DivZF2C    *dzfc = new DivZF2C( dzfcasmbl );

  SxCellSide *sxcf = new SxCellSide( sxcfasmbl );
  SyCellSide *sycf = new SyCellSide( sycfasmbl );
  SzCellSide *szcf = new SzCellSide( szcfasmbl );

  SxCell     *sxc  = new SxCell( sxcasmbl  );
  SyCell     *syc  = new SyCell( sycasmbl  );
  SzCell     *szc  = new SzCell( szcasmbl  );

  SpatialOpDatabase<InterpXC2F>::self().register_new_operator( rxcf );
  SpatialOpDatabase<InterpYC2F>::self().register_new_operator( rycf );
  SpatialOpDatabase<InterpZC2F>::self().register_new_operator( rzcf );

  SpatialOpDatabase<InterpXF2C>::self().register_new_operator( rxfc );
  SpatialOpDatabase<InterpYF2C>::self().register_new_operator( ryfc );
  SpatialOpDatabase<InterpZF2C>::self().register_new_operator( rzfc );

  SpatialOpDatabase<GradXC2F  >::self().register_new_operator( gxcf );
  SpatialOpDatabase<GradYC2F  >::self().register_new_operator( gycf );
  SpatialOpDatabase<GradZC2F  >::self().register_new_operator( gzcf );

  SpatialOpDatabase<GradXF2C  >::self().register_new_operator( gxfc );
  SpatialOpDatabase<GradYF2C  >::self().register_new_operator( gyfc );
  SpatialOpDatabase<GradZF2C  >::self().register_new_operator( gzfc );

  SpatialOpDatabase<DivXF2C   >::self().register_new_operator( dxfc );
  SpatialOpDatabase<DivYF2C   >::self().register_new_operator( dyfc );
  SpatialOpDatabase<DivZF2C   >::self().register_new_operator( dzfc );

  SpatialOpDatabase<SxCellSide>::self().register_new_operator( sxcf );
  SpatialOpDatabase<SyCellSide>::self().register_new_operator( sycf );
  SpatialOpDatabase<SzCellSide>::self().register_new_operator( szcf );

  SpatialOpDatabase<SxCell    >::self().register_new_operator( sxc );
  SpatialOpDatabase<SyCell    >::self().register_new_operator( syc );
  SpatialOpDatabase<SzCell    >::self().register_new_operator( szc );
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

  EpetraExt::RowMatrixToMatrixMarketFile( "Dx.mm", dxfc.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Gx.mm", gxcf.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Rx.mm", rxcf.get_linalg_mat(), "", "" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Sx.mm", sxc .get_linalg_mat(), "", "" );


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
  EpetraExt::RowMatrixToMatrixMarketFile( "Lx.mm", sxc.get_linalg_mat(), "", "" );

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
  tmp.reset(); tmp.add_field_contribution(x    ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "x.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(f    ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(xint ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "xx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(fint ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fintx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(gradF); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "dfdx.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(d2f  ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "d2fdx2.mm", tmpField.get_linalg_vec(), "", "" );
  tmp.reset(); tmp.add_field_contribution(fint2); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fintx2.mm", tmpField.get_linalg_vec(), "", "" );


  if( dim[1] > 1 ){

    // get the operators
    InterpYC2F & rycf = *SpatialOpDatabase<InterpYC2F>::self().retrieve_operator( dim );
    InterpYF2C & ryfc = *SpatialOpDatabase<InterpYF2C>::self().retrieve_operator( dim );
    GradYC2F   & gycf = *SpatialOpDatabase<GradYC2F  >::self().retrieve_operator( dim );
    DivYF2C    & dyfc = *SpatialOpDatabase<DivYF2C   >::self().retrieve_operator( dim );
    SyCell     & syc  = *SpatialOpDatabase<SyCell    >::self().retrieve_operator( dim );

    EpetraExt::RowMatrixToMatrixMarketFile( "Dy.mm", dyfc.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Gy.mm", gycf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Ry.mm", rycf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Sy.mm", syc.get_linalg_mat(), "", "" );


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
    
    EpetraExt::RowMatrixToMatrixMarketFile( "Ly.mm", syc.get_linalg_mat(), "", "" );

   // this is really a hack, but it is a way to eliminate ghost values from the output...
    tmp.reset(); tmp.add_field_contribution(y    ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "y.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(f    ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fy.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(yint ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "yy.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fy   ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "finty.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fy2  ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "finty2.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(dfdy ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "dfdy.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(d2f  ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "d2fdy2.mm", tmpField.get_linalg_vec(), "", "" );
 
  }


  if( dim[2] > 1 ){

    // get the operators
    InterpZC2F & rzcf = *SpatialOpDatabase<InterpZC2F>::self().retrieve_operator( dim );
    InterpZF2C & rzfc = *SpatialOpDatabase<InterpZF2C>::self().retrieve_operator( dim );
    GradZC2F   & gzcf = *SpatialOpDatabase<GradZC2F  >::self().retrieve_operator( dim );
    DivZF2C    & dzfc = *SpatialOpDatabase<DivZF2C   >::self().retrieve_operator( dim );
    SzCell     & szc  = *SpatialOpDatabase<SzCell    >::self().retrieve_operator( dim );

    EpetraExt::RowMatrixToMatrixMarketFile( "Dz.mm", dzfc.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Gz.mm", gzcf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Rz.mm", rzcf.get_linalg_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "Sz.mm", szc .get_linalg_mat(), "", "" );


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

    ix=0;
    for( int k=klo; k<khi ; ++k ){
      for( int j=jlo; j<jhi; ++j ){
	for( int i=ilo; i<ihi; ++i ){
	  const double err = std::abs( d2fa[ix] - d2f[ix] );
	  if( err > 1.0e-11 ) std::cout << "("<<i<<","<<j<<","<<k<<") err: " << err << std::endl;
	  assert( err < 1.0e-10 );
	  ++ix;
	}
      }
    }
    
    EpetraExt::RowMatrixToMatrixMarketFile( "Lz.mm", szc.get_linalg_mat(), "", "" );

   // this is really a hack, but it is a way to eliminate ghost values from the output...
    tmp.reset(); tmp.add_field_contribution(z    ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "z.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(f    ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(zint ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "zz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fz   ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fintz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(fz2  ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "fintz2.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(dfdz ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "dfdz.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(d2f  ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "d2fdz2.mm", tmpField.get_linalg_vec(), "", "" );
    tmp.reset(); tmp.add_field_contribution(d2fa ); tmpField=tmp;  EpetraExt::VectorToMatrixMarketFile( "d2fdz2a.mm", tmpField.get_linalg_vec(), "", "" );
 
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

  // get the operators
  GradXC2F   & gradX = *SpatialOpDatabase<GradXC2F  >::self().retrieve_operator( dim );
  GradYC2F   & gradY = *SpatialOpDatabase<GradYC2F  >::self().retrieve_operator( dim );
  GradZC2F   & gradZ = *SpatialOpDatabase<GradZC2F  >::self().retrieve_operator( dim );

  DivXF2C    & divX = *SpatialOpDatabase<DivXF2C   >::self().retrieve_operator( dim );
  DivYF2C    & divY = *SpatialOpDatabase<DivYF2C   >::self().retrieve_operator( dim );
  DivZF2C    & divZ = *SpatialOpDatabase<DivZF2C   >::self().retrieve_operator( dim );

  SxCell     & sxc  = *SpatialOpDatabase<SxCell    >::self().retrieve_operator( dim );
  SyCell     & syc  = *SpatialOpDatabase<SyCell    >::self().retrieve_operator( dim );
  SzCell     & szc  = *SpatialOpDatabase<SzCell    >::self().retrieve_operator( dim );


  divX.apply_to_op( gradX, sxc );
  divY.apply_to_op( gradY, syc );
  divZ.apply_to_op( gradZ, szc );

  LinearSystem & linSys = LinSysFactory::self().get_linsys( LinSysInfo(dim) );

  LHS & A = linSys.get_lhs();
  A.reset();
  A.add_op_contribution( sxc );
  A.add_op_contribution( syc );
  A.add_op_contribution( szc );
  RHS & b = linSys.get_rhs();
  b.reset( 1.0 );

  linSys.solve();

  CellFieldNoGhost tmpField( dim, NULL, InternalStorage );
  tmpField = linSys.get_soln_field();

  EpetraExt::VectorToMatrixMarketFile( "soln.mm", tmpField.get_linalg_vec(), "", "" );

  tmpField = b;
  EpetraExt::VectorToMatrixMarketFile( "rhs.mm", tmpField.get_linalg_vec(), "", "" );

  EpetraExt::RowMatrixToMatrixMarketFile( "A.mm", A.epetra_mat(), "", "" );

  return true;
}

//====================================================================

int main()
{
  vector<int> dim(3,1);
  dim[0]=10;
  dim[1]=13;
  dim[2]=31;

  build_ops( dim );

  test( dim );
  test_linsys( dim );
  return 0;
}
