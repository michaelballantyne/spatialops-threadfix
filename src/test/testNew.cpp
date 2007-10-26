#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
using namespace std;


#include <EpetraExt_VectorOut.h>
#include <EpetraExt_RowMatrixOut.h>

#include <FVStaggered.h>
#include <FVStaggeredBCTools.h>
#include <LinearSystem.h>

#include <Grid.h>
#include <Functions.h>

using namespace SpatialOps;
using namespace FVStaggered;

static int ssxid, ssyid, sszid;
static int xsxid, xsyid, xszid;
static int ysxid, ysyid, yszid;
static int zsxid, zsyid, zszid;

//--------------------------------------------------------------------

void build_ops( const std::vector<int>& dim,
		const std::vector<double>& spacing )
{
  using SpatialOps::SpatialOpDatabase;

  std::vector<double> area(3,1);
  area[0] = spacing[1]*spacing[2];
  area[1] = spacing[0]*spacing[2];
  area[2] = spacing[0]*spacing[1];
  const double vol = spacing[0]*spacing[1]*spacing[2];


  // Divergence operators 
  {
    DivSSurfXSVol::Assembler Dssxsva( dim, area[0], vol );
    DivSSurfYSVol::Assembler Dssysva( dim, area[1], vol );
    DivSSurfZSVol::Assembler Dsszsva( dim, area[2], vol );

    DivXSurfXXVol::Assembler Dxsxxva( dim, area[0], vol );
    DivXSurfYXVol::Assembler Dxsyxva( dim, area[1], vol );
    DivXSurfZXVol::Assembler Dxszxva( dim, area[2], vol );

    DivYSurfXYVol::Assembler Dysxyva( dim, area[0], vol );
    DivYSurfYYVol::Assembler Dysyyva( dim, area[1], vol );
    DivYSurfZYVol::Assembler Dyszyva( dim, area[2], vol );

    DivZSurfXZVol::Assembler Dzsxzva( dim, area[0], vol );
    DivZSurfYZVol::Assembler Dzsyzva( dim, area[1], vol );
    DivZSurfZZVol::Assembler Dzszzva( dim, area[2], vol );

    SpatialOpDatabase<DivSSurfXSVol>::self().register_new_operator( new DivSSurfXSVol(Dssxsva) );
    SpatialOpDatabase<DivSSurfYSVol>::self().register_new_operator( new DivSSurfYSVol(Dssysva) );
    SpatialOpDatabase<DivSSurfZSVol>::self().register_new_operator( new DivSSurfZSVol(Dsszsva) );

    SpatialOpDatabase<DivXSurfXXVol>::self().register_new_operator( new DivXSurfXXVol(Dxsxxva) );
    SpatialOpDatabase<DivXSurfYXVol>::self().register_new_operator( new DivXSurfYXVol(Dxsyxva) );
    SpatialOpDatabase<DivXSurfZXVol>::self().register_new_operator( new DivXSurfZXVol(Dxszxva) );

    SpatialOpDatabase<DivYSurfXYVol>::self().register_new_operator( new DivYSurfXYVol(Dysxyva) );
    SpatialOpDatabase<DivYSurfYYVol>::self().register_new_operator( new DivYSurfYYVol(Dysyyva) );
    SpatialOpDatabase<DivYSurfZYVol>::self().register_new_operator( new DivYSurfZYVol(Dyszyva) );

    SpatialOpDatabase<DivZSurfXZVol>::self().register_new_operator( new DivZSurfXZVol(Dzsxzva) );
    SpatialOpDatabase<DivZSurfYZVol>::self().register_new_operator( new DivZSurfYZVol(Dzsyzva) );
    SpatialOpDatabase<DivZSurfZZVol>::self().register_new_operator( new DivZSurfZZVol(Dzszzva) );
  }

  // gradient operators - diffusive fluxes
  {
    GradSVolSSurfX::Assembler Gsvssxa( spacing[0], dim );
    GradSVolSSurfY::Assembler Gsvssya( spacing[1], dim );
    GradSVolSSurfZ::Assembler Gsvssza( spacing[2], dim );

    GradXVolXSurfX::Assembler Gxvxsxa( spacing[0], dim );
    GradXVolXSurfY::Assembler Gxvxsya( spacing[1], dim );
    GradXVolXSurfZ::Assembler Gxvxsza( spacing[2], dim );

    GradYVolYSurfX::Assembler Gyvysxa( spacing[0], dim );
    GradYVolYSurfY::Assembler Gyvysya( spacing[1], dim );
    GradYVolYSurfZ::Assembler Gyvysza( spacing[2], dim );

    GradZVolZSurfX::Assembler Gzvzsxa( spacing[0], dim );
    GradZVolZSurfY::Assembler Gzvzsya( spacing[1], dim );
    GradZVolZSurfZ::Assembler Gzvzsza( spacing[2], dim );

    SpatialOpDatabase<GradSVolSSurfX>::self().register_new_operator( new GradSVolSSurfX(Gsvssxa) );
    SpatialOpDatabase<GradSVolSSurfY>::self().register_new_operator( new GradSVolSSurfY(Gsvssya) );
    SpatialOpDatabase<GradSVolSSurfZ>::self().register_new_operator( new GradSVolSSurfZ(Gsvssza) );

    SpatialOpDatabase<GradXVolXSurfX>::self().register_new_operator( new GradXVolXSurfX(Gxvxsxa) );
    SpatialOpDatabase<GradXVolXSurfY>::self().register_new_operator( new GradXVolXSurfY(Gxvxsya) );
    SpatialOpDatabase<GradXVolXSurfZ>::self().register_new_operator( new GradXVolXSurfZ(Gxvxsza) );

    SpatialOpDatabase<GradYVolYSurfX>::self().register_new_operator( new GradYVolYSurfX(Gyvysxa) );
    SpatialOpDatabase<GradYVolYSurfY>::self().register_new_operator( new GradYVolYSurfY(Gyvysya) );
    SpatialOpDatabase<GradYVolYSurfZ>::self().register_new_operator( new GradYVolYSurfZ(Gyvysza) );

    SpatialOpDatabase<GradZVolZSurfX>::self().register_new_operator( new GradZVolZSurfX(Gzvzsxa) );
    SpatialOpDatabase<GradZVolZSurfY>::self().register_new_operator( new GradZVolZSurfY(Gzvzsya) );
    SpatialOpDatabase<GradZVolZSurfZ>::self().register_new_operator( new GradZVolZSurfZ(Gzvzsza) );

  }

  // interpolant scalar volume to scalar surface (diffusivities)
  {
//     InterpSVolSSurf::Assembler Rsvssa( dim );
//     SpatialOpDatabase<InterpSVolSSurf>::self().register_new_operator( new InterpSVolSSurf(Rsvssa) );

    {
      InterpSVolSSurfX::Assembler Rsvssxa( dim );
      InterpSVolSSurfY::Assembler Rsvssya( dim );
      InterpSVolSSurfZ::Assembler Rsvssza( dim );

      SpatialOpDatabase<InterpSVolSSurfX>::self().register_new_operator( new InterpSVolSSurfX(Rsvssxa) );
      SpatialOpDatabase<InterpSVolSSurfY>::self().register_new_operator( new InterpSVolSSurfY(Rsvssya) );
      SpatialOpDatabase<InterpSVolSSurfZ>::self().register_new_operator( new InterpSVolSSurfZ(Rsvssza) );
    }
  }

  // interpolants - scalar volume to staggered surfaces (viscosity, dilatation)
  {
//     InterpSVolXSurfX::Assembler Rsvxsxa( dim );
//     InterpSVolXSurfY::Assembler Rsvxsya( dim );
//     InterpSVolXSurfZ::Assembler Rsvxsza( dim );
//     SpatialOpDatabase<InterpSVolXSurfX>::self().register_new_operator( new InterpSVolXSurfX(Rsvxsxa) );
//     SpatialOpDatabase<InterpSVolXSurfX>::self().register_new_operator( new InterpSVolXSurfX(Rsvxsxa) );
//     SpatialOpDatabase<InterpSVolXSurfX>::self().register_new_operator( new InterpSVolXSurfX(Rsvxsxa) );

//     InterpSVolYSurfX::Assembler Rsvysxa( dim );
//     InterpSVolYSurfY::Assembler Rsvysya( dim );
//     InterpSVolYSurfZ::Assembler Rsvysza( dim );
//     SpatialOpDatabase<InterpSVolYSurfX>::self().register_new_operator( new InterpSVolYSurfX(Rsvysxa) );
//     SpatialOpDatabase<InterpSVolYSurfY>::self().register_new_operator( new InterpSVolYSurfY(Rsvysya) );
//     SpatialOpDatabase<InterpSVolYSurfZ>::self().register_new_operator( new InterpSVolYSurfZ(Rsvysza) );

//     InterpSVolZSurfX::Assembler Rsvzsxa( dim );
//     InterpSVolZSurfY::Assembler Rsvzsya( dim );
//     InterpSVolZSurfZ::Assembler Rsvzsza( dim );
//     SpatialOpDatabase<InterpSVolZSurfX>::self().register_new_operator( new InterpSVolZSurfX(Rsvzsxa) );
//     SpatialOpDatabase<InterpSVolZSurfY>::self().register_new_operator( new InterpSVolZSurfY(Rsvzsya) );
//     SpatialOpDatabase<InterpSVolZSurfZ>::self().register_new_operator( new InterpSVolZSurfZ(Rsvzsza) );
  }

  // interpolants - scalar volume to staggered volume (density)
  {
    InterpSVolXVol::Assembler Rsvxva( dim );
    InterpSVolYVol::Assembler Rsvyva( dim );
    InterpSVolZVol::Assembler Rsvzva( dim );

    SpatialOpDatabase<InterpSVolXVol>::self().register_new_operator( new InterpSVolXVol(Rsvxva) );
    SpatialOpDatabase<InterpSVolYVol>::self().register_new_operator( new InterpSVolYVol(Rsvyva) );
    SpatialOpDatabase<InterpSVolZVol>::self().register_new_operator( new InterpSVolZVol(Rsvzva) );
  }

  // interpolants - staggered volume to staggered surfaces (advecting velocities)
  {
    InterpXVolYSurfX::Assembler Rxvysxa( dim );
    InterpXVolZSurfX::Assembler Rxvzsxa( dim );

    InterpYVolXSurfY::Assembler Ryvxsya( dim );
    InterpYVolZSurfY::Assembler Ryvzsya( dim );

    InterpZVolXSurfZ::Assembler Rzvxsza( dim );
    InterpZVolYSurfZ::Assembler Rzvysza( dim );

    SpatialOpDatabase<InterpXVolYSurfX>::self().register_new_operator( new InterpXVolYSurfX(Rxvysxa) );
    SpatialOpDatabase<InterpXVolZSurfX>::self().register_new_operator( new InterpXVolZSurfX(Rxvzsxa) );

    SpatialOpDatabase<InterpYVolXSurfY>::self().register_new_operator( new InterpYVolXSurfY(Ryvxsya) );
    SpatialOpDatabase<InterpYVolZSurfY>::self().register_new_operator( new InterpYVolZSurfY(Ryvzsya) );

    SpatialOpDatabase<InterpZVolXSurfZ>::self().register_new_operator( new InterpZVolXSurfZ(Rzvxsza) );
    SpatialOpDatabase<InterpZVolYSurfZ>::self().register_new_operator( new InterpZVolYSurfZ(Rzvysza) );

    // interpolants - volume to surface for staggered cells.
    {
//       InterpXVolXSurf::Assembler Rxvxsa( dim );
//       SpatialOpDatabase<InterpXVolXSurf>::self().register_new_operator( new InterpXVolXSurf(Rxvxsa) );

      InterpXVolXSurfX::Assembler Rxvxsxa( dim );
      InterpXVolXSurfY::Assembler Rxvxsya( dim );
      InterpXVolXSurfZ::Assembler Rxvxsza( dim );

      SpatialOpDatabase<InterpXVolXSurfX>::self().register_new_operator( new InterpXVolXSurfX(Rxvxsxa) );
      SpatialOpDatabase<InterpXVolXSurfY>::self().register_new_operator( new InterpXVolXSurfY(Rxvxsya) );
      SpatialOpDatabase<InterpXVolXSurfZ>::self().register_new_operator( new InterpXVolXSurfZ(Rxvxsza) );



//       InterpYVolYSurf::Assembler Ryvysa( dim );
//       SpatialOpDatabase<InterpYVolYSurf>::self().register_new_operator( new InterpYVolYSurf(Ryvysa) );

      InterpYVolYSurfX::Assembler Ryvysxa( dim );
      InterpYVolYSurfY::Assembler Ryvysya( dim );
      InterpYVolYSurfZ::Assembler Ryvysza( dim );

      SpatialOpDatabase<InterpYVolYSurfX>::self().register_new_operator( new InterpYVolYSurfX(Ryvysxa) );
      SpatialOpDatabase<InterpYVolYSurfY>::self().register_new_operator( new InterpYVolYSurfY(Ryvysya) );
      SpatialOpDatabase<InterpYVolYSurfZ>::self().register_new_operator( new InterpYVolYSurfZ(Ryvysza) );



//       InterpZVolZSurf::Assembler Rzvzsa( dim );
//       SpatialOpDatabase<InterpZVolZSurf>::self().register_new_operator( new InterpZVolZSurf(Rzvzsa) );

      InterpZVolZSurfX::Assembler Rzvzsxa( dim );
      InterpZVolZSurfY::Assembler Rzvzsya( dim );
      InterpZVolZSurfZ::Assembler Rzvzsza( dim );

      SpatialOpDatabase<InterpZVolZSurfX>::self().register_new_operator( new InterpZVolZSurfX(Rzvzsxa) );
      SpatialOpDatabase<InterpZVolZSurfY>::self().register_new_operator( new InterpZVolZSurfY(Rzvzsya) );
      SpatialOpDatabase<InterpZVolZSurfZ>::self().register_new_operator( new InterpZVolZSurfZ(Rzvzsza) );
    }
  }

  // interpolants - staggard volume to scalar surface  (advecting velocities)
  {
    InterpXVolSSurfX::Assembler Rxvssxa( dim );
    InterpYVolSSurfY::Assembler Ryvssya( dim );
    InterpZVolSSurfZ::Assembler Rzvssza( dim );

    SpatialOpDatabase<InterpXVolSSurfX>::self().register_new_operator( new InterpXVolSSurfX(Rxvssxa) );
    SpatialOpDatabase<InterpYVolSSurfY>::self().register_new_operator( new InterpYVolSSurfY(Ryvssya) );
    SpatialOpDatabase<InterpZVolSSurfZ>::self().register_new_operator( new InterpZVolSSurfZ(Rzvssza) );
  }


  // scratch operators
  {
    ScratchSVol::Assembler Ssvxa( dim, XDIR::value );
    ScratchSVol::Assembler Ssvya( dim, YDIR::value );
    ScratchSVol::Assembler Ssvza( dim, ZDIR::value );

    ScratchXVol::Assembler Sxvxa( dim, XDIR::value );
    ScratchXVol::Assembler Sxvya( dim, YDIR::value );
    ScratchXVol::Assembler Sxvza( dim, ZDIR::value );

    ScratchYVol::Assembler Syvxa( dim, XDIR::value );
    ScratchYVol::Assembler Syvya( dim, YDIR::value );
    ScratchYVol::Assembler Syvza( dim, ZDIR::value );

    ScratchZVol::Assembler Szvxa( dim, XDIR::value );
    ScratchZVol::Assembler Szvya( dim, YDIR::value );
    ScratchZVol::Assembler Szvza( dim, ZDIR::value );
    
    ssxid = SpatialOpDatabase<ScratchSVol>::self().register_new_operator( new ScratchSVol(Ssvxa) );
    ssyid = SpatialOpDatabase<ScratchSVol>::self().register_new_operator( new ScratchSVol(Ssvya) );
    sszid = SpatialOpDatabase<ScratchSVol>::self().register_new_operator( new ScratchSVol(Ssvza) );

    xsxid = SpatialOpDatabase<ScratchXVol>::self().register_new_operator( new ScratchXVol(Sxvxa) );
    xsyid = SpatialOpDatabase<ScratchXVol>::self().register_new_operator( new ScratchXVol(Sxvya) );
    xszid = SpatialOpDatabase<ScratchXVol>::self().register_new_operator( new ScratchXVol(Sxvza) );

    ysxid = SpatialOpDatabase<ScratchYVol>::self().register_new_operator( new ScratchYVol(Syvxa) );
    ysyid = SpatialOpDatabase<ScratchYVol>::self().register_new_operator( new ScratchYVol(Syvya) );
    yszid = SpatialOpDatabase<ScratchYVol>::self().register_new_operator( new ScratchYVol(Syvza) );

    zsxid = SpatialOpDatabase<ScratchZVol>::self().register_new_operator( new ScratchZVol(Szvxa) );
    zsyid = SpatialOpDatabase<ScratchZVol>::self().register_new_operator( new ScratchZVol(Szvya) );
    zszid = SpatialOpDatabase<ScratchZVol>::self().register_new_operator( new ScratchZVol(Szvza) );


    // scratch grad ops - need to test these.
 
    ScratchSVolSSurfX::Assembler Ssvssxa( spacing[0], dim );
    ScratchSVolSSurfY::Assembler Ssvssya( spacing[1], dim );
    ScratchSVolSSurfZ::Assembler Ssvssza( spacing[2], dim );

    ScratchXVolXSurfX::Assembler Sxvxsxa( spacing[0], dim );
    ScratchXVolXSurfY::Assembler Sxvxsya( spacing[1], dim );
    ScratchXVolXSurfZ::Assembler Sxvxsza( spacing[2], dim );

    ScratchYVolYSurfX::Assembler Syvysxa( spacing[0], dim );
    ScratchYVolYSurfY::Assembler Syvysya( spacing[1], dim );
    ScratchYVolYSurfZ::Assembler Syvysza( spacing[2], dim );

    ScratchZVolZSurfX::Assembler Szvzsxa( spacing[0], dim );
    ScratchZVolZSurfY::Assembler Szvzsya( spacing[1], dim );
    ScratchZVolZSurfZ::Assembler Szvzsza( spacing[2], dim );

    SpatialOpDatabase<ScratchSVolSSurfX>::self().register_new_operator( new ScratchSVolSSurfX(Ssvssxa) );
    SpatialOpDatabase<ScratchSVolSSurfY>::self().register_new_operator( new ScratchSVolSSurfY(Ssvssya) );
    SpatialOpDatabase<ScratchSVolSSurfZ>::self().register_new_operator( new ScratchSVolSSurfZ(Ssvssza) );

    SpatialOpDatabase<ScratchXVolXSurfX>::self().register_new_operator( new ScratchXVolXSurfX(Sxvxsxa) );
    SpatialOpDatabase<ScratchXVolXSurfY>::self().register_new_operator( new ScratchXVolXSurfY(Sxvxsya) );
    SpatialOpDatabase<ScratchXVolXSurfZ>::self().register_new_operator( new ScratchXVolXSurfZ(Sxvxsza) );

    SpatialOpDatabase<ScratchYVolYSurfX>::self().register_new_operator( new ScratchYVolYSurfX(Syvysxa) );
    SpatialOpDatabase<ScratchYVolYSurfY>::self().register_new_operator( new ScratchYVolYSurfY(Syvysya) );
    SpatialOpDatabase<ScratchYVolYSurfZ>::self().register_new_operator( new ScratchYVolYSurfZ(Syvysza) );

    SpatialOpDatabase<ScratchZVolZSurfX>::self().register_new_operator( new ScratchZVolZSurfX(Szvzsxa) );
    SpatialOpDatabase<ScratchZVolZSurfY>::self().register_new_operator( new ScratchZVolZSurfY(Szvzsya) );
    SpatialOpDatabase<ScratchZVolZSurfZ>::self().register_new_operator( new ScratchZVolZSurfZ(Szvzsza) );

    // bc operators
    ScratchSSurfX::Assembler Sssxa( dim, XDIR::value );
    ScratchSSurfY::Assembler Sssya( dim, YDIR::value );
    ScratchSSurfZ::Assembler Sssza( dim, ZDIR::value );

    SpatialOpDatabase<ScratchSSurfX>::self().register_new_operator( new ScratchSSurfX(Sssxa) );
    SpatialOpDatabase<ScratchSSurfY>::self().register_new_operator( new ScratchSSurfY(Sssya) );
    SpatialOpDatabase<ScratchSSurfZ>::self().register_new_operator( new ScratchSSurfZ(Sssza) );
  }
}

//--------------------------------------------------------------------

template<typename FieldT>
void compare( const FieldT& f1, const FieldT& f2,
	      double& maxAbsErr, double& maxRelErr,
	      double& avgAbsErr, double& avgRelErr )
{
  static const double SMALL = numeric_limits<double>::epsilon() * 10.0;

  maxAbsErr = 0.0;
  maxRelErr = 0.0;
  avgAbsErr = 0.0;
  avgRelErr = 0.0;

  int npts = 0;

  typename FieldT::const_interior_iterator if1 = f1.interior_begin();
  typename FieldT::const_interior_iterator if2 = f2.interior_begin();
  const typename FieldT::const_interior_iterator if1e= f1.interior_end();

  for( ; if1!=if1e; ++if1, ++if2 ){
    const double absErr = std::abs( *if1 - *if2 );
    const double relErr = absErr / (std::abs(*if1)+SMALL);
    maxAbsErr = std::max( maxAbsErr, absErr );
    maxRelErr = std::max( maxRelErr, relErr );
    avgRelErr += relErr;
    avgAbsErr += absErr;
    ++npts;
  }
  avgRelErr /= npts;
  avgAbsErr /= npts;
}

//--------------------------------------------------------------------

template<typename FieldT>
void report_errors( const FieldT& phiExact, const FieldT& phi )
{
  double absErr, relErr, avgAbsErr, avgRelErr;
  compare( phiExact, phi, absErr, relErr, avgAbsErr, avgRelErr );

  cout << scientific << setprecision(5)
       << setw(12) << absErr    << " |"
       << setw(12) << relErr    << " |"
       << setw(12) << avgAbsErr << " |"
       << setw(12) << avgRelErr << " |"
       << endl;
}

//--------------------------------------------------------------------

template<typename OpT, typename FuncType1, typename FuncType2>
void test_grad_op( const Grid& grid,
		   const FuncType1& funcPhi,
		   const FuncType2& funcFPhi )
{
  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const std::vector<int>& dim = grid.extent();

  if( get_n_tot<SrcFieldT>(dim) == 1 || get_n_tot<DestFieldT>(dim) == 1 ) return;

  SrcFieldT  phi      ( get_n_tot<SrcFieldT >(dim),
			get_entries_per_comp<SrcFieldT >(dim),
			get_ghost_set<SrcFieldT >(dim),
			NULL );
  DestFieldT fphi     ( get_n_tot<DestFieldT>(dim),
			get_entries_per_comp<DestFieldT>(dim),
			get_ghost_set<DestFieldT>(dim),
			NULL );
  DestFieldT fphiExact( get_n_tot<DestFieldT>(dim),
			get_entries_per_comp<DestFieldT>(dim),
			get_ghost_set<DestFieldT>(dim),
			NULL );

  switch( FuncType2::FieldType::Location::Dir::value ){
  case XDIR::value :
    funcPhi.evaluate( phi );
    funcFPhi.dx( fphiExact );
    break;
  case YDIR::value :
    funcPhi.evaluate( phi );
    funcFPhi.dy( fphiExact );
    break;
  case ZDIR::value :
    funcPhi.evaluate( phi );
    funcFPhi.dz( fphiExact );
    break;
  default:
    assert(1);
  }

  const OpT* const op = SpatialOpDatabase<OpT>::self().retrieve_operator();

  op->apply_to_field( phi, fphi );

  report_errors( fphiExact, fphi );
}
//--------------------------------------------------------------------

template<typename OpT, typename FuncType1, typename FuncType2>
void test_interp_op( const Grid& grid,
		     const FuncType1& funcPhi,
		     const FuncType2& funcFPhi )
{
  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const std::vector<int>& dim = grid.extent();

  if( get_n_tot<SrcFieldT>(dim) == 1 || get_n_tot<DestFieldT>(dim) == 1 ) return;

  SrcFieldT  phi      ( get_n_tot<SrcFieldT >(dim),
			get_entries_per_comp<SrcFieldT >(dim),
			get_ghost_set<SrcFieldT >(dim),
			NULL );
  DestFieldT fphi     ( get_n_tot<DestFieldT>(dim),
			get_entries_per_comp<DestFieldT>(dim),
			get_ghost_set<DestFieldT>(dim),
			NULL );
  DestFieldT fphiExact( get_n_tot<DestFieldT>(dim),
			get_entries_per_comp<DestFieldT>(dim),
			get_ghost_set<DestFieldT>(dim),
			NULL );

  funcPhi.evaluate( phi );
  funcFPhi.evaluate( fphiExact );

  const OpT* const op = SpatialOpDatabase<OpT>::self().retrieve_operator();
  op->apply_to_field( phi, fphi );

  report_errors( fphiExact, fphi );
}

//--------------------------------------------------------------------

template<typename OpT, typename FuncType1, typename FuncType2>
void test_div_op( const Grid& grid,
		  const FuncType1& funcPhi,
		  const FuncType2& funcFPhi )
{
  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const std::vector<int>& dim = grid.extent();

  if( get_n_tot<SrcFieldT>(dim) == 1 || get_n_tot<DestFieldT>(dim) == 1 ) return;

  SrcFieldT  phi      ( get_n_tot<SrcFieldT >(dim),
			get_entries_per_comp<SrcFieldT >(dim),
			get_ghost_set<SrcFieldT >(dim),
			NULL );
  DestFieldT fphi     ( get_n_tot<DestFieldT>(dim),
			get_entries_per_comp<DestFieldT>(dim),
			get_ghost_set<DestFieldT>(dim),
			NULL );
  DestFieldT fphiExact( get_n_tot<DestFieldT>(dim),
			get_entries_per_comp<DestFieldT>(dim),
			get_ghost_set<DestFieldT>(dim),
			NULL );

  switch( FuncType1::FieldType::Location::Dir::value ){
  case XDIR::value :
    funcPhi.dx( phi );
    funcFPhi.d2x( fphiExact );
    break;
  case YDIR::value :
    funcPhi.dy( phi );
    funcFPhi.d2y( fphiExact );
    break;
  case ZDIR::value :
    funcPhi.dz( phi );
    funcFPhi.d2z( fphiExact );
    break;
  default:
    assert(1);
  }

  const OpT* const op = SpatialOpDatabase<OpT>::self().retrieve_operator();

  op->apply_to_field( phi, fphi );

  report_errors( fphiExact, fphi );

}

//--------------------------------------------------------------------

template<typename OpT, typename Dir>
bool test_bc_helper( const vector<int>&dim,
		     const int ii,
		     const int jj,
		     const int kk,
		     const double bcVal )
{
  using namespace SpatialOps;
  using namespace FVStaggered;

  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const OpT& op = *SpatialOpDatabase<OpT>::self().retrieve_operator();

  SrcFieldT   f( get_n_tot<SrcFieldT >(dim), get_entries_per_comp<SrcFieldT >(dim), get_ghost_set<SrcFieldT >(dim), NULL );
  DestFieldT df( get_n_tot<DestFieldT>(dim), get_entries_per_comp<DestFieldT>(dim), get_ghost_set<DestFieldT>(dim), NULL );

  int icnt=0;
  for( typename SrcFieldT::iterator ifld=f.begin(); ifld!=f.end(); ++ifld,++icnt ) *ifld = icnt;

  // assign the BC.
  assign_bc_point<OpT,Dir>( op, ii, jj, kk, dim, bcVal, f );

  // calculate the dest field
  op.apply_to_field( f, df );

  // verify that the BC was set properly - this is a bit of a hack.
  IndexTriplet st(ii,jj,kk), dt(ii,jj,kk);
  bc_adjust_ijk<SrcFieldT,DestFieldT,Dir>( st, dt, dim );
  
  const int ix = ijk2flat<DestFieldT,DestFieldT::Location::IsSurface,Dir>::value( dim, dt );

  const double abserr = abs(df[ix]-bcVal);
  const double relerr = abserr/abs(bcVal);

  return ( abserr<1.0e-10 && relerr<1.0e-8 );
}

//--------------------------------------------------------------------

void test_bc( const Grid& g )
{
  using namespace SpatialOps;
  using namespace FVStaggered;

  const vector<int>& dim = g.extent();

  cout << endl << "Testing BC setting stuff:" << endl;

  bool isOkay = true;;

  if( dim[0]>1 ){
    // X BCs - Left side
    cout << "  X Dir, (-) side ... ";
    int i=0;
    for( int j=0; j<dim[1]; ++j ){
      for( int k=0; k<dim[2]; ++k ){
	const bool result1 = test_bc_helper<GradSVolSSurfX, XDIR>( dim, i,j,k, 1.2345  );
	const bool result2 = test_bc_helper<InterpSVolSSurfX,XDIR>( dim, i,j,k, 123.456 );
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
	const bool result1 = test_bc_helper<GradSVolSSurfX,  XDIR>( dim, i,j,k, 5.4321 );
	const bool result2 = test_bc_helper<InterpSVolSSurfX,XDIR>( dim, i,j,k, 123.456 );
	if( !result1 || !result2 ) isOkay=false;
      }
    }
    if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;
  }

  if( dim[1]>1 ){
    // Y BCs - Left side
    cout << "  Y Dir, (-) side ... ";
    int j=0;
    for( int i=0; i<dim[0]; ++i ){
      for( int k=0; k<dim[2]; ++k ){
	const bool result1 = test_bc_helper<GradSVolSSurfY, YDIR >( dim, i,j,k, 1.2345 );
	const bool result2 = test_bc_helper<InterpSVolSSurfY,YDIR>( dim, i,j,k, 123.456 );
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
	const bool result1 = test_bc_helper<GradSVolSSurfY,  YDIR>( dim, i,j,k, 5.4321 );
	const bool result2 = test_bc_helper<InterpSVolSSurfY,YDIR>( dim, i,j,k, 123.456 );
	if( !result1 || !result2 ) isOkay=false;
      }
    }
    if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;
  }

  if( dim[2]>1 ){
    // Z BCs - Left side
    cout << "  Z Dir, (-) side ... ";
    int k=0;
    for( int i=0; i<dim[0]; ++i ){
      for( int j=0; j<dim[1]; ++j ){
	const bool result1 = test_bc_helper<GradSVolSSurfZ,  ZDIR>( dim, i,j,k, 1.2345 );
	const bool result2 = test_bc_helper<InterpSVolSSurfZ,ZDIR>( dim, i,j,k, 123.456 );
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
	const bool result1 = test_bc_helper<GradSVolSSurfZ,  ZDIR>( dim, i,j,k, 5.4321 );
	const bool result2 = test_bc_helper<InterpSVolSSurfZ,ZDIR>( dim, i,j,k, 123.456 );
	if( !result1 || !result2 ) isOkay=false;
      }
    }
    if( isOkay ) cout << "PASS" << endl;  else cout << "FAIL" << endl;
  }

  cout << endl;
}

//--------------------------------------------------------------------

void test_ops()
{
  ScratchSVolSSurfX& Sx = *SpatialOpDatabase<ScratchSVolSSurfX>::self().retrieve_operator();
  GradSVolSSurfX&    Gx = *SpatialOpDatabase<GradSVolSSurfX>::self().retrieve_operator();
  Sx.reset_entries(1.0);

//   EpetraExt::RowMatrixToMatrixMarketFile( "Sx.mm", Sx.get_linalg_mat(), "", "" );
//   EpetraExt::RowMatrixToMatrixMarketFile( "Gx.mm", Gx.get_linalg_mat(), "", "" );

  Sx += Gx;
  Sx -= Gx;

  ScratchSVolSSurfY& Sy = *SpatialOpDatabase<ScratchSVolSSurfY>::self().retrieve_operator();
  GradSVolSSurfY&    Gy = *SpatialOpDatabase<GradSVolSSurfY>::self().retrieve_operator();
  Sy += Gy;
}

//--------------------------------------------------------------------

int main()
{
  vector<int> dim(3,1);
  dim[0] = 3 ;
  dim[1] = 3 ;
  dim[2] = 1 ;

  cout << "interior nx = "; cin >> dim[0];
  cout << "interior ny = "; cin >> dim[1];
  cout << "interior nz = "; cin >> dim[2];
  cout << endl;

  std::vector<double> length(3,1);
  std::vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    if( dim[i]>1 ) spacing[i] = length[i]/( get_nx<SVolField>(dim[i]-1) );
  }

  build_ops( dim, spacing );
  const Grid grid( dim, spacing );

  test_bc( grid );
  test_ops();

  // Scalar-Volume to scalar face gradients and laplacians
  {
    SVolField phi( get_n_tot<SVolField>(dim), get_entries_per_comp<SVolField>(dim), get_ghost_set<SVolField>(dim), NULL );

    // sin function
    const SinFun<SVolField  > fun     ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );
    const SinFun<SSurfXField> gradFunX( grid.xcoord_sxsurf(), grid.ycoord_sxsurf(), grid.zcoord_sxsurf() );
    const SinFun<SSurfYField> gradFunY( grid.xcoord_sysurf(), grid.ycoord_sysurf(), grid.zcoord_sysurf() );
    const SinFun<SSurfZField> gradFunZ( grid.xcoord_szsurf(), grid.ycoord_szsurf(), grid.zcoord_szsurf() );
    const SinFun<SVolField  > divFun  ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );
    const SinFun<SSurfXField> interpX ( grid.xcoord_sxsurf(), grid.ycoord_sxsurf(), grid.zcoord_sxsurf()  );
    const SinFun<SSurfYField> interpY ( grid.xcoord_sysurf(), grid.ycoord_sysurf(), grid.zcoord_sysurf()  );
    const SinFun<SSurfZField> interpZ ( grid.xcoord_szsurf(), grid.ycoord_szsurf(), grid.zcoord_szsurf()  );

    cout << "=====================================================" << endl
	 << "Interpolant scalar volume -> scalar surfaces" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpSVolSSurfX>( grid, fun, interpX );
    test_interp_op<InterpSVolSSurfY>( grid, fun, interpY );
    test_interp_op<InterpSVolSSurfZ>( grid, fun, interpZ );
    cout << "=====================================================" << endl << endl;


    cout << "=====================================================" << endl
	 << "Gradient scalar volume -> scalar surfaces" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradSVolSSurfX>( grid, fun, gradFunX );
    test_grad_op<GradSVolSSurfY>( grid, fun, gradFunY );
    test_grad_op<GradSVolSSurfZ>( grid, fun, gradFunZ );
    cout << "=====================================================" << endl << endl;


    cout << "=====================================================" << endl
	 << "Divergence scalar surfaces -> scalar volume" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_div_op<DivSSurfXSVol>( grid, gradFunX, divFun );
    test_div_op<DivSSurfYSVol>( grid, gradFunY, divFun );
    test_div_op<DivSSurfZSVol>( grid, gradFunZ, divFun );
    cout << "=====================================================" << endl << endl;
  }


  {
    const SinFun<SVolField> svolfun( grid.xcoord_svol(), grid.ycoord_svol(), grid.zcoord_svol() );
    const SinFun<XVolField> xvolfun( grid.xcoord_xvol(), grid.ycoord_xvol(), grid.zcoord_xvol() );
    const SinFun<YVolField> yvolfun( grid.xcoord_yvol(), grid.ycoord_yvol(), grid.zcoord_yvol() );
    const SinFun<ZVolField> zvolfun( grid.xcoord_zvol(), grid.ycoord_zvol(), grid.zcoord_zvol() );

    cout << "=====================================================" << endl
	 << "Interpolate scalar volume to staggered volumes" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpSVolXVol>( grid, svolfun, xvolfun );
    test_interp_op<InterpSVolYVol>( grid, svolfun, yvolfun );
    test_interp_op<InterpSVolZVol>( grid, svolfun, zvolfun );
    cout << "=====================================================" << endl << endl;
  }


  // X-volume to x-surface face component gradients
  if( dim[0]>1 ){

    // jcs problems with the y-divergence and z-divergence operators

    // sin function
    const SinFun<XVolField  > sinfun ( grid.xcoord_xvol(),   grid.ycoord_xvol(),   grid.zcoord_xvol()   );
    const SinFun<XSurfXField> gradX  ( grid.xcoord_xxsurf(), grid.ycoord_xxsurf(), grid.zcoord_xxsurf() );
    const SinFun<XSurfYField> gradY  ( grid.xcoord_xysurf(), grid.ycoord_xysurf(), grid.zcoord_xysurf() );
    const SinFun<XSurfZField> gradZ  ( grid.xcoord_xzsurf(), grid.ycoord_xzsurf(), grid.zcoord_xzsurf() );
    const SinFun<XVolField  > divFun ( grid.xcoord_xvol(),   grid.ycoord_xvol(),   grid.zcoord_xvol()   );
    const SinFun<XSurfXField> interpX( grid.xcoord_xxsurf(), grid.ycoord_xxsurf(), grid.zcoord_xxsurf() );
    const SinFun<XSurfYField> interpY( grid.xcoord_xysurf(), grid.ycoord_xysurf(), grid.zcoord_xysurf() );
    const SinFun<XSurfZField> interpZ( grid.xcoord_xzsurf(), grid.ycoord_xzsurf(), grid.zcoord_xzsurf() );
//     const SinFun<XSurfField > interp ( grid.xcoord_xsurf(),  grid.ycoord_xsurf(),  grid.zcoord_xsurf()  );

    cout << "=====================================================" << endl
	 << "Gradient x-volume -> x-surface" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradXVolXSurfX>( grid, sinfun, gradX );
    test_grad_op<GradXVolXSurfY>( grid, sinfun, gradY );
    test_grad_op<GradXVolXSurfZ>( grid, sinfun, gradZ );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
	 << "Divergence x-surface -> x-volume" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_div_op<DivXSurfXXVol>( grid, gradX, divFun );
    test_div_op<DivXSurfYXVol>( grid, gradY, divFun );
    test_div_op<DivXSurfZXVol>( grid, gradZ, divFun );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
	 << "Interpolate x-volume -> x-surfaces" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpXVolXSurfX>( grid, sinfun, interpX );
    test_interp_op<InterpXVolXSurfY>( grid, sinfun, interpY );
    test_interp_op<InterpXVolXSurfZ>( grid, sinfun, interpZ );
//     test_interp_op<InterpXVolXSurf >( grid, sinfun, interp  );
    cout << "=====================================================" << endl << endl;
  }


  // Y-volume to y-surface face component gradients
  if( dim[1]>1 ){

    // jcs we seem to have a problem with the x-divergence and z-divergence operators.

    // sin function
    const SinFun<YVolField  > sinfun ( grid.xcoord_yvol(),   grid.ycoord_yvol(),   grid.zcoord_yvol()   );
    const SinFun<YSurfXField> gradX  ( grid.xcoord_yxsurf(), grid.ycoord_yxsurf(), grid.zcoord_yxsurf() );
    const SinFun<YSurfYField> gradY  ( grid.xcoord_yysurf(), grid.ycoord_yysurf(), grid.zcoord_yysurf() );
    const SinFun<YSurfZField> gradZ  ( grid.xcoord_yzsurf(), grid.ycoord_yzsurf(), grid.zcoord_yzsurf() );
    const SinFun<YVolField  > divFun ( grid.xcoord_yvol(),   grid.ycoord_yvol(),   grid.zcoord_yvol()   );
    const SinFun<YSurfXField> interpX( grid.xcoord_yxsurf(), grid.ycoord_yxsurf(), grid.zcoord_yxsurf() );
    const SinFun<YSurfYField> interpY( grid.xcoord_yysurf(), grid.ycoord_yysurf(), grid.zcoord_yysurf() );
    const SinFun<YSurfZField> interpZ( grid.xcoord_yzsurf(), grid.ycoord_yzsurf(), grid.zcoord_yzsurf() );
    //    const SinFun<YSurfField > interp ( grid.xcoord_ysurf(),  grid.ycoord_ysurf(),  grid.zcoord_ysurf()  );

    cout << "=====================================================" << endl
	 << "Gradient y-volume -> y-surface" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradYVolYSurfX>( grid, sinfun, gradX );
    test_grad_op<GradYVolYSurfY>( grid, sinfun, gradY );
    test_grad_op<GradYVolYSurfZ>( grid, sinfun, gradZ );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
	 << "Divergence y-surface -> y-volume" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_div_op<DivYSurfXYVol>( grid, gradX, divFun );
    test_div_op<DivYSurfYYVol>( grid, gradY, divFun );
    test_div_op<DivYSurfZYVol>( grid, gradZ, divFun );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
	 << "Interpolate y-volume -> y-surfaces" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpYVolYSurfX>( grid, sinfun, interpX );
    test_interp_op<InterpYVolYSurfY>( grid, sinfun, interpY );
    test_interp_op<InterpYVolYSurfZ>( grid, sinfun, interpZ );
//     test_interp_op<InterpYVolYSurf >( grid, sinfun, interp  );
    cout << "=====================================================" << endl << endl;
  }


  // Z-volume to z-surface face component gradients
  if( dim[2]>1 ){

    // jcs problems with the x-gradient and z-gradient operators.
    // jcs problems with all of the divergence operators.
    // jcs problems with all of the interpolant operators

    // sin function
    const SinFun<ZVolField  > sinfun ( grid.xcoord_zvol(),   grid.ycoord_zvol(),   grid.zcoord_zvol()   );
    const SinFun<ZSurfXField> gradX  ( grid.xcoord_zxsurf(), grid.ycoord_zxsurf(), grid.zcoord_zxsurf() );
    const SinFun<ZSurfYField> gradY  ( grid.xcoord_zysurf(), grid.ycoord_zysurf(), grid.zcoord_zysurf() );
    const SinFun<ZSurfZField> gradZ  ( grid.xcoord_zzsurf(), grid.ycoord_zzsurf(), grid.zcoord_zzsurf() );
    const SinFun<ZVolField  > divFun ( grid.xcoord_zvol(),   grid.ycoord_zvol(),   grid.zcoord_zvol()   );
    const SinFun<ZSurfXField> interpX( grid.xcoord_zxsurf(), grid.ycoord_zxsurf(), grid.zcoord_zxsurf() );
    const SinFun<ZSurfYField> interpY( grid.xcoord_zysurf(), grid.ycoord_zysurf(), grid.zcoord_zysurf() );
    const SinFun<ZSurfZField> interpZ( grid.xcoord_zzsurf(), grid.ycoord_zzsurf(), grid.zcoord_zzsurf() );
//     const SinFun<ZSurfField > interp ( grid.xcoord_zsurf(),  grid.ycoord_zsurf(),  grid.zcoord_zsurf()  );

    cout << "=====================================================" << endl
	 << "Gradient z-volume -> z-surface" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradZVolZSurfX>( grid, sinfun, gradX );
    test_grad_op<GradZVolZSurfY>( grid, sinfun, gradY );
    test_grad_op<GradZVolZSurfZ>( grid, sinfun, gradZ );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
	 << "Divergence z-surface -> z-volume" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_div_op<DivZSurfXZVol>( grid, gradX, divFun );
    test_div_op<DivZSurfYZVol>( grid, gradY, divFun );
    test_div_op<DivZSurfZZVol>( grid, gradZ, divFun );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
	 << "Interpolate z-volume -> z-surfaces" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpZVolZSurfX>( grid, sinfun, interpX );
    test_interp_op<InterpZVolZSurfY>( grid, sinfun, interpY );
    test_interp_op<InterpZVolZSurfZ>( grid, sinfun, interpZ );
//     test_interp_op<InterpZVolZSurf >( grid, sinfun, interp  );
    cout << "=====================================================" << endl << endl;
  }




  {
    cout << endl << "Scalar volume scratch operators" << endl
	 << " max abs err | max rel err | avg abs err | avg rel err |" << endl
	 << "-------------|-------------|-------------|-------------|" << endl;
    ScratchSVol* const Ssx = SpatialOpDatabase<ScratchSVol>::self().retrieve_operator( ssxid );
    ScratchSVol* const Ssy = SpatialOpDatabase<ScratchSVol>::self().retrieve_operator( ssyid );
    ScratchSVol* const Ssz = SpatialOpDatabase<ScratchSVol>::self().retrieve_operator( sszid );

    const GradSVolSSurfX* const Gsx = SpatialOpDatabase<GradSVolSSurfX>::self().retrieve_operator();
    const GradSVolSSurfY* const Gsy = SpatialOpDatabase<GradSVolSSurfY>::self().retrieve_operator();
    const GradSVolSSurfZ* const Gsz = SpatialOpDatabase<GradSVolSSurfZ>::self().retrieve_operator();

    const DivSSurfXSVol * const Dsx = SpatialOpDatabase<DivSSurfXSVol>::self().retrieve_operator();
    const DivSSurfYSVol * const Dsy = SpatialOpDatabase<DivSSurfYSVol>::self().retrieve_operator();
    const DivSSurfZSVol * const Dsz = SpatialOpDatabase<DivSSurfZSVol>::self().retrieve_operator();

    Dsx->apply_to_op( *Gsx, *Ssx );
    Dsy->apply_to_op( *Gsy, *Ssy );
    Dsz->apply_to_op( *Gsz, *Ssz );

    const SinFun<SVolField  > fun     ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );
    const SinFun<SVolField  > divFun  ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );

    SVolField phi       ( get_n_tot<SVolField>(dim), get_entries_per_comp<SVolField>(dim), get_ghost_set<SVolField>(dim), NULL );
    SVolField d2phi     ( get_n_tot<SVolField>(dim), get_entries_per_comp<SVolField>(dim), get_ghost_set<SVolField>(dim), NULL );
    SVolField d2phiExact( get_n_tot<SVolField>(dim), get_entries_per_comp<SVolField>(dim), get_ghost_set<SVolField>(dim), NULL );

    fun.evaluate( phi );

    if( dim[0]>1 ){
      Ssx->apply_to_field( phi, d2phi );
      divFun.d2x( d2phiExact );
      report_errors( d2phiExact, d2phi );
    }
    if( dim[1]>1 ){
      Ssy->apply_to_field( phi, d2phi );
      divFun.d2y( d2phiExact );
      report_errors( d2phiExact, d2phi );
    }
    if( dim[2]>1 ){
      Ssz->apply_to_field( phi, d2phi );
      divFun.d2z( d2phiExact );
      report_errors( d2phiExact, d2phi );
    }
  }

  {
    cout << endl << "x-volume scratch operators ... ";
    ScratchXVol* const Sxx = SpatialOpDatabase<ScratchXVol>::self().retrieve_operator( xsxid );
    ScratchXVol* const Sxy = SpatialOpDatabase<ScratchXVol>::self().retrieve_operator( xsyid );
    ScratchXVol* const Sxz = SpatialOpDatabase<ScratchXVol>::self().retrieve_operator( xszid );

    const GradXVolXSurfX* const Gsx = SpatialOpDatabase<GradXVolXSurfX>::self().retrieve_operator();
    const GradXVolXSurfY* const Gsy = SpatialOpDatabase<GradXVolXSurfY>::self().retrieve_operator();
    const GradXVolXSurfZ* const Gsz = SpatialOpDatabase<GradXVolXSurfZ>::self().retrieve_operator();

    const DivXSurfXXVol * const Dsx = SpatialOpDatabase<DivXSurfXXVol>::self().retrieve_operator();
    const DivXSurfYXVol * const Dsy = SpatialOpDatabase<DivXSurfYXVol>::self().retrieve_operator();
    const DivXSurfZXVol * const Dsz = SpatialOpDatabase<DivXSurfZXVol>::self().retrieve_operator();

    Dsx->apply_to_op( *Gsx, *Sxx );
    Dsy->apply_to_op( *Gsy, *Sxy );
    Dsz->apply_to_op( *Gsz, *Sxz );

//     EpetraExt::RowMatrixToMatrixMarketFile( "Sxx.mm", Sxx->get_linalg_mat(), "", "" );
//     EpetraExt::RowMatrixToMatrixMarketFile( "Sxy.mm", Sxy->get_linalg_mat(), "", "" );
//     EpetraExt::RowMatrixToMatrixMarketFile( "Sxz.mm", Sxz->get_linalg_mat(), "", "" );

    cout << "done" << endl;
  }

}
