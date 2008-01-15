#include <FVStaggered.h>
#include <buildOps.h>

using namespace SpatialOps;
using namespace FVStaggered;

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
    InterpSVolXSurfX::Assembler Rsvxsxa( dim );
    InterpSVolXSurfY::Assembler Rsvxsya( dim );
    InterpSVolXSurfZ::Assembler Rsvxsza( dim );
    SpatialOpDatabase<InterpSVolXSurfX>::self().register_new_operator( new InterpSVolXSurfX(Rsvxsxa) );
    SpatialOpDatabase<InterpSVolXSurfY>::self().register_new_operator( new InterpSVolXSurfY(Rsvxsya) );
    SpatialOpDatabase<InterpSVolXSurfZ>::self().register_new_operator( new InterpSVolXSurfZ(Rsvxsza) );

    InterpSVolYSurfX::Assembler Rsvysxa( dim );
    InterpSVolYSurfY::Assembler Rsvysya( dim );
    InterpSVolYSurfZ::Assembler Rsvysza( dim );
    SpatialOpDatabase<InterpSVolYSurfX>::self().register_new_operator( new InterpSVolYSurfX(Rsvysxa) );
    SpatialOpDatabase<InterpSVolYSurfY>::self().register_new_operator( new InterpSVolYSurfY(Rsvysya) );
    SpatialOpDatabase<InterpSVolYSurfZ>::self().register_new_operator( new InterpSVolYSurfZ(Rsvysza) );

    InterpSVolZSurfX::Assembler Rsvzsxa( dim );
    InterpSVolZSurfY::Assembler Rsvzsya( dim );
    InterpSVolZSurfZ::Assembler Rsvzsza( dim );
    SpatialOpDatabase<InterpSVolZSurfX>::self().register_new_operator( new InterpSVolZSurfX(Rsvzsxa) );
    SpatialOpDatabase<InterpSVolZSurfY>::self().register_new_operator( new InterpSVolZSurfY(Rsvzsya) );
    SpatialOpDatabase<InterpSVolZSurfZ>::self().register_new_operator( new InterpSVolZSurfZ(Rsvzsza) );
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
    
    SpatialOpDatabase<ScratchSVol>::self().register_new_operator( new ScratchSVol(Ssvxa) );
    SpatialOpDatabase<ScratchSVol>::self().register_new_operator( new ScratchSVol(Ssvya) );
    SpatialOpDatabase<ScratchSVol>::self().register_new_operator( new ScratchSVol(Ssvza) );

    SpatialOpDatabase<ScratchXVol>::self().register_new_operator( new ScratchXVol(Sxvxa) );
    SpatialOpDatabase<ScratchXVol>::self().register_new_operator( new ScratchXVol(Sxvya) );
    SpatialOpDatabase<ScratchXVol>::self().register_new_operator( new ScratchXVol(Sxvza) );

    SpatialOpDatabase<ScratchYVol>::self().register_new_operator( new ScratchYVol(Syvxa) );
    SpatialOpDatabase<ScratchYVol>::self().register_new_operator( new ScratchYVol(Syvya) );
    SpatialOpDatabase<ScratchYVol>::self().register_new_operator( new ScratchYVol(Syvza) );

    SpatialOpDatabase<ScratchZVol>::self().register_new_operator( new ScratchZVol(Szvxa) );
    SpatialOpDatabase<ScratchZVol>::self().register_new_operator( new ScratchZVol(Szvya) );
    SpatialOpDatabase<ScratchZVol>::self().register_new_operator( new ScratchZVol(Szvza) );


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
