#include <FVStaggered.h>
#include <buildOps.h>

using namespace SpatialOps;
using namespace FVStaggered;

void build_ops( const std::vector<int>& dim,
		const std::vector<double>& spacing,
		const std::vector<bool>& bcPlus )
{
  using SpatialOps::SpatialOpDatabase;

  std::vector<double> area(3,1);
  area[0] = spacing[1]*spacing[2];
  area[1] = spacing[0]*spacing[2];
  area[2] = spacing[0]*spacing[1];
  const double vol = spacing[0]*spacing[1]*spacing[2];

  //--------------------------------------------------------
  // Divergence operators 
  //--------------------------------------------------------
  {
    DivSSurfXSVol::Assembler Dssxsva( dim, area[0], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivSSurfYSVol::Assembler Dssysva( dim, area[1], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivSSurfZSVol::Assembler Dsszsva( dim, area[2], vol, bcPlus[0], bcPlus[1], bcPlus[2] );

    DivXSurfXXVol::Assembler Dxsxxva( dim, area[0], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivXSurfYXVol::Assembler Dxsyxva( dim, area[1], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivXSurfZXVol::Assembler Dxszxva( dim, area[2], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
						       											
    DivYSurfXYVol::Assembler Dysxyva( dim, area[0], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivYSurfYYVol::Assembler Dysyyva( dim, area[1], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivYSurfZYVol::Assembler Dyszyva( dim, area[2], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
						       											
    DivZSurfXZVol::Assembler Dzsxzva( dim, area[0], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivZSurfYZVol::Assembler Dzsyzva( dim, area[1], vol, bcPlus[0], bcPlus[1], bcPlus[2] );
    DivZSurfZZVol::Assembler Dzszzva( dim, area[2], vol, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<DivSSurfXSVol>::self().register_new_operator( new DivSSurfXSVol(Dssxsva) ); //SpatialOpDatabase<DivSSurfXSVol>::self().retrieve_operator()->write_matlab("Dsxsv");
    SpatialOpDatabase<DivSSurfYSVol>::self().register_new_operator( new DivSSurfYSVol(Dssysva) ); //SpatialOpDatabase<DivSSurfYSVol>::self().retrieve_operator()->write_matlab("Dsysv");
    SpatialOpDatabase<DivSSurfZSVol>::self().register_new_operator( new DivSSurfZSVol(Dsszsva) );

    SpatialOpDatabase<DivXSurfXXVol>::self().register_new_operator( new DivXSurfXXVol(Dxsxxva) ); //SpatialOpDatabase<DivXSurfXXVol>::self().retrieve_operator()->write_matlab("Dxsxxv");
    SpatialOpDatabase<DivXSurfYXVol>::self().register_new_operator( new DivXSurfYXVol(Dxsyxva) ); //SpatialOpDatabase<DivXSurfYXVol>::self().retrieve_operator()->write_matlab("Dxsyxv");
    SpatialOpDatabase<DivXSurfZXVol>::self().register_new_operator( new DivXSurfZXVol(Dxszxva) );

    SpatialOpDatabase<DivYSurfXYVol>::self().register_new_operator( new DivYSurfXYVol(Dysxyva) ); //SpatialOpDatabase<DivYSurfXYVol>::self().retrieve_operator()->write_matlab("Dysxyv");
    SpatialOpDatabase<DivYSurfYYVol>::self().register_new_operator( new DivYSurfYYVol(Dysyyva) ); //SpatialOpDatabase<DivYSurfXYVol>::self().retrieve_operator()->write_matlab("Dysyyv");
    SpatialOpDatabase<DivYSurfZYVol>::self().register_new_operator( new DivYSurfZYVol(Dyszyva) ); //SpatialOpDatabase<DivYSurfXYVol>::self().retrieve_operator()->write_matlab("Dyszyv");

    SpatialOpDatabase<DivZSurfXZVol>::self().register_new_operator( new DivZSurfXZVol(Dzsxzva) );
    SpatialOpDatabase<DivZSurfYZVol>::self().register_new_operator( new DivZSurfYZVol(Dzsyzva) );
    SpatialOpDatabase<DivZSurfZZVol>::self().register_new_operator( new DivZSurfZZVol(Dzszzva) ); //SpatialOpDatabase<DivZSurfZZVol>::self().retrieve_operator()->write_matlab("Dzszzv");
  }

  //--------------------------------------------------------
  // gradient operators - diffusive fluxes
  //--------------------------------------------------------
  {
    GradSVolSSurfX::Assembler Gsvssxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradSVolSSurfY::Assembler Gsvssya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradSVolSSurfZ::Assembler Gsvssza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
						      											
    GradXVolXSurfX::Assembler Gxvxsxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradXVolXSurfY::Assembler Gxvxsya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradXVolXSurfZ::Assembler Gxvxsza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
						      											
    GradYVolYSurfX::Assembler Gyvysxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradYVolYSurfY::Assembler Gyvysya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradYVolYSurfZ::Assembler Gyvysza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    GradZVolZSurfX::Assembler Gzvzsxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradZVolZSurfY::Assembler Gzvzsya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradZVolZSurfZ::Assembler Gzvzsza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<GradSVolSSurfX>::self().register_new_operator( new GradSVolSSurfX(Gsvssxa) ); //SpatialOpDatabase<GradSVolSSurfX>::self().retrieve_operator()->write_matlab("Gssx");
    SpatialOpDatabase<GradSVolSSurfY>::self().register_new_operator( new GradSVolSSurfY(Gsvssya) ); //SpatialOpDatabase<GradSVolSSurfY>::self().retrieve_operator()->write_matlab("Gssy");
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

  //--------------------------------------------------------
  // interpolant scalar volume to scalar surface (diffusivities)
  //--------------------------------------------------------
  {
    InterpSVolSSurfX::Assembler Rsvssxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolSSurfY::Assembler Rsvssya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolSSurfZ::Assembler Rsvssza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<InterpSVolSSurfX>::self().register_new_operator( new InterpSVolSSurfX(Rsvssxa) );
    SpatialOpDatabase<InterpSVolSSurfY>::self().register_new_operator( new InterpSVolSSurfY(Rsvssya) );
    SpatialOpDatabase<InterpSVolSSurfZ>::self().register_new_operator( new InterpSVolSSurfZ(Rsvssza) );
  }

  //--------------------------------------------------------
  // interpolants - scalar volume to staggered surfaces (viscosity, dilatation)
  //--------------------------------------------------------
  {
    InterpSVolXSurfX::Assembler Rsvxsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolXSurfY::Assembler Rsvxsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolXSurfZ::Assembler Rsvxsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    SpatialOpDatabase<InterpSVolXSurfX>::self().register_new_operator( new InterpSVolXSurfX(Rsvxsxa) );  //SpatialOpDatabase<InterpSVolXSurfX>::self().retrieve_operator()->write_matlab("Rsvxsx");
    SpatialOpDatabase<InterpSVolXSurfY>::self().register_new_operator( new InterpSVolXSurfY(Rsvxsya) );  //SpatialOpDatabase<InterpSVolXSurfY>::self().retrieve_operator()->write_matlab("Rsvxsy");
    SpatialOpDatabase<InterpSVolXSurfZ>::self().register_new_operator( new InterpSVolXSurfZ(Rsvxsza) );  //SpatialOpDatabase<InterpSVolXSurfZ>::self().retrieve_operator()->write_matlab("Rsvxsz");

    InterpSVolYSurfX::Assembler Rsvysxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolYSurfY::Assembler Rsvysya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolYSurfZ::Assembler Rsvysza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    SpatialOpDatabase<InterpSVolYSurfX>::self().register_new_operator( new InterpSVolYSurfX(Rsvysxa) );  //SpatialOpDatabase<InterpSVolYSurfX>::self().retrieve_operator()->write_matlab("Rsvysx");
    SpatialOpDatabase<InterpSVolYSurfY>::self().register_new_operator( new InterpSVolYSurfY(Rsvysya) );  //SpatialOpDatabase<InterpSVolYSurfY>::self().retrieve_operator()->write_matlab("Rsvysy");
    SpatialOpDatabase<InterpSVolYSurfZ>::self().register_new_operator( new InterpSVolYSurfZ(Rsvysza) );  //SpatialOpDatabase<InterpSVolYSurfZ>::self().retrieve_operator()->write_matlab("Rsvysz");

    InterpSVolZSurfX::Assembler Rsvzsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolZSurfY::Assembler Rsvzsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolZSurfZ::Assembler Rsvzsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    SpatialOpDatabase<InterpSVolZSurfX>::self().register_new_operator( new InterpSVolZSurfX(Rsvzsxa) );  //SpatialOpDatabase<InterpSVolZSurfX>::self().retrieve_operator()->write_matlab("Rsvzsx");
    SpatialOpDatabase<InterpSVolZSurfY>::self().register_new_operator( new InterpSVolZSurfY(Rsvzsya) );  //SpatialOpDatabase<InterpSVolZSurfY>::self().retrieve_operator()->write_matlab("Rsvzsy");
    SpatialOpDatabase<InterpSVolZSurfZ>::self().register_new_operator( new InterpSVolZSurfZ(Rsvzsza) );  //SpatialOpDatabase<InterpSVolZSurfZ>::self().retrieve_operator()->write_matlab("Rsvzsz");
  }

  //--------------------------------------------------------
  // interpolants - scalar volume to staggered volume (density)
  //--------------------------------------------------------
  {
    InterpSVolXVol::Assembler Rsvxva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolYVol::Assembler Rsvyva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolZVol::Assembler Rsvzva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<InterpSVolXVol>::self().register_new_operator( new InterpSVolXVol(Rsvxva) );  //SpatialOpDatabase<InterpSVolXVol>::self().retrieve_operator()->write_matlab("Rsvxv");
    SpatialOpDatabase<InterpSVolYVol>::self().register_new_operator( new InterpSVolYVol(Rsvyva) );  //SpatialOpDatabase<InterpSVolYVol>::self().retrieve_operator()->write_matlab("Rsvyv");
    SpatialOpDatabase<InterpSVolZVol>::self().register_new_operator( new InterpSVolZVol(Rsvzva) );  //SpatialOpDatabase<InterpSVolZVol>::self().retrieve_operator()->write_matlab("Rsvzv");
  }

  //--------------------------------------------------------
  // interpolants - staggered volume to staggered surfaces (advecting velocities)
  //--------------------------------------------------------
  {
    InterpXVolYSurfX::Assembler Rxvysxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpXVolZSurfX::Assembler Rxvzsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    InterpYVolXSurfY::Assembler Ryvxsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpYVolZSurfY::Assembler Ryvzsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    InterpZVolXSurfZ::Assembler Rzvxsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpZVolYSurfZ::Assembler Rzvysza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<InterpXVolYSurfX>::self().register_new_operator( new InterpXVolYSurfX(Rxvysxa) );   //SpatialOpDatabase<InterpXVolYSurfX>::self().retrieve_operator()->write_matlab("Rxvysx");
    SpatialOpDatabase<InterpXVolZSurfX>::self().register_new_operator( new InterpXVolZSurfX(Rxvzsxa) );	  //

    SpatialOpDatabase<InterpYVolXSurfY>::self().register_new_operator( new InterpYVolXSurfY(Ryvxsya) );
    SpatialOpDatabase<InterpYVolZSurfY>::self().register_new_operator( new InterpYVolZSurfY(Ryvzsya) );

    SpatialOpDatabase<InterpZVolXSurfZ>::self().register_new_operator( new InterpZVolXSurfZ(Rzvxsza) );
    SpatialOpDatabase<InterpZVolYSurfZ>::self().register_new_operator( new InterpZVolYSurfZ(Rzvysza) );

    // interpolants - volume to surface for staggered cells.
    {
      InterpXVolXSurfX::Assembler Rxvxsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpXVolXSurfY::Assembler Rxvxsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpXVolXSurfZ::Assembler Rxvxsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

      SpatialOpDatabase<InterpXVolXSurfX>::self().register_new_operator( new InterpXVolXSurfX(Rxvxsxa) );  //SpatialOpDatabase<InterpXVolXSurfX>::self().retrieve_operator()->write_matlab("Rxvxsx");
      SpatialOpDatabase<InterpXVolXSurfY>::self().register_new_operator( new InterpXVolXSurfY(Rxvxsya) );  //SpatialOpDatabase<InterpXVolXSurfY>::self().retrieve_operator()->write_matlab("Rxvxsy");
      SpatialOpDatabase<InterpXVolXSurfZ>::self().register_new_operator( new InterpXVolXSurfZ(Rxvxsza) );  //SpatialOpDatabase<InterpXVolXSurfZ>::self().retrieve_operator()->write_matlab("Rxvxsz");


      InterpYVolYSurfX::Assembler Ryvysxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpYVolYSurfY::Assembler Ryvysya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpYVolYSurfZ::Assembler Ryvysza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

      SpatialOpDatabase<InterpYVolYSurfX>::self().register_new_operator( new InterpYVolYSurfX(Ryvysxa) );
      SpatialOpDatabase<InterpYVolYSurfY>::self().register_new_operator( new InterpYVolYSurfY(Ryvysya) );
      SpatialOpDatabase<InterpYVolYSurfZ>::self().register_new_operator( new InterpYVolYSurfZ(Ryvysza) );


      InterpZVolZSurfX::Assembler Rzvzsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpZVolZSurfY::Assembler Rzvzsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpZVolZSurfZ::Assembler Rzvzsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

      SpatialOpDatabase<InterpZVolZSurfX>::self().register_new_operator( new InterpZVolZSurfX(Rzvzsxa) );  //SpatialOpDatabase<InterpZVolZSurfX>::self().retrieve_operator()->write_matlab("Rzvzsx");
      SpatialOpDatabase<InterpZVolZSurfY>::self().register_new_operator( new InterpZVolZSurfY(Rzvzsya) );  //
      SpatialOpDatabase<InterpZVolZSurfZ>::self().register_new_operator( new InterpZVolZSurfZ(Rzvzsza) );  //SpatialOpDatabase<InterpZVolZSurfZ>::self().retrieve_operator()->write_matlab("Rzvzsz");
    }
  }

  //--------------------------------------------------------
  // interpolants - staggard volume to scalar surface  (advecting velocities)
  //--------------------------------------------------------
  {
    InterpXVolSSurfX::Assembler Rxvssxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpYVolSSurfY::Assembler Ryvssya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpZVolSSurfZ::Assembler Rzvssza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<InterpXVolSSurfX>::self().register_new_operator( new InterpXVolSSurfX(Rxvssxa) );
    SpatialOpDatabase<InterpYVolSSurfY>::self().register_new_operator( new InterpYVolSSurfY(Ryvssya) );
    SpatialOpDatabase<InterpZVolSSurfZ>::self().register_new_operator( new InterpZVolSSurfZ(Rzvssza) );
  }


  //--------------------------------------------------------
  // scratch operators
  //--------------------------------------------------------
  {
    ScratchSVol::Assembler Ssvxa( dim, XDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSVol::Assembler Ssvya( dim, YDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSVol::Assembler Ssvza( dim, ZDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );

    ScratchXVol::Assembler Sxvxa( dim, XDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchXVol::Assembler Sxvya( dim, YDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchXVol::Assembler Sxvza( dim, ZDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );

    ScratchYVol::Assembler Syvxa( dim, XDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchYVol::Assembler Syvya( dim, YDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchYVol::Assembler Syvza( dim, ZDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );

    ScratchZVol::Assembler Szvxa( dim, XDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchZVol::Assembler Szvya( dim, YDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchZVol::Assembler Szvza( dim, ZDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    
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
 
    ScratchSVolSSurfX::Assembler Ssvssxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSVolSSurfY::Assembler Ssvssya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSVolSSurfZ::Assembler Ssvssza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    ScratchXVolXSurfX::Assembler Sxvxsxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchXVolXSurfY::Assembler Sxvxsya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchXVolXSurfZ::Assembler Sxvxsza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    ScratchYVolYSurfX::Assembler Syvysxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchYVolYSurfY::Assembler Syvysya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchYVolYSurfZ::Assembler Syvysza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    ScratchZVolZSurfX::Assembler Szvzsxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchZVolZSurfY::Assembler Szvzsya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchZVolZSurfZ::Assembler Szvzsza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );

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
    ScratchSSurfX::Assembler Sssxa( dim, XDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSSurfY::Assembler Sssya( dim, YDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSSurfZ::Assembler Sssza( dim, ZDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );

    SpatialOpDatabase<ScratchSSurfX>::self().register_new_operator( new ScratchSSurfX(Sssxa) );
    SpatialOpDatabase<ScratchSSurfY>::self().register_new_operator( new ScratchSSurfY(Sssya) );
    SpatialOpDatabase<ScratchSSurfZ>::self().register_new_operator( new ScratchSSurfZ(Sssza) );
  }
}
