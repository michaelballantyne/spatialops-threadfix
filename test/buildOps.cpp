#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FVStaggered.h>

#include "buildOps.h"

using namespace SpatialOps;
using namespace structured;

void build_ops( const IntVec& dim,
                const std::vector<double>& spacing,
                const std::vector<bool>& bcPlus,
                OperatorDatabase& opDB )
{
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

    opDB.register_new_operator<DivSSurfXSVol>( Dssxsva ); //opDB.retrieve_operator<DivSSurfXSVol>()->write_matlab("Dsxsv");
    opDB.register_new_operator<DivSSurfYSVol>( Dssysva ); //opDB.retrieve_operator<DivSSurfYSVol>()->write_matlab("Dsysv");
    opDB.register_new_operator<DivSSurfZSVol>( Dsszsva );

    opDB.register_new_operator<DivXSurfXXVol>( Dxsxxva ); //opDB.retrieve_operator<DivXSurfXXVol>()->write_matlab("Dxsxxv");
    opDB.register_new_operator<DivXSurfYXVol>( Dxsyxva ); //opDB.retrieve_operator<DivXSurfYXVol>()->write_matlab("Dxsyxv");
    opDB.register_new_operator<DivXSurfZXVol>( Dxszxva );

    opDB.register_new_operator<DivYSurfXYVol>( Dysxyva ); //opDB.retrieve_operator<DivYSurfXYVol>()->write_matlab("Dysxyv");
    opDB.register_new_operator<DivYSurfYYVol>( Dysyyva ); //opDB.retrieve_operator<DivYSurfXYVol>()->write_matlab("Dysyyv");
    opDB.register_new_operator<DivYSurfZYVol>( Dyszyva ); //opDB.retrieve_operator<DivYSurfXYVol>()->write_matlab("Dyszyv");

    opDB.register_new_operator<DivZSurfXZVol>( Dzsxzva );
    opDB.register_new_operator<DivZSurfYZVol>( Dzsyzva );
    opDB.register_new_operator<DivZSurfZZVol>( Dzszzva ); //opDB.retrieve_operator<DivZSurfZZVol>()->write_matlab("Dzszzv");
  }

  //--------------------------------------------------------
  // gradient operators - diffusive fluxes
  //--------------------------------------------------------
  {
    GradSVolSSurfX::Assembler Gsvssxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradSVolSSurfY::Assembler Gsvssya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradSVolSSurfZ::Assembler Gsvssza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradSVolSSurfX>( Gsvssxa ); //opDB.retrieve_operator<GradSVolSSurfX>()->write_matlab("Gssx");
    opDB.register_new_operator<GradSVolSSurfY>( Gsvssya ); //opDB.retrieve_operator<GradSVolSSurfY>()->write_matlab("Gssy");
    opDB.register_new_operator<GradSVolSSurfZ>( Gsvssza );


    GradXVolXSurfX::Assembler Gxvxsxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradXVolXSurfY::Assembler Gxvxsya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradXVolXSurfZ::Assembler Gxvxsza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradXVolXSurfX>( Gxvxsxa );
    opDB.register_new_operator<GradXVolXSurfY>( Gxvxsya );
    opDB.register_new_operator<GradXVolXSurfZ>( Gxvxsza );


    GradYVolYSurfX::Assembler Gyvysxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradYVolYSurfY::Assembler Gyvysya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradYVolYSurfZ::Assembler Gyvysza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradYVolYSurfX>( new GradYVolYSurfX(Gyvysxa) );
    opDB.register_new_operator<GradYVolYSurfY>( new GradYVolYSurfY(Gyvysya) );
    opDB.register_new_operator<GradYVolYSurfZ>( new GradYVolYSurfZ(Gyvysza) );


    GradZVolZSurfX::Assembler Gzvzsxa( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradZVolZSurfY::Assembler Gzvzsya( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradZVolZSurfZ::Assembler Gzvzsza( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradZVolZSurfX>( new GradZVolZSurfX(Gzvzsxa) );
    opDB.register_new_operator<GradZVolZSurfY>( new GradZVolZSurfY(Gzvzsya) );
    opDB.register_new_operator<GradZVolZSurfZ>( new GradZVolZSurfZ(Gzvzsza) );


    GradSVolXVol  ::Assembler Gsvxva ( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );  
    GradSVolYVol  ::Assembler Gsvyva ( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradSVolZVol  ::Assembler Gsvzva ( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradSVolXVol>( new GradSVolXVol(Gsvxva) );
    opDB.register_new_operator<GradSVolYVol>( new GradSVolYVol(Gsvyva) );
    opDB.register_new_operator<GradSVolZVol>( new GradSVolZVol(Gsvzva) );


    GradXVolYSurfX::Assembler Gxvysxa( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradXVolZSurfX::Assembler Gxvzsxa( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradXVolYSurfX>( new GradXVolYSurfX(Gxvysxa) );
    opDB.register_new_operator<GradXVolZSurfX>( new GradXVolZSurfX(Gxvzsxa) );


    GradYVolXSurfY::Assembler Gyvxsya( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradYVolZSurfY::Assembler Gyvzsya( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradYVolXSurfY>( new GradYVolXSurfY(Gyvxsya) );
    opDB.register_new_operator<GradYVolZSurfY>( new GradYVolZSurfY(Gyvzsya) );


    GradZVolXSurfZ::Assembler Gzvxsza( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradZVolYSurfZ::Assembler Gzvysza( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradZVolXSurfZ>( new GradZVolXSurfZ(Gzvxsza) );
    opDB.register_new_operator<GradZVolYSurfZ>( new GradZVolYSurfZ(Gzvysza) );

    // dilatation
    GradXVolSVol::Assembler Gxvsva( spacing[0], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradYVolSVol::Assembler Gyvsva( spacing[1], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    GradZVolSVol::Assembler Gzvsva( spacing[2], dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<GradXVolSVol>( new GradXVolSVol(Gxvsva) );
    opDB.register_new_operator<GradYVolSVol>( new GradYVolSVol(Gyvsva) );
    opDB.register_new_operator<GradZVolSVol>( new GradZVolSVol(Gzvsva) );
  }

  //--------------------------------------------------------
  // interpolant scalar volume to scalar surface (diffusivities)
  //--------------------------------------------------------
  {
    InterpSVolSSurfX::Assembler Rsvssxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolSSurfY::Assembler Rsvssya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolSSurfZ::Assembler Rsvssza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    opDB.register_new_operator<InterpSVolSSurfX>( new InterpSVolSSurfX(Rsvssxa) );
    opDB.register_new_operator<InterpSVolSSurfY>( new InterpSVolSSurfY(Rsvssya) );
    opDB.register_new_operator<InterpSVolSSurfZ>( new InterpSVolSSurfZ(Rsvssza) );
  }

  //--------------------------------------------------------
  // interpolants - scalar volume to staggered surfaces (viscosity, dilatation)
  //--------------------------------------------------------
  {
    InterpSVolXSurfX::Assembler Rsvxsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolXSurfY::Assembler Rsvxsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolXSurfZ::Assembler Rsvxsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<InterpSVolXSurfX>( new InterpSVolXSurfX(Rsvxsxa) );  //opDB.retrieve_operator<InterpSVolXSurfX>()->write_matlab("Rsvxsx");
    opDB.register_new_operator<InterpSVolXSurfY>( new InterpSVolXSurfY(Rsvxsya) );  //opDB.retrieve_operator<InterpSVolXSurfY>()->write_matlab("Rsvxsy");
    opDB.register_new_operator<InterpSVolXSurfZ>( new InterpSVolXSurfZ(Rsvxsza) );  //opDB.retrieve_operator<InterpSVolXSurfZ>()->write_matlab("Rsvxsz");

    InterpSVolYSurfX::Assembler Rsvysxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolYSurfY::Assembler Rsvysya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolYSurfZ::Assembler Rsvysza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<InterpSVolYSurfX>( new InterpSVolYSurfX(Rsvysxa) );  //opDB.retrieve_operator<InterpSVolYSurfX>()->write_matlab("Rsvysx");
    opDB.register_new_operator<InterpSVolYSurfY>( new InterpSVolYSurfY(Rsvysya) );  //opDB.retrieve_operator<InterpSVolYSurfY>()->write_matlab("Rsvysy");
    opDB.register_new_operator<InterpSVolYSurfZ>( new InterpSVolYSurfZ(Rsvysza) );  //opDB.retrieve_operator<InterpSVolYSurfZ>()->write_matlab("Rsvysz");

    InterpSVolZSurfX::Assembler Rsvzsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolZSurfY::Assembler Rsvzsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolZSurfZ::Assembler Rsvzsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    opDB.register_new_operator<InterpSVolZSurfX>( new InterpSVolZSurfX(Rsvzsxa) );  //opDB.retrieve_operator<InterpSVolZSurfX>()->write_matlab("Rsvzsx");
    opDB.register_new_operator<InterpSVolZSurfY>( new InterpSVolZSurfY(Rsvzsya) );  //opDB.retrieve_operator<InterpSVolZSurfY>()->write_matlab("Rsvzsy");
    opDB.register_new_operator<InterpSVolZSurfZ>( new InterpSVolZSurfZ(Rsvzsza) );  //opDB.retrieve_operator<InterpSVolZSurfZ>()->write_matlab("Rsvzsz");
  }

  //--------------------------------------------------------
  // interpolants - scalar volume to staggered volume (density)
  //--------------------------------------------------------
  {
    InterpSVolXVol::Assembler Rsvxva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolYVol::Assembler Rsvyva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSVolZVol::Assembler Rsvzva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    opDB.register_new_operator<InterpSVolXVol>( new InterpSVolXVol(Rsvxva) );  //opDB.retrieve_operator<InterpSVolXVol>()->write_matlab("Rsvxv");
    opDB.register_new_operator<InterpSVolYVol>( new InterpSVolYVol(Rsvyva) );  //opDB.retrieve_operator<InterpSVolYVol>()->write_matlab("Rsvyv");
    opDB.register_new_operator<InterpSVolZVol>( new InterpSVolZVol(Rsvzva) );  //opDB.retrieve_operator<InterpSVolZVol>()->write_matlab("Rsvzv");
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

    opDB.register_new_operator<InterpXVolYSurfX>( new InterpXVolYSurfX(Rxvysxa) );  //opDB.retrieve_operator<InterpXVolYSurfX>()->write_matlab("Rxvysx");
    opDB.register_new_operator<InterpXVolZSurfX>( new InterpXVolZSurfX(Rxvzsxa) );  //opDB.retrieve_operator<InterpXVolZSurfX>()->write_matlab("Rxvzsx");

    opDB.register_new_operator<InterpYVolXSurfY>( new InterpYVolXSurfY(Ryvxsya) );  //opDB.retrieve_operator<InterpYVolXSurfY>()->write_matlab("Ryvxsy");
    opDB.register_new_operator<InterpYVolZSurfY>( new InterpYVolZSurfY(Ryvzsya) );

    opDB.register_new_operator<InterpZVolXSurfZ>( new InterpZVolXSurfZ(Rzvxsza) );
    opDB.register_new_operator<InterpZVolYSurfZ>( new InterpZVolYSurfZ(Rzvysza) );

    // interpolants - volume to surface for staggered cells.
    {
      InterpXVolXSurfX::Assembler Rxvxsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpXVolXSurfY::Assembler Rxvxsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpXVolXSurfZ::Assembler Rxvxsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

      opDB.register_new_operator<InterpXVolXSurfX>( new InterpXVolXSurfX(Rxvxsxa) );  //opDB.retrieve_operator<InterpXVolXSurfX>()->write_matlab("Rxvxsx");
      opDB.register_new_operator<InterpXVolXSurfY>( new InterpXVolXSurfY(Rxvxsya) );  //opDB.retrieve_operator<InterpXVolXSurfY>()->write_matlab("Rxvxsy");
      opDB.register_new_operator<InterpXVolXSurfZ>( new InterpXVolXSurfZ(Rxvxsza) );  //opDB.retrieve_operator<InterpXVolXSurfZ>()->write_matlab("Rxvxsz");


      InterpYVolYSurfX::Assembler Ryvysxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpYVolYSurfY::Assembler Ryvysya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpYVolYSurfZ::Assembler Ryvysza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

      opDB.register_new_operator<InterpYVolYSurfX>( new InterpYVolYSurfX(Ryvysxa) );
      opDB.register_new_operator<InterpYVolYSurfY>( new InterpYVolYSurfY(Ryvysya) );
      opDB.register_new_operator<InterpYVolYSurfZ>( new InterpYVolYSurfZ(Ryvysza) );


      InterpZVolZSurfX::Assembler Rzvzsxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpZVolZSurfY::Assembler Rzvzsya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
      InterpZVolZSurfZ::Assembler Rzvzsza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

      opDB.register_new_operator<InterpZVolZSurfX>( new InterpZVolZSurfX(Rzvzsxa) );  //opDB.retrieve_operator<InterpZVolZSurfX>()->write_matlab("Rzvzsx");
      opDB.register_new_operator<InterpZVolZSurfY>( new InterpZVolZSurfY(Rzvzsya) );  //
      opDB.register_new_operator<InterpZVolZSurfZ>( new InterpZVolZSurfZ(Rzvzsza) );  //opDB.retrieve_operator<InterpZVolZSurfZ>()->write_matlab("Rzvzsz");
    }
  }

  //--------------------------------------------------------
  // interpolants - staggard volume to scalar surface  (advecting velocities)
  //--------------------------------------------------------
  {
    InterpXVolSSurfX::Assembler Rxvssxa( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpYVolSSurfY::Assembler Ryvssya( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpZVolSSurfZ::Assembler Rzvssza( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    opDB.register_new_operator<InterpXVolSSurfX>( new InterpXVolSSurfX(Rxvssxa) );
    opDB.register_new_operator<InterpYVolSSurfY>( new InterpYVolSSurfY(Ryvssya) );
    opDB.register_new_operator<InterpZVolSSurfZ>( new InterpZVolSSurfZ(Rzvssza) );
  }

  //--------------------------------------------------------
  // scalar surface to staggered volumes (pressure gradients)
  //--------------------------------------------------------
  {
    InterpSSurfXXVol::Assembler Rsxxva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSSurfYYVol::Assembler Rsyyva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );
    InterpSSurfZZVol::Assembler Rszzva( dim, bcPlus[0], bcPlus[1], bcPlus[2] );

    opDB.register_new_operator<InterpSSurfXXVol>( new InterpSSurfXXVol(Rsxxva) );
    opDB.register_new_operator<InterpSSurfYYVol>( new InterpSSurfYYVol(Rsyyva) );
    opDB.register_new_operator<InterpSSurfZZVol>( new InterpSSurfZZVol(Rszzva) );
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
    
    opDB.register_new_operator<ScratchSVol>( new ScratchSVol(Ssvxa) );
    opDB.register_new_operator<ScratchSVol>( new ScratchSVol(Ssvya) );
    opDB.register_new_operator<ScratchSVol>( new ScratchSVol(Ssvza) );

    opDB.register_new_operator<ScratchXVol>( new ScratchXVol(Sxvxa) );
    opDB.register_new_operator<ScratchXVol>( new ScratchXVol(Sxvya) );
    opDB.register_new_operator<ScratchXVol>( new ScratchXVol(Sxvza) );

    opDB.register_new_operator<ScratchYVol>( new ScratchYVol(Syvxa) );
    opDB.register_new_operator<ScratchYVol>( new ScratchYVol(Syvya) );
    opDB.register_new_operator<ScratchYVol>( new ScratchYVol(Syvza) );

    opDB.register_new_operator<ScratchZVol>( new ScratchZVol(Szvxa) );
    opDB.register_new_operator<ScratchZVol>( new ScratchZVol(Szvya) );
    opDB.register_new_operator<ScratchZVol>( new ScratchZVol(Szvza) );


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

    opDB.register_new_operator<ScratchSVolSSurfX>( new ScratchSVolSSurfX(Ssvssxa) );
    opDB.register_new_operator<ScratchSVolSSurfY>( new ScratchSVolSSurfY(Ssvssya) );
    opDB.register_new_operator<ScratchSVolSSurfZ>( new ScratchSVolSSurfZ(Ssvssza) );

    opDB.register_new_operator<ScratchXVolXSurfX>( new ScratchXVolXSurfX(Sxvxsxa) );
    opDB.register_new_operator<ScratchXVolXSurfY>( new ScratchXVolXSurfY(Sxvxsya) );
    opDB.register_new_operator<ScratchXVolXSurfZ>( new ScratchXVolXSurfZ(Sxvxsza) );

    opDB.register_new_operator<ScratchYVolYSurfX>( new ScratchYVolYSurfX(Syvysxa) );
    opDB.register_new_operator<ScratchYVolYSurfY>( new ScratchYVolYSurfY(Syvysya) );
    opDB.register_new_operator<ScratchYVolYSurfZ>( new ScratchYVolYSurfZ(Syvysza) );

    opDB.register_new_operator<ScratchZVolZSurfX>( new ScratchZVolZSurfX(Szvzsxa) );
    opDB.register_new_operator<ScratchZVolZSurfY>( new ScratchZVolZSurfY(Szvzsya) );
    opDB.register_new_operator<ScratchZVolZSurfZ>( new ScratchZVolZSurfZ(Szvzsza) );

    // bc operators
    ScratchSSurfX::Assembler Sssxa( dim, XDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSSurfY::Assembler Sssya( dim, YDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );
    ScratchSSurfZ::Assembler Sssza( dim, ZDIR::value, bcPlus[0], bcPlus[1], bcPlus[2] );

    opDB.register_new_operator<ScratchSSurfX>( new ScratchSSurfX(Sssxa) );
    opDB.register_new_operator<ScratchSSurfY>( new ScratchSSurfY(Sssya) );
    opDB.register_new_operator<ScratchSSurfZ>( new ScratchSSurfZ(Sssza) );
  }
}
