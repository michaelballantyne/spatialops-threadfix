#ifndef FVTEST_GRID_H
#define FVTEST_GRID_H

#include <spatialops/structured/FVStaggered.h>
#include <spatialops/OperatorDatabase.h> //jcs shove into FVStaggered.h?

namespace SpatialOps{
namespace structured{


  /**
   *  @class Grid
   *  @author James C. Sutherland
   *  @date   August, 2007
   *
   *  @brief Provides coordinate information for a uniform, 3-D
   *  staggered mesh configuration.  Allows you to obtain the
   *  coordinates for any field location on a staggered uniform mesh.
   */
class Grid
{
public:
  Grid( const std::vector<int>& dim,
        const std::vector<double>& spacing,
        const std::vector<bool>& bcPlusFlag,
        const OperatorDatabase& opDB );
  ~Grid();

  const std::vector<int>& extent() const{return dim_;}

  void write() const;


  // scalar-cell coordinates

  const SVolField& xcoord_svol() const{return svx_;}
  const SVolField& ycoord_svol() const{return svy_;}
  const SVolField& zcoord_svol() const{return svz_;}

  const SSurfXField& xcoord_sxsurf() const{return ssxx_;}
  const SSurfXField& ycoord_sxsurf() const{return ssxy_;}
  const SSurfXField& zcoord_sxsurf() const{return ssxz_;}

  const SSurfYField& xcoord_sysurf() const{return ssyx_;}
  const SSurfYField& ycoord_sysurf() const{return ssyy_;}
  const SSurfYField& zcoord_sysurf() const{return ssyz_;}

  const SSurfZField& xcoord_szsurf() const{return sszx_;}
  const SSurfZField& ycoord_szsurf() const{return sszy_;}
  const SSurfZField& zcoord_szsurf() const{return sszz_;}


  // x-cell coordinates

  const XVolField& xcoord_xvol() const{return xvx_;}
  const XVolField& ycoord_xvol() const{return xvy_;}
  const XVolField& zcoord_xvol() const{return xvz_;}

  const XSurfXField& xcoord_xxsurf() const{return xsxx_;}
  const XSurfXField& ycoord_xxsurf() const{return xsxy_;}
  const XSurfXField& zcoord_xxsurf() const{return xsxz_;}

  const XSurfYField& xcoord_xysurf() const{return xsyx_;}
  const XSurfYField& ycoord_xysurf() const{return xsyy_;}
  const XSurfYField& zcoord_xysurf() const{return xsyz_;}

  const XSurfZField& xcoord_xzsurf() const{return xszx_;}
  const XSurfZField& ycoord_xzsurf() const{return xszy_;}
  const XSurfZField& zcoord_xzsurf() const{return xszz_;}


  // y-cell coordinates

  const YVolField& xcoord_yvol() const{return yvx_;}
  const YVolField& ycoord_yvol() const{return yvy_;}
  const YVolField& zcoord_yvol() const{return yvz_;}

  const YSurfXField& xcoord_yxsurf() const{return ysxx_;}
  const YSurfXField& ycoord_yxsurf() const{return ysxy_;}
  const YSurfXField& zcoord_yxsurf() const{return ysxz_;}

  const YSurfYField& xcoord_yysurf() const{return ysyx_;}
  const YSurfYField& ycoord_yysurf() const{return ysyy_;}
  const YSurfYField& zcoord_yysurf() const{return ysyz_;}

  const YSurfZField& xcoord_yzsurf() const{return yszx_;}
  const YSurfZField& ycoord_yzsurf() const{return yszy_;}
  const YSurfZField& zcoord_yzsurf() const{return yszz_;}


  // z-cell coordinates

  const ZVolField& xcoord_zvol() const{return zvx_;}
  const ZVolField& ycoord_zvol() const{return zvy_;}
  const ZVolField& zcoord_zvol() const{return zvz_;}

  const ZSurfXField& xcoord_zxsurf() const{return zsxx_;}
  const ZSurfXField& ycoord_zxsurf() const{return zsxy_;}
  const ZSurfXField& zcoord_zxsurf() const{return zsxz_;}

  const ZSurfYField& xcoord_zysurf() const{return zsyx_;}
  const ZSurfYField& ycoord_zysurf() const{return zsyy_;}
  const ZSurfYField& zcoord_zysurf() const{return zsyz_;}

  const ZSurfZField& xcoord_zzsurf() const{return zszx_;}
  const ZSurfZField& ycoord_zzsurf() const{return zszy_;}
  const ZSurfZField& zcoord_zzsurf() const{return zszz_;}

private:

  const std::vector<int> dim_;

  SVolField  svx_, svy_, svz_;
  SSurfXField ssxx_, ssxy_, ssxz_;
  SSurfYField ssyx_, ssyy_, ssyz_;
  SSurfZField sszx_, sszy_, sszz_;

  XVolField  xvx_, xvy_, xvz_;
  XSurfXField xsxx_, xsxy_, xsxz_;
  XSurfYField xsyx_, xsyy_, xsyz_;
  XSurfZField xszx_, xszy_, xszz_;

  YVolField  yvx_, yvy_, yvz_;
  YSurfXField ysxx_, ysxy_, ysxz_;
  YSurfYField ysyx_, ysyy_, ysyz_;
  YSurfZField yszx_, yszy_, yszz_;

  ZVolField  zvx_, zvy_, zvz_;
  ZSurfXField zsxx_, zsxy_, zsxz_;
  ZSurfYField zsyx_, zsyy_, zsyz_;
  ZSurfZField zszx_, zszy_, zszz_;

};

//====================================================================

Grid::Grid( const std::vector<int>& dim,
            const std::vector<double>& spacing,
            const std::vector<bool>& bcPlusFlag,
            const OperatorDatabase& opDB )
  : dim_( dim ),

    svx_( get_n_tot<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    svy_( get_n_tot<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    svz_( get_n_tot<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

    ssxx_( get_n_tot<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ssxy_( get_n_tot<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ssxz_( get_n_tot<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ssyx_( get_n_tot<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ssyy_( get_n_tot<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ssyz_( get_n_tot<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    sszx_( get_n_tot<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    sszy_( get_n_tot<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    sszz_( get_n_tot<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),


    xvx_( get_n_tot<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xvy_( get_n_tot<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xvz_( get_n_tot<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),


    xsxx_( get_n_tot<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xsxy_( get_n_tot<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xsxz_( get_n_tot<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xsyx_( get_n_tot<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xsyy_( get_n_tot<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xsyz_( get_n_tot<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xszx_( get_n_tot<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xszy_( get_n_tot<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    xszz_( get_n_tot<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

    yvx_( get_n_tot<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    yvy_( get_n_tot<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    yvz_( get_n_tot<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

    ysxx_( get_n_tot<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ysxy_( get_n_tot<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ysxz_( get_n_tot<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ysyx_( get_n_tot<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ysyy_( get_n_tot<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    ysyz_( get_n_tot<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    yszx_( get_n_tot<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    yszy_( get_n_tot<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    yszz_( get_n_tot<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),


    zvx_( get_n_tot<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zvy_( get_n_tot<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zvz_( get_n_tot<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

    zsxx_( get_n_tot<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zsxy_( get_n_tot<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zsxz_( get_n_tot<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zsyx_( get_n_tot<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zsyy_( get_n_tot<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zsyz_( get_n_tot<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zszx_( get_n_tot<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zszy_( get_n_tot<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
    zszz_( get_n_tot<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), get_ghost_set<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL )
{
  using namespace SpatialOps;
  using namespace structured;

  {
    SVolField::iterator isvx=svx_.begin();
    SVolField::iterator isvy=svy_.begin();
    SVolField::iterator isvz=svz_.begin();

    const int ihi = get_nx<SVolField>(dim,bcPlusFlag[0]);
    const int jhi = get_ny<SVolField>(dim,bcPlusFlag[1]);
    const int khi = get_nz<SVolField>(dim,bcPlusFlag[2]);

    const int ngxm = dim[0]>1 ? SVolField::Ghost::NGHOST : 0;
    const int ngym = dim[1]>1 ? SVolField::Ghost::NGHOST : 0;
    const int ngzm = dim[2]>1 ? SVolField::Ghost::NGHOST : 0;

    // set the mesh for the scalar cell
    for( int k=0; k<khi; ++k ){
      const double z = spacing[2]*(double(k)+0.5-ngzm);
      for( int j=0; j<jhi; ++j ){
        const double y = spacing[1]*(double(j)+0.5-ngym);
        for( int i=0; i<ihi; ++i ){
          const double x = spacing[0]*(double(i)+0.5-ngxm);
          *isvx = x;
          *isvy = y;
          *isvz = z;
          ++isvx; ++isvy; ++isvz;
        }
      }
    }

  }


  if( dim[0]>1 ){
    const int ngxm = dim[0]>1 ? XVolField::Ghost::NGHOST : 0;
    const int ngym = dim[1]>1 ? XVolField::Ghost::NGHOST : 0;
    const int ngzm = dim[2]>1 ? XVolField::Ghost::NGHOST : 0;

    XVolField::iterator ixvx=xvx_.begin();
    XVolField::iterator ixvy=xvy_.begin();
    XVolField::iterator ixvz=xvz_.begin();
    const int ihi = get_nx<XVolField>(dim,bcPlusFlag[0]);
    const int jhi = get_ny<XVolField>(dim,bcPlusFlag[1]);
    const int khi = get_nz<XVolField>(dim,bcPlusFlag[2]);

    for( int k=0; k<khi; ++k ){
      const double z = spacing[2]*(double(k)+0.5-ngzm);
      for( int j=0; j<jhi; ++j ){
        const double y = spacing[1]*(double(j)+0.5-ngym);
        for( int i=0; i<ihi; ++i ){
          const double x = spacing[0]*(i-ngxm); // offset in x.
          *ixvx = x;
          *ixvy = y;
          *ixvz = z;
          ++ixvx; ++ixvy; ++ixvz;
        }
      }
    }

  }

  if( dim[1]>1 ){
    const int ngxm = dim[0]>1 ? YVolField::Ghost::NGHOST : 0;
    const int ngym = dim[1]>1 ? YVolField::Ghost::NGHOST : 0;
    const int ngzm = dim[2]>1 ? YVolField::Ghost::NGHOST : 0;

    YVolField::iterator iyvx=yvx_.begin();
    YVolField::iterator iyvy=yvy_.begin();
    YVolField::iterator iyvz=yvz_.begin();
    const int ihi = get_nx<YVolField>(dim,bcPlusFlag[0]);
    const int jhi = get_ny<YVolField>(dim,bcPlusFlag[1]);
    const int khi = get_nz<YVolField>(dim,bcPlusFlag[2]);

    for( int k=0; k<khi; ++k ){
      const double z = spacing[2]*(double(k)+0.5-ngzm);
      for( int j=0; j<jhi; ++j ){
        const double y = spacing[1]*(j-ngym); // offset in y.
        for( int i=0; i<ihi; ++i ){
          const double x = spacing[0]*(double(i)+0.5-ngxm);
          *iyvx = x;
          *iyvy = y;
          *iyvz = z;
          ++iyvx; ++iyvy; ++iyvz;
        }
      }
    }

  }

  if( dim[2]>1 ){
    const int ngxm = dim[0]>1 ? ZVolField::Ghost::NGHOST : 0;
    const int ngym = dim[1]>1 ? ZVolField::Ghost::NGHOST : 0;
    const int ngzm = dim[2]>1 ? ZVolField::Ghost::NGHOST : 0;

    ZVolField::iterator izvx=zvx_.begin();
    ZVolField::iterator izvy=zvy_.begin();
    ZVolField::iterator izvz=zvz_.begin();
    const int ihi = get_nx<ZVolField>(dim,bcPlusFlag[0]);
    const int jhi = get_ny<ZVolField>(dim,bcPlusFlag[1]);
    const int khi = get_nz<ZVolField>(dim,bcPlusFlag[2]);

    for( int k=0; k<khi; ++k ){
      const double z = spacing[2]*(k-ngzm); // offset in z
      for( int j=0; j<jhi; ++j ){
        const double y = spacing[1]*(double(j)+0.5-ngym);
        for( int i=0; i<ihi; ++i ){
          const double x = spacing[0]*(double(i)+0.5-ngxm);
          *izvx = x;
          *izvy = y;
          *izvz = z;
          ++izvx; ++izvy; ++izvz;
        }
      }
    }

  }


  // scalar cell grid coordinates
  {
    const InterpSVolSSurfX* const Rsvssx = opDB.retrieve_operator<InterpSVolSSurfX>();
    const InterpSVolSSurfY* const Rsvssy = opDB.retrieve_operator<InterpSVolSSurfY>();
    const InterpSVolSSurfZ* const Rsvssz = opDB.retrieve_operator<InterpSVolSSurfZ>();
    Rsvssx->apply_to_field( svx_, ssxx_ );
    Rsvssx->apply_to_field( svy_, ssxy_ );
    Rsvssx->apply_to_field( svz_, ssxz_ );

    Rsvssy->apply_to_field( svx_, ssyx_ );
    Rsvssy->apply_to_field( svy_, ssyy_ );
    Rsvssy->apply_to_field( svz_, ssyz_ );

    Rsvssz->apply_to_field( svx_, sszx_ );
    Rsvssz->apply_to_field( svy_, sszy_ );
    Rsvssz->apply_to_field( svz_, sszz_ );

  }

  // x-cell grid coordinates
  {
    const InterpXVolXSurfX* const Rxvxsx = opDB.retrieve_operator<InterpXVolXSurfX>();
    const InterpXVolXSurfY* const Rxvxsy = opDB.retrieve_operator<InterpXVolXSurfY>();
    const InterpXVolXSurfZ* const Rxvxsz = opDB.retrieve_operator<InterpXVolXSurfZ>();
    Rxvxsx->apply_to_field( xvx_, xsxx_ );
    Rxvxsx->apply_to_field( xvy_, xsxy_ );
    Rxvxsx->apply_to_field( xvz_, xsxz_ );

    Rxvxsy->apply_to_field( xvx_, xsyx_ );
    Rxvxsy->apply_to_field( xvy_, xsyy_ );
    Rxvxsy->apply_to_field( xvz_, xsyz_ );

    Rxvxsz->apply_to_field( xvx_, xszx_ );
    Rxvxsz->apply_to_field( xvy_, xszy_ );
    Rxvxsz->apply_to_field( xvz_, xszz_ );
  }

  // y-cell grid coordinates
  {
    const InterpYVolYSurfX* const Ryvysx = opDB.retrieve_operator<InterpYVolYSurfX>();
    const InterpYVolYSurfY* const Ryvysy = opDB.retrieve_operator<InterpYVolYSurfY>();
    const InterpYVolYSurfZ* const Ryvysz = opDB.retrieve_operator<InterpYVolYSurfZ>();
    Ryvysx->apply_to_field( yvx_, ysxx_ );
    Ryvysx->apply_to_field( yvy_, ysxy_ );
    Ryvysx->apply_to_field( yvz_, ysxz_ );

    Ryvysy->apply_to_field( yvx_, ysyx_ );
    Ryvysy->apply_to_field( yvy_, ysyy_ );
    Ryvysy->apply_to_field( yvz_, ysyz_ );

    Ryvysz->apply_to_field( yvx_, yszx_ );
    Ryvysz->apply_to_field( yvy_, yszy_ );
    Ryvysz->apply_to_field( yvz_, yszz_ );
  }

  // z-cell grid coordinates
  {
    const InterpZVolZSurfX* const Rzvzsx = opDB.retrieve_operator<InterpZVolZSurfX>();
    const InterpZVolZSurfY* const Rzvzsy = opDB.retrieve_operator<InterpZVolZSurfY>();
    const InterpZVolZSurfZ* const Rzvzsz = opDB.retrieve_operator<InterpZVolZSurfZ>();
    Rzvzsx->apply_to_field( zvx_, zsxx_ );
    Rzvzsx->apply_to_field( zvy_, zsxy_ );
    Rzvzsx->apply_to_field( zvz_, zsxz_ );

    Rzvzsy->apply_to_field( zvx_, zsyx_ );
    Rzvzsy->apply_to_field( zvy_, zsyy_ );
    Rzvzsy->apply_to_field( zvz_, zsyz_ );

    Rzvzsz->apply_to_field( zvx_, zszx_ );
    Rzvzsz->apply_to_field( zvy_, zszy_ );
    Rzvzsz->apply_to_field( zvz_, zszz_ );
  }


}
//--------------------------------------------------------------------
Grid::~Grid()
{
}

void
Grid::write() const
{

  if( dim_[0]>1 ){
    svx_.write_matlab("svx");
    ssxx_.write_matlab("ssxx");

    xvx_.write_matlab("xvx");
    xsxx_.write_matlab("xsxx");
    xsyx_.write_matlab("xsyx");
    xszx_.write_matlab("xszx");

    if( dim_[1]>1 ){
      ssxy_.write_matlab("ssxy");
      xvy_.write_matlab("xvy");
      xsxy_.write_matlab("xsxy");
      xsyy_.write_matlab("xsyy");
      xszy_.write_matlab("xszy");
    }
    if( dim_[2]>1 ){
      ssxz_.write_matlab("ssxz");
      xvz_.write_matlab("xvz");
      xsxz_.write_matlab("xsxz");
      xsyz_.write_matlab("xsyz");
      xszz_.write_matlab("xszz");
    }
  }

  if( dim_[1]>1 ){
    if( dim_[0]>1 ){
      ssyx_.write_matlab("ssyx");

      yvx_.write_matlab("yvx");
      ysxx_.write_matlab("ysxx");
      ysyx_.write_matlab("ysyx");
      yszx_.write_matlab("yszx");
    }

    svy_.write_matlab("svy");
    ssyy_.write_matlab("ssyy");

    yvy_.write_matlab("yvy");
    ysxy_.write_matlab("ysxy");
    ysyy_.write_matlab("ysyy");
    yszy_.write_matlab("yszy");

    if( dim_[2]>1 ){
      ssyz_.write_matlab("ssyz");
      yvz_.write_matlab("yvz");
      ysxz_.write_matlab("ysxz");
      ysyz_.write_matlab("ysyz");
      yszz_.write_matlab("yszz");
    }
  }

  if( dim_[2]>1 ){
    if( dim_[0]>1 ){
      sszx_.write_matlab("sszx");
      zvx_.write_matlab("zvx");
      zsxx_.write_matlab("zsxx");
      zsyx_.write_matlab("zsyx");
      zszx_.write_matlab("zszx");
    }

    if( dim_[1]>1 ){
      sszy_.write_matlab("sszy");
      zvy_.write_matlab("zvy");
      zsxy_.write_matlab("zsxy");
      zsyy_.write_matlab("zsyy");
      zszy_.write_matlab("zszy");
    }

    svz_.write_matlab("svz");
    sszz_.write_matlab("sszz");

    yvz_.write_matlab("yvz");
    ysxz_.write_matlab("ysxz");
    ysyz_.write_matlab("ysyz");
    yszz_.write_matlab("yszz");
  }

}

}
}

#endif
