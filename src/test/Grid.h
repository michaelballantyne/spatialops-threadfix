#ifndef FVTEST_GRID_H
#define FVTEST_GRID_H

#include <FVStaggered.h>

namespace SpatialOps{
namespace FVStaggered{


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
	const std::vector<double>& spacing );
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
	    const std::vector<double>& spacing )
  : dim_( dim ),

    svx_( get_n_tot<SVolField >(dim), get_ghost_set<SVolField >(dim), NULL ),
    svy_( get_n_tot<SVolField >(dim), get_ghost_set<SVolField >(dim), NULL ),
    svz_( get_n_tot<SVolField >(dim), get_ghost_set<SVolField >(dim), NULL ),

    ssxx_( get_n_tot<SSurfXField>(dim), get_ghost_set<SSurfXField>(dim), NULL ),
    ssxy_( get_n_tot<SSurfXField>(dim), get_ghost_set<SSurfXField>(dim), NULL ),
    ssxz_( get_n_tot<SSurfXField>(dim), get_ghost_set<SSurfXField>(dim), NULL ),
    ssyx_( get_n_tot<SSurfYField>(dim), get_ghost_set<SSurfYField>(dim), NULL ),
    ssyy_( get_n_tot<SSurfYField>(dim), get_ghost_set<SSurfYField>(dim), NULL ),
    ssyz_( get_n_tot<SSurfYField>(dim), get_ghost_set<SSurfYField>(dim), NULL ),
    sszx_( get_n_tot<SSurfZField>(dim), get_ghost_set<SSurfZField>(dim), NULL ),
    sszy_( get_n_tot<SSurfZField>(dim), get_ghost_set<SSurfZField>(dim), NULL ),
    sszz_( get_n_tot<SSurfZField>(dim), get_ghost_set<SSurfZField>(dim), NULL ),


    xvx_( get_n_tot<XVolField >(dim), get_ghost_set<XVolField >(dim), NULL ),
    xvy_( get_n_tot<XVolField >(dim), get_ghost_set<XVolField >(dim), NULL ),
    xvz_( get_n_tot<XVolField >(dim), get_ghost_set<XVolField >(dim), NULL ),


    xsxx_( get_n_tot<XSurfXField>(dim), get_ghost_set<XSurfXField>(dim), NULL ),
    xsxy_( get_n_tot<XSurfXField>(dim), get_ghost_set<XSurfXField>(dim), NULL ),
    xsxz_( get_n_tot<XSurfXField>(dim), get_ghost_set<XSurfXField>(dim), NULL ),
    xsyx_( get_n_tot<XSurfYField>(dim), get_ghost_set<XSurfYField>(dim), NULL ),
    xsyy_( get_n_tot<XSurfYField>(dim), get_ghost_set<XSurfYField>(dim), NULL ),
    xsyz_( get_n_tot<XSurfYField>(dim), get_ghost_set<XSurfYField>(dim), NULL ),
    xszx_( get_n_tot<XSurfZField>(dim), get_ghost_set<XSurfZField>(dim), NULL ),
    xszy_( get_n_tot<XSurfZField>(dim), get_ghost_set<XSurfZField>(dim), NULL ),
    xszz_( get_n_tot<XSurfZField>(dim), get_ghost_set<XSurfZField>(dim), NULL ),

    yvx_( get_n_tot<YVolField >(dim), get_ghost_set<YVolField >(dim), NULL ),
    yvy_( get_n_tot<YVolField >(dim), get_ghost_set<YVolField >(dim), NULL ),
    yvz_( get_n_tot<YVolField >(dim), get_ghost_set<YVolField >(dim), NULL ),

    ysxx_( get_n_tot<YSurfXField>(dim), get_ghost_set<YSurfXField>(dim), NULL ),
    ysxy_( get_n_tot<YSurfXField>(dim), get_ghost_set<YSurfXField>(dim), NULL ),
    ysxz_( get_n_tot<YSurfXField>(dim), get_ghost_set<YSurfXField>(dim), NULL ),
    ysyx_( get_n_tot<YSurfYField>(dim), get_ghost_set<YSurfYField>(dim), NULL ),
    ysyy_( get_n_tot<YSurfYField>(dim), get_ghost_set<YSurfYField>(dim), NULL ),
    ysyz_( get_n_tot<YSurfYField>(dim), get_ghost_set<YSurfYField>(dim), NULL ),
    yszx_( get_n_tot<YSurfZField>(dim), get_ghost_set<YSurfZField>(dim), NULL ),
    yszy_( get_n_tot<YSurfZField>(dim), get_ghost_set<YSurfZField>(dim), NULL ),
    yszz_( get_n_tot<YSurfZField>(dim), get_ghost_set<YSurfZField>(dim), NULL ),


    zvx_( get_n_tot<ZVolField >(dim), get_ghost_set<ZVolField >(dim), NULL ),
    zvy_( get_n_tot<ZVolField >(dim), get_ghost_set<ZVolField >(dim), NULL ),
    zvz_( get_n_tot<ZVolField >(dim), get_ghost_set<ZVolField >(dim), NULL ),

    zsxx_( get_n_tot<ZSurfXField>(dim), get_ghost_set<ZSurfXField>(dim), NULL ),
    zsxy_( get_n_tot<ZSurfXField>(dim), get_ghost_set<ZSurfXField>(dim), NULL ),
    zsxz_( get_n_tot<ZSurfXField>(dim), get_ghost_set<ZSurfXField>(dim), NULL ),
    zsyx_( get_n_tot<ZSurfYField>(dim), get_ghost_set<ZSurfYField>(dim), NULL ),
    zsyy_( get_n_tot<ZSurfYField>(dim), get_ghost_set<ZSurfYField>(dim), NULL ),
    zsyz_( get_n_tot<ZSurfYField>(dim), get_ghost_set<ZSurfYField>(dim), NULL ),
    zszx_( get_n_tot<ZSurfZField>(dim), get_ghost_set<ZSurfZField>(dim), NULL ),
    zszy_( get_n_tot<ZSurfZField>(dim), get_ghost_set<ZSurfZField>(dim), NULL ),
    zszz_( get_n_tot<ZSurfZField>(dim), get_ghost_set<ZSurfZField>(dim), NULL )
{
  using namespace SpatialOps;
  using namespace FVStaggered;

  {
    SVolField::iterator isvx=svx_.begin();
    SVolField::iterator isvy=svy_.begin();
    SVolField::iterator isvz=svz_.begin();

    const int ihi = get_nx<SVolField>(dim[0]);
    const int jhi = get_ny<SVolField>(dim[1]);
    const int khi = get_nz<SVolField>(dim[2]);

    const int ngxm = dim[0]>1 ? SVolField::Ghost::NM : 0;
    const int ngym = dim[1]>1 ? SVolField::Ghost::NM : 0;
    const int ngzm = dim[2]>1 ? SVolField::Ghost::NM : 0;

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
    const int ngxm = dim[0]>1 ? XVolField::Ghost::NM : 0;
    const int ngym = dim[1]>1 ? XVolField::Ghost::NM : 0;
    const int ngzm = dim[2]>1 ? XVolField::Ghost::NM : 0;

    XVolField::iterator ixvx=xvx_.begin();
    XVolField::iterator ixvy=xvy_.begin();
    XVolField::iterator ixvz=xvz_.begin();
    const int ihi = get_nx<XVolField>(dim[0]);
    const int jhi = get_ny<XVolField>(dim[1]);
    const int khi = get_nz<XVolField>(dim[2]);

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
    const int ngxm = dim[0]>1 ? YVolField::Ghost::NM : 0;
    const int ngym = dim[1]>1 ? YVolField::Ghost::NM : 0;
    const int ngzm = dim[2]>1 ? YVolField::Ghost::NM : 0;

    YVolField::iterator iyvx=yvx_.begin();
    YVolField::iterator iyvy=yvy_.begin();
    YVolField::iterator iyvz=yvz_.begin();
    const int ihi = get_nx<YVolField>(dim[0]);
    const int jhi = get_ny<YVolField>(dim[1]);
    const int khi = get_nz<YVolField>(dim[2]);

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
    const int ngxm = dim[0]>1 ? ZVolField::Ghost::NM : 0;
    const int ngym = dim[1]>1 ? ZVolField::Ghost::NM : 0;
    const int ngzm = dim[2]>1 ? ZVolField::Ghost::NM : 0;

    ZVolField::iterator izvx=zvx_.begin();
    ZVolField::iterator izvy=zvy_.begin();
    ZVolField::iterator izvz=zvz_.begin();
    const int ihi = get_nx<ZVolField>(dim[0]);
    const int jhi = get_ny<ZVolField>(dim[1]);
    const int khi = get_nz<ZVolField>(dim[2]);

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

    // interpolate to get face values.
//     const InterpSVolSSurf* const Rsvss = SpatialOpDatabase<InterpSVolSSurf>::self().retrieve_operator();
//     EpetraExt::RowMatrixToMatrixMarketFile( "Rsvss.mm", Rsvss->get_linalg_mat(), "", "" );

//     Rsvss->apply_to_field( svx_, ssx_ );
//     Rsvss->apply_to_field( svy_, ssy_ );
//     Rsvss->apply_to_field( svz_, ssz_ );

//     // assign surface values to vector components (for convenience)
//     ssxx_ = ssx_;
//     ssxy_ = ssy_;
//     ssxz_ = ssz_;
    
//     ssyx_ = ssx_;
//     ssyy_ = ssy_;
//     ssyz_ = ssz_;

//     sszx_ = ssx_;
//     sszy_ = ssy_;
//     sszz_ = ssz_;


    const InterpSVolSSurfX* const Rsvssx = SpatialOpDatabase<InterpSVolSSurfX>::self().retrieve_operator();
    const InterpSVolSSurfY* const Rsvssy = SpatialOpDatabase<InterpSVolSSurfY>::self().retrieve_operator();
    const InterpSVolSSurfZ* const Rsvssz = SpatialOpDatabase<InterpSVolSSurfZ>::self().retrieve_operator();
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
//     const InterpSVolXVol* const Rsvxv = SpatialOpDatabase<InterpSVolXVol>::self().retrieve_operator();
//     Rsvxv->apply_to_field( svx_, xvx_ );
//     Rsvxv->apply_to_field( svy_, xvy_ );
//     Rsvxv->apply_to_field( svz_, xvz_ );
  

    const InterpXVolXSurfX* const Rxvxsx = SpatialOpDatabase<InterpXVolXSurfX>::self().retrieve_operator();
    const InterpXVolXSurfY* const Rxvxsy = SpatialOpDatabase<InterpXVolXSurfY>::self().retrieve_operator();
    const InterpXVolXSurfZ* const Rxvxsz = SpatialOpDatabase<InterpXVolXSurfZ>::self().retrieve_operator();
    Rxvxsx->apply_to_field( xvx_, xsxx_ );
    Rxvxsx->apply_to_field( xvy_, xsxy_ );
    Rxvxsx->apply_to_field( xvz_, xsxz_ );

    Rxvxsy->apply_to_field( xvx_, xsyx_ );
    Rxvxsy->apply_to_field( xvy_, xsyy_ );
    Rxvxsy->apply_to_field( xvz_, xsyz_ );

    Rxvxsz->apply_to_field( xvx_, xszx_ );
    Rxvxsz->apply_to_field( xvy_, xszy_ );
    Rxvxsz->apply_to_field( xvz_, xszz_ );

//     const InterpSVolXSurf* const Rsvxs = SpatialOpDatabase<InterpSVolXSurf>::self().retrieve_operator();
//     Rsvxs->apply_to_field( svx_, xsx_ );
//     Rsvxs->apply_to_field( svy_, xsy_ );
//     Rsvxs->apply_to_field( svz_, xsz_ );

    /*
    // assign surface values to vector components (for convenience)
    xsxx_ = xsx_;    xsxy_ = xsy_;    xsxz_ = xsz_;
    xsyx_ = xsx_;    xsyy_ = xsy_;    xsyz_ = xsz_;
    xszx_ = xsx_;    xszy_ = xsy_;    xszz_ = xsz_;
    */
  }

  // y-cell grid coordinates
  {
//     const InterpSVolYVol* const Rsvyv = SpatialOpDatabase<InterpSVolYVol>::self().retrieve_operator();
//     Rsvyv->apply_to_field( svx_, yvx_ );
//     Rsvyv->apply_to_field( svy_, yvy_ );
//     Rsvyv->apply_to_field( svz_, yvz_ );


    const InterpYVolYSurfX* const Ryvysx = SpatialOpDatabase<InterpYVolYSurfX>::self().retrieve_operator();
    const InterpYVolYSurfY* const Ryvysy = SpatialOpDatabase<InterpYVolYSurfY>::self().retrieve_operator();
    const InterpYVolYSurfZ* const Ryvysz = SpatialOpDatabase<InterpYVolYSurfZ>::self().retrieve_operator();
    Ryvysx->apply_to_field( yvx_, ysxx_ );
    Ryvysx->apply_to_field( yvy_, ysxy_ );
    Ryvysx->apply_to_field( yvz_, ysxz_ );

    Ryvysy->apply_to_field( yvx_, ysyx_ );
    Ryvysy->apply_to_field( yvy_, ysyy_ );
    Ryvysy->apply_to_field( yvz_, ysyz_ );

    Ryvysz->apply_to_field( yvx_, yszx_ );
    Ryvysz->apply_to_field( yvy_, yszy_ );
    Ryvysz->apply_to_field( yvz_, yszz_ );


//     const InterpSVolYSurf* const Rsvys = SpatialOpDatabase<InterpSVolYSurf>::self().retrieve_operator();
//     Rsvys->apply_to_field( svx_, ysx_ );
//     Rsvys->apply_to_field( svy_, ysy_ );
//     Rsvys->apply_to_field( svy_, ysz_ );
  }

  // z-cell grid coordinates
  {
//     const InterpSVolZVol* const Rsvzv = SpatialOpDatabase<InterpSVolZVol>::self().retrieve_operator();
//     Rsvzv->apply_to_field( svx_, zvx_ );
//     Rsvzv->apply_to_field( svy_, zvy_ );
//     Rsvzv->apply_to_field( svz_, zvz_ );


    const InterpZVolZSurfX* const Rzvzsx = SpatialOpDatabase<InterpZVolZSurfX>::self().retrieve_operator();
    const InterpZVolZSurfY* const Rzvzsy = SpatialOpDatabase<InterpZVolZSurfY>::self().retrieve_operator();
    const InterpZVolZSurfZ* const Rzvzsz = SpatialOpDatabase<InterpZVolZSurfZ>::self().retrieve_operator();
    Rzvzsx->apply_to_field( zvx_, zsxx_ );
    Rzvzsx->apply_to_field( zvy_, zsxy_ );
    Rzvzsx->apply_to_field( zvz_, zsxz_ );

    Rzvzsy->apply_to_field( zvx_, zsyx_ );
    Rzvzsy->apply_to_field( zvy_, zsyy_ );
    Rzvzsy->apply_to_field( zvz_, zsyz_ );

    Rzvzsz->apply_to_field( zvx_, zszx_ );
    Rzvzsz->apply_to_field( zvy_, zszy_ );
    Rzvzsz->apply_to_field( zvz_, zszz_ );

//     const InterpSVolZSurf* const Rsvzs = SpatialOpDatabase<InterpSVolZSurf>::self().retrieve_operator();
//     Rsvzs->apply_to_field( svx_, zsx_ );
//     Rsvzs->apply_to_field( svy_, zsy_ );
//     Rsvzs->apply_to_field( svz_, zsz_ );
  }


}
//--------------------------------------------------------------------
Grid::~Grid()
{
}

#include <EpetraExt_VectorOut.h>
#include <EpetraExt_RowMatrixOut.h>

void
Grid::write() const
{
  /*
  EpetraExt::VectorToMatrixMarketFile( "svx.mm", svx_.get_linalg_vec(), "", "" );

  EpetraExt::VectorToMatrixMarketFile( "xvx.mm", xvx_.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "yvx.mm", yvx_.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "zvx.mm", zvx_.get_linalg_vec(), "", "" );

  EpetraExt::VectorToMatrixMarketFile( "svy.mm", svy_.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "xvy.mm", xvy_.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "yvy.mm", yvy_.get_linalg_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "zvy.mm", zvy_.get_linalg_vec(), "", "" );
  */
}

}
}

#endif
