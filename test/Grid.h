#ifndef FVTEST_GRID_H
#define FVTEST_GRID_H

#include <spatialops/structured/FVStaggered.h>

namespace SpatialOps{

  class OperatorDatabase;

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
  Grid( const IntVec& dim,
        const std::vector<double>& spacing,
        const std::vector<bool>& bcPlusFlag,
        const OperatorDatabase& opDB );
  ~Grid();

  const IntVec& extent() const{return dim_;}

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

  const IntVec dim_;

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

} // namespace structured
} // namespace SpatialOps

#endif
