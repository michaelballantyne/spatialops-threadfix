#include "Grid.h"

#include <spatialops/WriteMatlab.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FVTools.h>

#include <spatialops/structured/matrix/FVStaggeredInterpolant.h>
#include <spatialops/structured/matrix/FVStaggeredGradient.h>
#include <spatialops/structured/matrix/FVStaggeredDivergence.h>
#include <spatialops/structured/matrix/FVStaggeredScratch.h>
#include <spatialops/structured/matrix/FVTopHatFilter.h>
#include <spatialops/structured/matrix/FVRestrictOp.h>
#include <spatialops/structured/matrix/FVStaggeredOperatorTypes.h>

namespace SpatialOps{
namespace structured{

  Grid::Grid( const IntVec& dim,
              const std::vector<double>& spacing,
              const std::vector<bool>& bcPlusFlag,
              const OperatorDatabase& opDB )
    : dim_( dim ),

      svx_( get_window_with_ghost<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      svy_( get_window_with_ghost<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      svz_( get_window_with_ghost<SVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

      ssxx_( get_window_with_ghost<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ssxy_( get_window_with_ghost<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ssxz_( get_window_with_ghost<SSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ssyx_( get_window_with_ghost<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ssyy_( get_window_with_ghost<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ssyz_( get_window_with_ghost<SSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      sszx_( get_window_with_ghost<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      sszy_( get_window_with_ghost<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      sszz_( get_window_with_ghost<SSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),


      xvx_( get_window_with_ghost<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xvy_( get_window_with_ghost<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xvz_( get_window_with_ghost<XVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),


      xsxx_( get_window_with_ghost<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xsxy_( get_window_with_ghost<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xsxz_( get_window_with_ghost<XSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xsyx_( get_window_with_ghost<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xsyy_( get_window_with_ghost<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xsyz_( get_window_with_ghost<XSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xszx_( get_window_with_ghost<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xszy_( get_window_with_ghost<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      xszz_( get_window_with_ghost<XSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

      yvx_( get_window_with_ghost<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      yvy_( get_window_with_ghost<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      yvz_( get_window_with_ghost<YVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

      ysxx_( get_window_with_ghost<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ysxy_( get_window_with_ghost<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ysxz_( get_window_with_ghost<YSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ysyx_( get_window_with_ghost<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ysyy_( get_window_with_ghost<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      ysyz_( get_window_with_ghost<YSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      yszx_( get_window_with_ghost<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      yszy_( get_window_with_ghost<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      yszz_( get_window_with_ghost<YSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),


      zvx_( get_window_with_ghost<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zvy_( get_window_with_ghost<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zvz_( get_window_with_ghost<ZVolField >(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),

      zsxx_( get_window_with_ghost<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zsxy_( get_window_with_ghost<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zsxz_( get_window_with_ghost<ZSurfXField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zsyx_( get_window_with_ghost<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zsyy_( get_window_with_ghost<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zsyz_( get_window_with_ghost<ZSurfYField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zszx_( get_window_with_ghost<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zszy_( get_window_with_ghost<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL ),
      zszz_( get_window_with_ghost<ZSurfZField>(dim,bcPlusFlag[0],bcPlusFlag[1],bcPlusFlag[2]), NULL )
  {
    using namespace SpatialOps;
    using namespace structured;

    {
      SVolField::iterator isvx=svx_.begin();
      SVolField::iterator isvy=svy_.begin();
      SVolField::iterator isvz=svz_.begin();

      const int ihi = get_nx_with_ghost<SVolField>(dim[0],bcPlusFlag[0]);
      const int jhi = get_ny_with_ghost<SVolField>(dim[1],bcPlusFlag[1]);
      const int khi = get_nz_with_ghost<SVolField>(dim[2],bcPlusFlag[2]);

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
      const int ihi = get_nx_with_ghost<XVolField>(dim[0],bcPlusFlag[0]);
      const int jhi = get_ny_with_ghost<XVolField>(dim[1],bcPlusFlag[1]);
      const int khi = get_nz_with_ghost<XVolField>(dim[2],bcPlusFlag[2]);

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
      const int ihi = get_nx_with_ghost<YVolField>(dim[0],bcPlusFlag[0]);
      const int jhi = get_ny_with_ghost<YVolField>(dim[1],bcPlusFlag[1]);
      const int khi = get_nz_with_ghost<YVolField>(dim[2],bcPlusFlag[2]);

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
      const int ihi = get_nx_with_ghost<ZVolField>(dim[0],bcPlusFlag[0]);
      const int jhi = get_ny_with_ghost<ZVolField>(dim[1],bcPlusFlag[1]);
      const int khi = get_nz_with_ghost<ZVolField>(dim[2],bcPlusFlag[2]);

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
      write_matlab(svx_,"svx");
      write_matlab(ssxx_,"ssxx");

      write_matlab(xvx_,"xvx");
      write_matlab(xsxx_,"xsxx");
      write_matlab(xsyx_,"xsyx");
      write_matlab(xszx_,"xszx");

      if( dim_[1]>1 ){
        write_matlab(ssxy_,"ssxy");
        write_matlab(xvy_ ,"xvy");
        write_matlab(xsxy_,"xsxy");
        write_matlab(xsyy_,"xsyy");
        write_matlab(xszy_,"xszy");
      }
      if( dim_[2]>1 ){
        write_matlab(ssxz_,"ssxz");
        write_matlab(xvz_ ,"xvz");
        write_matlab(xsxz_,"xsxz");
        write_matlab(xsyz_,"xsyz");
        write_matlab(xszz_,"xszz");
      }
    }

    if( dim_[1]>1 ){
      if( dim_[0]>1 ){
        write_matlab(ssyx_,"ssyx");

        write_matlab(yvx_ ,"yvx");
        write_matlab(ysxx_,"ysxx");
        write_matlab(ysyx_,"ysyx");
        write_matlab(yszx_,"yszx");
      }

      write_matlab(svy_ ,"svy");
      write_matlab(ssyy_,"ssyy");

      write_matlab(yvy_ ,"yvy");
      write_matlab(ysxy_,"ysxy");
      write_matlab(ysyy_,"ysyy");
      write_matlab(yszy_,"yszy");

      if( dim_[2]>1 ){
        write_matlab(ssyz_,"ssyz");
        write_matlab(yvz_ ,"yvz");
        write_matlab(ysxz_,"ysxz");
        write_matlab(ysyz_,"ysyz");
        write_matlab(yszz_,"yszz");
      }
    }

    if( dim_[2]>1 ){
      if( dim_[0]>1 ){
        write_matlab(sszx_,"sszx");
        write_matlab(zvx_ ,"zvx");
        write_matlab(zsxx_,"zsxx");
        write_matlab(zsyx_,"zsyx");
        write_matlab(zszx_,"zszx");
      }

      if( dim_[1]>1 ){
        write_matlab(sszy_,"sszy");
        write_matlab(zvy_ ,"zvy");
        write_matlab(zsxy_,"zsxy");
        write_matlab(zsyy_,"zsyy");
        write_matlab(zszy_,"zszy");
      }

      write_matlab(svz_ ,"svz");
      write_matlab(sszz_,"sszz");

      write_matlab(yvz_ ,"yvz");
      write_matlab(ysxz_,"ysxz");
      write_matlab(ysyz_,"ysyz");
      write_matlab(yszz_,"yszz");
    }

  }


} // namespace structured
} // namespace SpatialOps
