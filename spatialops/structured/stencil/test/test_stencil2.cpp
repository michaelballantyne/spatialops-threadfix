#include <spatialops/SpatialOpsTools.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/OperatorDatabase.h>

#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/stencil/StencilBuilder.h>
#include <spatialops/structured/stencil/Stencil2.h>
#include <spatialops/structured/stencil/FDStencil2.h>

#include "test_stencil_helper.h"
#include <test/TestHelper.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;

#include <stdexcept>

//--------------------------------------------------------------------

template< typename Vol >
bool
run_variants( const IntVec npts,
              const bool* bcPlus )
{
  TestHelper status(true);

  typedef typename FaceTypes<Vol>::XFace  XFace;
  typedef typename FaceTypes<Vol>::YFace  YFace;
  typedef typename FaceTypes<Vol>::ZFace  ZFace;

  const double length = 10.0;

  if( npts[0]>1 ){
    status( run_convergence<Interpolant,Vol,XFace,XDIR>( npts, bcPlus, length, 2.0 ), "x-Interpolant" );
    status( run_convergence<Gradient,   Vol,XFace,XDIR>( npts, bcPlus, length, 2.0 ), "x-Gradient"    );
    status( run_convergence<Divergence, XFace,Vol,XDIR>( npts, bcPlus, length, 2.0 ), "x-Divergence"  );
  }

  if( npts[1]>1 ){
    status( run_convergence<Interpolant,Vol,YFace,YDIR>( npts, bcPlus, length, 2.0 ), "y-Interpolant" );
    status( run_convergence<Gradient,   Vol,YFace,YDIR>( npts, bcPlus, length, 2.0 ), "y-Gradient"    );
    status( run_convergence<Divergence, YFace,Vol,YDIR>( npts, bcPlus, length, 2.0 ), "y-Divergence"  );
  }

  if( npts[2]>1 ){
    status( run_convergence<Interpolant,Vol,ZFace,ZDIR>( npts, bcPlus, length, 2.0 ), "z-Interpolant" );
    status( run_convergence<Gradient,   Vol,ZFace,ZDIR>( npts, bcPlus, length, 2.0 ), "z-Gradient"    );
    status( run_convergence<Divergence, ZFace,Vol,ZDIR>( npts, bcPlus, length, 2.0 ), "z-Divergence"  );
  }
  return status.ok();
}

//--------------------------------------------------------------------

#define TEST_EXTENTS( SRC, DEST,                                                                                                 \
                      DirT,                                                                                                      \
                      S2Ox,  S2Oy,  S2Oz,                                                                                        \
                      DOx,   DOy,   DOz,                                                                                         \
                      ULx,   ULy,   ULz,                                                                                         \
                      name )                                                                                                     \
    {                                                                                                                            \
      using std::string;                                                                                                         \
      typedef ExtentsAndOffsets<SRC,DEST> Extents;                                                                               \
      status( IsSameType< Extents::Dir,        DirT                          >::result, string(name) + string(" dir"       ) );  \
      status( IsSameType< Extents::Src2Offset, IndexTriplet<S2Ox,S2Oy,S2Oz > >::result, string(name) + string(" s2 offset" ) );  \
      status( IsSameType< Extents::DestOffset, IndexTriplet<DOx, DOy, DOz  > >::result, string(name) + string(" d  offset" ) );  \
      status( IsSameType< Extents::Src1Extent, IndexTriplet<ULx, ULy, ULz  > >::result, string(name) + string(" UB"        ) );  \
    }
//status( IsSameType< Extents::UpperLoopBCAug, IndexTriplet<ULBCx, ULBCy, ULBCz> >::result, string(name) + string(" UB Aug."   ) );

//-------------------------------------------------------------------

bool test_compile_time()
{
  using namespace SpatialOps;
  using namespace structured;
  using namespace s2detail;

  TestHelper status(false);

  status( IsSameType< ActiveDir< YVolField, SVolField   >::type, YDIR >::result, "YVol->SVol (y)" );
  status( IsSameType< ActiveDir< YVolField, XSurfYField >::type, XDIR >::result, "YVol->ZSY  (x)" );
  status( IsSameType< ActiveDir< YVolField, ZSurfYField >::type, ZDIR >::result, "YVol->ZSY  (z)" );

  status( IsSameType< ActiveDir< ZVolField, SVolField   >::type, ZDIR >::result, "ZVol->SVol (z)" );
  status( IsSameType< ActiveDir< ZVolField, XSurfZField >::type, XDIR >::result, "ZVol->XSZ  (x)" );
  status( IsSameType< ActiveDir< ZVolField, YSurfZField >::type, YDIR >::result, "ZVol->YSZ  (y)" );

  status( IsSameType< ActiveDir< SVolField, XVolField >::type, XDIR >::result, "SVol->XVol (x)" );
  status( IsSameType< ActiveDir< SVolField, YVolField >::type, YDIR >::result, "SVol->YVol (y)" );
  status( IsSameType< ActiveDir< SVolField, ZVolField >::type, ZDIR >::result, "SVol->ZVol (z)" );

  TEST_EXTENTS( SVolField, SSurfXField, XDIR,  1,0,0,  1,0,0,  -1, 0, 0,  "SVol->SSX" )
  TEST_EXTENTS( SVolField, SSurfYField, YDIR,  0,1,0,  0,1,0,   0,-1, 0,  "SVol->SSY" )
  TEST_EXTENTS( SVolField, SSurfZField, ZDIR,  0,0,1,  0,0,1,   0, 0,-1,  "SVol->SSZ" )
  TEST_EXTENTS( SSurfXField, SVolField, XDIR,  1,0,0,  0,0,0,  -1, 0, 0,  "SSX->SVol" )
  TEST_EXTENTS( SSurfYField, SVolField, YDIR,  0,1,0,  0,0,0,   0,-1, 0,  "SSY->SVol" )
  TEST_EXTENTS( SSurfZField, SVolField, ZDIR,  0,0,1,  0,0,0,   0, 0,-1,  "SSZ->SVol" )

  TEST_EXTENTS( XVolField, XSurfXField, XDIR,  1,0,0,  0,0,0,  -1, 0, 0,  "XVol->XSX" )
  TEST_EXTENTS( XVolField, XSurfYField, YDIR,  0,1,0,  0,1,0,   0,-1, 0,  "XVol->XSY" )
  TEST_EXTENTS( XVolField, XSurfZField, ZDIR,  0,0,1,  0,0,1,   0, 0,-1,  "XVol->XSZ" )

  TEST_EXTENTS( XSurfXField, XVolField, XDIR,  1,0,0,  1,0,0,  -1, 0, 0,  "XSX->XVol" )
  TEST_EXTENTS( XSurfYField, XVolField, YDIR,  0,1,0,  0,0,0,   0,-1, 0,  "XSY->XVol" )
  TEST_EXTENTS( XSurfZField, XVolField, ZDIR,  0,0,1,  0,0,0,   0, 0,-1,  "XSZ->XVol" )

  TEST_EXTENTS( YVolField, YSurfXField, XDIR,  1,0,0,  1,0,0,  -1, 0, 0,  "YVol->YSX" )
  TEST_EXTENTS( YVolField, YSurfYField, YDIR,  0,1,0,  0,0,0,   0,-1, 0,  "YVol->YSY" )
  TEST_EXTENTS( YVolField, YSurfZField, ZDIR,  0,0,1,  0,0,1,   0, 0,-1,  "YVol->YSZ" )
  TEST_EXTENTS( YSurfXField, YVolField, XDIR,  1,0,0,  0,0,0,  -1, 0, 0,  "YSX->YVol" )
  TEST_EXTENTS( YSurfYField, YVolField, YDIR,  0,1,0,  0,1,0,   0,-1, 0,  "YSY->YVol" )
  TEST_EXTENTS( YSurfZField, YVolField, ZDIR,  0,0,1,  0,0,0,   0, 0,-1,  "YSZ->YVol" )

  TEST_EXTENTS( ZVolField, ZSurfXField, XDIR,  1,0,0,  1,0,0,  -1, 0, 0,  "ZVol->ZSX" )
  TEST_EXTENTS( ZVolField, ZSurfYField, YDIR,  0,1,0,  0,1,0,   0,-1, 0,  "ZVol->ZSY" )
  TEST_EXTENTS( ZVolField, ZSurfZField, ZDIR,  0,0,1,  0,0,0,   0, 0,-1,  "ZVol->ZSZ" )
  TEST_EXTENTS( ZSurfXField, ZVolField, XDIR,  1,0,0,  0,0,0,  -1, 0, 0,  "ZSX->ZVol" )
  TEST_EXTENTS( ZSurfYField, ZVolField, YDIR,  0,1,0,  0,0,0,   0,-1, 0,  "ZSY->ZVol" )
  TEST_EXTENTS( ZSurfZField, ZVolField, ZDIR,  0,0,1,  0,0,1,   0, 0,-1,  "ZSZ->ZVol" )

  TEST_EXTENTS( XVolField, SVolField,   XDIR,  1,0,0,  0,0,0,  -1, 0, 0,  "XVol->SVol" )
  TEST_EXTENTS( XVolField, YSurfXField, YDIR,  0,1,0,  0,1,0,   0,-1, 0,  "XVol->YSX"  )
  TEST_EXTENTS( XVolField, ZSurfXField, ZDIR,  0,0,1,  0,0,1,   0, 0,-1,  "XVol->ZSX"  )

  return status.ok();
}

//--------------------------------------------------------------------

int main( int iarg, char* carg[] )
{ int nx, ny, nz;
  bool bcplus[] = { false, false, false };
  {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message\n" )
      ( "nx",   po::value<int>(&nx)->default_value(11), "number of points in x-dir for base mesh" )
      ( "ny",   po::value<int>(&ny)->default_value(11), "number of points in y-dir for base mesh" )
      ( "nz",   po::value<int>(&nz)->default_value(11), "number of points in z-dir for base mesh" )
      ( "bcx",  "physical boundary on +x side?" )
      ( "bcy",  "physical boundary on +y side?" )
      ( "bcz",  "physical boundary on +z side?" );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    if( args.count("bcx") ) bcplus[0] = true;
    if( args.count("bcy") ) bcplus[1] = true;
    if( args.count("bcz") ) bcplus[2] = true;

    if( args.count("help") ){
      cout << desc << endl
           << "Examples:" << endl
           << " test_stencil2 --nx 5 --ny 10 --nz 3 --bcx" << endl
           << " test_stencil2 --bcx --bcy --bcz" << endl
           << " test_stencil2 --nx 50 --bcz" << endl
           << endl;
      return -1;
    }
  }

  TestHelper status( true );
  const IntVec npts(nx,ny,nz);
  {
    const std::string bcx = bcplus[0] ? "ON" : "OFF";
    const std::string bcy = bcplus[1] ? "ON" : "OFF";
    const std::string bcz = bcplus[2] ? "ON" : "OFF";
    cout << "Run information: " << endl
         << "  bcx    : " << bcx << endl
         << "  bcy    : " << bcy << endl
         << "  bcz    : " << bcz << endl
         << "  domain : " << npts << endl
         << endl;
  }

  status( test_compile_time(), "Compile time type introspection tests" );
  cout << endl;
  const double length = 10.0;

//   status( run_convergence<Interpolant,SVolField,SSurfXField,XDIR>(npts,bcplus,length,2.0) );
//   status( run_convergence<Interpolant,SVolField,SSurfYField,YDIR>(npts,bcplus,length,2.0) );
//   status( run_convergence<Divergence,SSurfXField,SVolField,XDIR>(npts,bcplus,length,2.0) );
//   status( run_convergence<Divergence,SSurfYField,SVolField,YDIR>(npts,bcplus,length,2.0) );
//   status( run_convergence<Divergence,SSurfYField,SVolField,YDIR>(npts,bcplus,length,2.0) );
//   status( run_convergence<Interpolant,XVolField,XSurfXField,XDIR>(npts,bcplus,length,2.0) );
//   status( run_convergence<Divergence,XSurfYField,XVolField,YDIR>(npts,bcplus,length,2.0) );
//return 0;

  try{
#ifdef ENABLE_CUDA
    status( run_variants< SVolField >( npts, bcplus ), "SVol operators" );
    if( npts[0] > 1 ) status( run_variants< XVolField >( npts, bcplus ), "XVol operators" );
    if( npts[1] > 1 ) status( run_variants< YVolField >( npts, bcplus ), "YVol operators" );
    if( npts[2] > 1 ) status( run_variants< ZVolField >( npts, bcplus ), "ZVol operators" );

    if( npts[0]>1 & npts[1]>1 ) status( run_convergence< Interpolant, XVolField,   YSurfXField, YDIR >( npts, bcplus, length, 2.0 ), "InterpXVolYSurfX" );
    if( npts[0]>1 & npts[2]>1 ) status( run_convergence< Interpolant, XVolField,   ZSurfXField, ZDIR >( npts, bcplus, length, 2.0 ), "InterpXVolZSurfX" );

    if( npts[0]>1 & npts[1]>1 ) status( run_convergence< Gradient,    XVolField,   YSurfXField, YDIR >( npts, bcplus, length, 2.0 ), "GradXVolYSurfX" );
    if( npts[0]>1 & npts[2]>1 ) status( run_convergence< Gradient,    XVolField,   ZSurfXField, ZDIR >( npts, bcplus, length, 2.0 ), "GradXVolZSurfX" );

    if( npts[1]>1 & npts[0]>1 ) status( run_convergence< Interpolant, YVolField,   XSurfYField, XDIR >( npts, bcplus, length, 2.0 ), "InterpYVolXSurfY" );
    if( npts[1]>1 & npts[2]>1 ) status( run_convergence< Interpolant, YVolField,   ZSurfYField, ZDIR >( npts, bcplus, length, 2.0 ), "InterpYVolZSurfY" );

    if( npts[1]>1 & npts[0]>1 ) status( run_convergence< Gradient,    YVolField,   XSurfYField, XDIR >( npts, bcplus, length, 2.0 ), "GradYVolXSurfY" );
    if( npts[1]>1 & npts[2]>1 ) status( run_convergence< Gradient,    YVolField,   ZSurfYField, ZDIR >( npts, bcplus, length, 2.0 ), "GradYVolZSurfY" );

    if( npts[2]>1 & npts[0]>1 ) status( run_convergence< Interpolant, ZVolField,   XSurfZField, XDIR >( npts, bcplus, length, 2.0 ), "InterpZVolXSurfZ" );
    if( npts[2]>1 & npts[1]>1 ) status( run_convergence< Interpolant, ZVolField,   YSurfZField, YDIR >( npts, bcplus, length, 2.0 ), "InterpZVolYSurfZ" );

    if( npts[2]>1 & npts[0]>1 ) status( run_convergence< Gradient,    ZVolField,   XSurfZField, XDIR >( npts, bcplus, length, 2.0 ), "GradZVolXSurfZ" );
    if( npts[2]>1 & npts[1]>1 ) status( run_convergence< Gradient,    ZVolField,   YSurfZField, YDIR >( npts, bcplus, length, 2.0 ), "GradZVolYSurfZ" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, SVolField,   XVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpSVolXVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, SVolField,   YVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpSVolYVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, SVolField,   ZVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpSVolZVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, XVolField,   SVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpXVolSVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, YVolField,   SVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpYVolSVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, ZVolField,   SVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpZVolSVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, XSurfXField, XVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpXSXXVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, XSurfYField, XVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpXSYXVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, XSurfZField, XVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpXSZXVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, YSurfXField, YVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpYSXYVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, YSurfYField, YVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpYSYYVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, YSurfZField, YVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpYSZYVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, ZSurfXField, ZVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpZSXZVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, ZSurfYField, ZVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpZSYZVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, ZSurfZField, ZVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpZSZZVol" );
#else
    status( run_variants< SVolField >( npts, bcplus ), "SVol operators" );
    if( npts[0] > 1 ) status( run_variants< XVolField >( npts, bcplus ), "XVol operators" );
    if( npts[1] > 1 ) status( run_variants< YVolField >( npts, bcplus ), "YVol operators" );
    if( npts[2] > 1 ) status( run_variants< ZVolField >( npts, bcplus ), "ZVol operators" );

    if( npts[0]>1 & npts[1]>1 ) status( run_convergence< Interpolant, XVolField,   YSurfXField, YDIR >( npts, bcplus, length, 2.0 ), "InterpXVolYSurfX" );
    if( npts[0]>1 & npts[2]>1 ) status( run_convergence< Interpolant, XVolField,   ZSurfXField, ZDIR >( npts, bcplus, length, 2.0 ), "InterpXVolZSurfX" );

    if( npts[0]>1 & npts[1]>1 ) status( run_convergence< Gradient,    XVolField,   YSurfXField, YDIR >( npts, bcplus, length, 2.0 ), "GradXVolYSurfX" );
    if( npts[0]>1 & npts[2]>1 ) status( run_convergence< Gradient,    XVolField,   ZSurfXField, ZDIR >( npts, bcplus, length, 2.0 ), "GradXVolZSurfX" );

    if( npts[1]>1 & npts[0]>1 ) status( run_convergence< Interpolant, YVolField,   XSurfYField, XDIR >( npts, bcplus, length, 2.0 ), "InterpYVolXSurfY" );
    if( npts[1]>1 & npts[2]>1 ) status( run_convergence< Interpolant, YVolField,   ZSurfYField, ZDIR >( npts, bcplus, length, 2.0 ), "InterpYVolZSurfY" );

    if( npts[1]>1 & npts[0]>1 ) status( run_convergence< Gradient,    YVolField,   XSurfYField, XDIR >( npts, bcplus, length, 2.0 ), "GradYVolXSurfY" );
    if( npts[1]>1 & npts[2]>1 ) status( run_convergence< Gradient,    YVolField,   ZSurfYField, ZDIR >( npts, bcplus, length, 2.0 ), "GradYVolZSurfY" );

    if( npts[2]>1 & npts[0]>1 ) status( run_convergence< Interpolant, ZVolField,   XSurfZField, XDIR >( npts, bcplus, length, 2.0 ), "InterpZVolXSurfZ" );
    if( npts[2]>1 & npts[1]>1 ) status( run_convergence< Interpolant, ZVolField,   YSurfZField, YDIR >( npts, bcplus, length, 2.0 ), "InterpZVolYSurfZ" );

    if( npts[2]>1 & npts[0]>1 ) status( run_convergence< Gradient,    ZVolField,   XSurfZField, XDIR >( npts, bcplus, length, 2.0 ), "GradZVolXSurfZ" );
    if( npts[2]>1 & npts[1]>1 ) status( run_convergence< Gradient,    ZVolField,   YSurfZField, YDIR >( npts, bcplus, length, 2.0 ), "GradZVolYSurfZ" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, SVolField,   XVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpSVolXVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, SVolField,   YVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpSVolYVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, SVolField,   ZVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpSVolZVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, XVolField,   SVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpXVolSVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, YVolField,   SVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpYVolSVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, ZVolField,   SVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpZVolSVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, XSurfXField, XVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpXSXXVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, XSurfYField, XVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpXSYXVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, XSurfZField, XVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpXSZXVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, YSurfXField, YVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpYSXYVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, YSurfYField, YVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpYSYYVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, YSurfZField, YVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpYSZYVol" );

    if( npts[0]>1 )             status( run_convergence< Interpolant, ZSurfXField, ZVolField,   XDIR >( npts, bcplus, length, 2.0 ), "InterpZSXZVol" );
    if( npts[1]>1 )             status( run_convergence< Interpolant, ZSurfYField, ZVolField,   YDIR >( npts, bcplus, length, 2.0 ), "InterpZSYZVol" );
    if( npts[2]>1 )             status( run_convergence< Interpolant, ZSurfZField, ZVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "InterpZSZZVol" );

    if( npts[0]>1 )             status( run_convergence< InterpolantX, SVolField,  SVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Interp SVol->SVol (X)" );
    if( npts[1]>1 )             status( run_convergence< InterpolantY, SVolField,  SVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Interp SVol->SVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< InterpolantZ, SVolField,  SVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Interp SVol->SVol (Z)" );
    if( npts[0]>1 )             status( run_convergence< InterpolantX, XVolField,  XVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Interp XVol->XVol (X)" );
    if( npts[1]>1 )             status( run_convergence< InterpolantY, XVolField,  XVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Interp XVol->XVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< InterpolantZ, XVolField,  XVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Interp XVol->XVol (Z)" );
    if( npts[0]>1 )             status( run_convergence< InterpolantX, YVolField,  YVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Interp YVol->YVol (X)" );
    if( npts[1]>1 )             status( run_convergence< InterpolantY, YVolField,  YVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Interp YVol->YVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< InterpolantZ, YVolField,  YVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Interp YVol->YVol (Z)" );
    if( npts[0]>1 )             status( run_convergence< InterpolantX, ZVolField,  ZVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Interp ZVol->ZVol (X)" );
    if( npts[1]>1 )             status( run_convergence< InterpolantY, ZVolField,  ZVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Interp ZVol->ZVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< InterpolantZ, ZVolField,  ZVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Interp ZVol->ZVol (Z)" );

    if( npts[0]>1 )             status( run_convergence< GradientX,    SVolField,  SVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Grad SVol->SVol (X)" );
    if( npts[1]>1 )             status( run_convergence< GradientY,    SVolField,  SVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Grad SVol->SVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< GradientZ,    SVolField,  SVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Grad SVol->SVol (Z)" );
    if( npts[0]>1 )             status( run_convergence< GradientX,    XVolField,  XVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Grad XVol->XVol (X)" );
    if( npts[1]>1 )             status( run_convergence< GradientY,    XVolField,  XVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Grad XVol->XVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< GradientZ,    XVolField,  XVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Grad XVol->XVol (Z)" );
    if( npts[0]>1 )             status( run_convergence< GradientX,    YVolField,  YVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Grad YVol->YVol (X)" );
    if( npts[1]>1 )             status( run_convergence< GradientY,    YVolField,  YVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Grad YVol->YVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< GradientZ,    YVolField,  YVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Grad YVol->YVol (Z)" );
    if( npts[0]>1 )             status( run_convergence< GradientX,    ZVolField,  ZVolField,   XDIR >( npts, bcplus, length, 2.0 ), "Grad ZVol->ZVol (X)" );
    if( npts[1]>1 )             status( run_convergence< GradientY,    ZVolField,  ZVolField,   YDIR >( npts, bcplus, length, 2.0 ), "Grad ZVol->ZVol (Y)" );
    if( npts[2]>1 )             status( run_convergence< GradientZ,    ZVolField,  ZVolField,   ZDIR >( npts, bcplus, length, 2.0 ), "Grad ZVol->ZVol (Z)" );
#endif
    if( status.ok() ){
      cout << "ALL TESTS PASSED :)" << endl;
      return 0;
    }
  }
  catch( std::runtime_error& e ){
    cout << e.what() << endl;
  }

  cout << "******************************" << endl
       << " At least one test FAILED! :(" << endl
       << "******************************" << endl;
  return -1;
}

//--------------------------------------------------------------------
