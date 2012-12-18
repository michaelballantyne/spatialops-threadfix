#include <spatialops/SpatialOpsTools.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/OperatorDatabase.h>

#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/stencil/StencilBuilder.h>

#include "ReferenceStencil.h"
#include <test/TestHelper.h>

#include <spatialops/FieldExpressions.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;

#include <stdexcept>
using std::cout;
using std::endl;

//--------------------------------------------------------------------

template<typename Field>
void initialize_field(Field & f) {
    typename Field::iterator i = f.begin();
    MemoryWindow mw = f.window_with_ghost();
    int xLength = mw.extent(0);
    int yLength = mw.extent(1);
    int zLength = mw.extent(2);

    for(int x = 1; x <= xLength; x++) {
        for(int y = 1; y <= yLength; y++) {
            for(int z = 1; z <= zLength; z++) {
                *i = sin(x +
                         y * xLength +
                         z * xLength * yLength);
                i++;
            };
        };
    };
};

//--------------------------------------------------------------------

template<typename OpType, typename SrcType, typename DestType>
bool test_stencil2(OperatorDatabase & opdb,
                   SrcType  const & src,
                   DestType       & ref,
                   DestType       & test) {
    //zero out result fields:
    ref <<= 0.0;
    test <<= 0.0;

    //get operator:
    typedef typename SpatialOps::structured::OperatorTypeBuilder<OpType,SrcType,DestType>::type Op;
    const Op* const op = opdb.retrieve_operator<Op>();

    //run reference:
    ref_stencil2_apply_to_field(op->get_minus_coef(),
                                op->get_plus_coef(),
                                src,
                                ref);

    //run operator:
    op->apply_to_field(src, test);

    return (test == ref);
};

//--------------------------------------------------------------------

template<typename OpType, typename SrcType, typename DestType>
bool test_stencil4(OperatorDatabase & opdb,
                   SrcType  const & src,
                   DestType       & ref,
                   DestType       & test) {
    //zero out result fields:
    ref <<= 0.0;
    test <<= 0.0;

    //get operator:
    typedef typename SpatialOps::structured::OperatorTypeBuilder<OpType,SrcType,DestType>::type Op;
    const Op* const op = opdb.retrieve_operator<Op>();

    //run reference:
    ref_stencil4_apply_to_field(op->get_coef1(),
                                op->get_coef2(),
                                op->get_coef3(),
                                op->get_coef4(),
                                src,
                                ref);

    //run operator:
    op->apply_to_field(src, test);

    return (test == ref);
};

//--------------------------------------------------------------------

template<typename OpType, typename FieldType>
bool test_fd_stencil2(OperatorDatabase & opdb,
                      FieldType const & src,
                      FieldType       & ref,
                      FieldType       & test) {
    //zero out result fields:
    ref <<= 0.0;
    test <<= 0.0;

    //get operator:
    typedef typename SpatialOps::structured::OperatorTypeBuilder<OpType,FieldType,FieldType>::type Op;
    const Op* const op = opdb.retrieve_operator<Op>();

    //run reference:
    ref_fd_stencil2_apply_to_field<OpType,FieldType>(op->get_minus_coef(),
                                                     op->get_plus_coef(),
                                                     src,
                                                     ref);

    //run operator:
    op->apply_to_field(src, test);

    return (test == ref);
};

//--------------------------------------------------------------------

template<typename VolField, typename SurfXField, typename SurfYField, typename SurfZField>
inline bool test_basic_stencils(const IntVec npts, OperatorDatabase & opdb,
                                const VolField & srcVol,  const SurfXField & srcSurfX,  const SurfYField & srcSurfY,  const SurfZField & srcSurfZ,
                                      VolField & refVol,        SurfXField & refSurfX,        SurfYField & refSurfY,        SurfZField & refSurfZ,
                                      VolField & testVol,       SurfXField & testSurfX,       SurfYField & testSurfY,       SurfZField & testSurfZ) {
    TestHelper status(true);

    if( npts[0] > 1) {
        status( test_stencil2<Interpolant, VolField,   SurfXField>(opdb, srcVol,   refSurfX, testSurfX), "Interpolant VolField -> SurfXField (2)" );
        status( test_stencil2<Gradient,    VolField,   SurfXField>(opdb, srcVol,   refSurfX, testSurfX), "Gradient    VolField -> SurfXField (2)" );
        status( test_stencil2<Divergence,  SurfXField, VolField>  (opdb, srcSurfX, refVol,   testVol),   "Divergence  SurfXField -> VolField (2)" );
    };

    if( npts[1] > 1) {
        status( test_stencil2<Interpolant, VolField,   SurfYField>(opdb, srcVol,   refSurfY, testSurfY), "Interpolant VolField -> SurfYField (2)" );
        status( test_stencil2<Gradient,    VolField,   SurfYField>(opdb, srcVol,   refSurfY, testSurfY), "Gradient    VolField -> SurfYField (2)" );
        status( test_stencil2<Divergence,  SurfYField, VolField>  (opdb, srcSurfY, refVol,   testVol),   "Divergence  SurfYField -> VolField (2)" );
    };

    if( npts[2] > 1) {
        status( test_stencil2<Interpolant, VolField,   SurfZField>(opdb, srcVol,   refSurfZ, testSurfZ), "Interpolant VolField -> SurfZField (2)" );
        status( test_stencil2<Gradient,    VolField,   SurfZField>(opdb, srcVol,   refSurfZ, testSurfZ), "Gradient    VolField -> SurfZField (2)" );
        status( test_stencil2<Divergence,  SurfZField, VolField>  (opdb, srcSurfZ, refVol,   testVol),   "Divergence  SurfZField -> VolField (2)" );
    };

    return status.ok();
};

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
//-------------------------------------------------------------------

bool test_compile_time()
{
  using namespace SpatialOps;
  using namespace structured;
  using namespace RefStencil2Detail;

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
{
  int nx, ny, nz;
  bool bc[] = { false, false, false };

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

    if( args.count("bcx") ) bc[0] = true;
    if( args.count("bcy") ) bc[1] = true;
    if( args.count("bcz") ) bc[2] = true;

    if( args.count("help") ){
      cout << desc << endl
           << "Examples:" << endl
           << " test_stencil --nx 5 --ny 10 --nz 3 --bcx" << endl
           << " test_stencil --bcx --bcy --bcz" << endl
           << " test_stencil --nx 50 --bcz" << endl
           << endl;
      return -1;
    }
  }

  TestHelper status( true );
  const IntVec npts(nx,ny,nz);

  {
    const std::string bcx = bc[0] ? "ON" : "OFF";
    const std::string bcy = bc[1] ? "ON" : "OFF";
    const std::string bcz = bc[2] ? "ON" : "OFF";
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

  const MemoryWindow mwSVol    = get_window_with_ghost<SVolField>  (npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwSSurfX  = get_window_with_ghost<SSurfXField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwSSurfY  = get_window_with_ghost<SSurfYField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwSSurfZ  = get_window_with_ghost<SSurfZField>(npts, bc[0], bc[1], bc[2]);

  const MemoryWindow mwXVol    = get_window_with_ghost<XVolField>  (npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwXSurfX  = get_window_with_ghost<XSurfXField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwXSurfY  = get_window_with_ghost<XSurfYField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwXSurfZ  = get_window_with_ghost<XSurfZField>(npts, bc[0], bc[1], bc[2]);

  const MemoryWindow mwYVol    = get_window_with_ghost<YVolField>  (npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwYSurfX  = get_window_with_ghost<YSurfXField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwYSurfY  = get_window_with_ghost<YSurfYField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwYSurfZ  = get_window_with_ghost<YSurfZField>(npts, bc[0], bc[1], bc[2]);

  const MemoryWindow mwZVol    = get_window_with_ghost<ZVolField>  (npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwZSurfX  = get_window_with_ghost<ZSurfXField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwZSurfY  = get_window_with_ghost<ZSurfYField>(npts, bc[0], bc[1], bc[2]);
  const MemoryWindow mwZSurfZ  = get_window_with_ghost<ZSurfZField>(npts, bc[0], bc[1], bc[2]);

  SVolField   srcSVol  (mwSVol,   NULL), refSVol  (mwSVol,   NULL), testSVol  (mwSVol,   NULL);
  SSurfXField srcSSurfX(mwSSurfX, NULL), refSSurfX(mwSSurfX, NULL), testSSurfX(mwSSurfX, NULL);
  SSurfYField srcSSurfY(mwSSurfY, NULL), refSSurfY(mwSSurfY, NULL), testSSurfY(mwSSurfY, NULL);
  SSurfZField srcSSurfZ(mwSSurfZ, NULL), refSSurfZ(mwSSurfZ, NULL), testSSurfZ(mwSSurfZ, NULL);

  XVolField   srcXVol  (mwXVol,   NULL), refXVol  (mwXVol,   NULL), testXVol  (mwXVol,   NULL);
  XSurfXField srcXSurfX(mwXSurfX, NULL), refXSurfX(mwXSurfX, NULL), testXSurfX(mwXSurfX, NULL);
  XSurfYField srcXSurfY(mwXSurfY, NULL), refXSurfY(mwXSurfY, NULL), testXSurfY(mwXSurfY, NULL);
  XSurfZField srcXSurfZ(mwXSurfZ, NULL), refXSurfZ(mwXSurfZ, NULL), testXSurfZ(mwXSurfZ, NULL);

  YVolField   srcYVol  (mwYVol,   NULL), refYVol  (mwYVol,   NULL), testYVol  (mwYVol,   NULL);
  YSurfXField srcYSurfX(mwYSurfX, NULL), refYSurfX(mwYSurfX, NULL), testYSurfX(mwYSurfX, NULL);
  YSurfYField srcYSurfY(mwYSurfY, NULL), refYSurfY(mwYSurfY, NULL), testYSurfY(mwYSurfY, NULL);
  YSurfZField srcYSurfZ(mwYSurfZ, NULL), refYSurfZ(mwYSurfZ, NULL), testYSurfZ(mwYSurfZ, NULL);

  ZVolField   srcZVol  (mwZVol,   NULL), refZVol  (mwZVol,   NULL), testZVol  (mwZVol,   NULL);
  ZSurfXField srcZSurfX(mwZSurfX, NULL), refZSurfX(mwZSurfX, NULL), testZSurfX(mwZSurfX, NULL);
  ZSurfYField srcZSurfY(mwZSurfY, NULL), refZSurfY(mwZSurfY, NULL), testZSurfY(mwZSurfY, NULL);
  ZSurfZField srcZSurfZ(mwZSurfZ, NULL), refZSurfZ(mwZSurfZ, NULL), testZSurfZ(mwZSurfZ, NULL);

  initialize_field(srcSVol);
  initialize_field(srcSSurfX);
  initialize_field(srcSSurfY);
  initialize_field(srcSSurfZ);

  initialize_field(srcXVol);
  initialize_field(srcXSurfX);
  initialize_field(srcXSurfY);
  initialize_field(srcXSurfZ);

  initialize_field(srcYVol);
  initialize_field(srcYSurfX);
  initialize_field(srcYSurfY);
  initialize_field(srcYSurfZ);

  initialize_field(srcZVol);
  initialize_field(srcZSurfX);
  initialize_field(srcZSurfY);
  initialize_field(srcZSurfZ);

  OperatorDatabase opdb;
  build_stencils( npts[0], npts[1], npts[2], length, length, length, opdb );

  try{
      //Stencil2 tests:

      status( test_basic_stencils(npts, opdb, srcSVol, srcSSurfX, srcSSurfY, srcSSurfZ, refSVol, refSSurfX, refSSurfY, refSSurfZ, testSVol, testSSurfX, testSSurfY, testSSurfZ), "SVol operators");
      if( npts[0]>1 ) { status( test_basic_stencils(npts, opdb, srcXVol, srcXSurfX, srcXSurfY, srcXSurfZ, refXVol, refXSurfX, refXSurfY, refXSurfZ, testXVol, testXSurfX, testXSurfY, testXSurfZ), "XVol operators"); };
      if( npts[1]>1 ) { status( test_basic_stencils(npts, opdb, srcYVol, srcYSurfX, srcYSurfY, srcYSurfZ, refYVol, refYSurfX, refYSurfY, refYSurfZ, testYVol, testYSurfX, testYSurfY, testYSurfZ), "YVol operators"); };
      if( npts[2]>1 ) { status( test_basic_stencils(npts, opdb, srcZVol, srcZSurfX, srcZSurfY, srcZSurfZ, refZVol, refZSurfX, refZSurfY, refZSurfZ, testZVol, testZSurfX, testZSurfY, testZSurfZ), "ZVol operators"); };


      if( npts[0]>1 & npts[1]>1 ) status( test_stencil2<Interpolant, XVolField, YSurfXField>(opdb, srcXVol, refYSurfX, testYSurfX), "Interpolant XVolField -> YSurfXField (2)" );
      if( npts[0]>1 & npts[1]>1 ) status( test_stencil2<Gradient,    XVolField, YSurfXField>(opdb, srcXVol, refYSurfX, testYSurfX), "Gradient    XVolField -> YSurfXField (2)" );

      if( npts[0]>1 & npts[2]>1 ) status( test_stencil2<Interpolant, XVolField, ZSurfXField>(opdb, srcXVol, refZSurfX, testZSurfX), "Interpolant XVolField -> ZSurfXField (2)" );
      if( npts[0]>1 & npts[2]>1 ) status( test_stencil2<Gradient,    XVolField, ZSurfXField>(opdb, srcXVol, refZSurfX, testZSurfX), "Gradient    XVolField -> ZSurfXField (2)" );

      if( npts[1]>1 & npts[0]>1 ) status( test_stencil2<Interpolant, YVolField, XSurfYField>(opdb, srcYVol, refXSurfY, testXSurfY), "Interpolant YVolField -> XSurfYField (2)" );
      if( npts[1]>1 & npts[0]>1 ) status( test_stencil2<Gradient,    YVolField, XSurfYField>(opdb, srcYVol, refXSurfY, testXSurfY), "Gradient    YVolField -> XSurfYField (2)" );

      if( npts[1]>1 & npts[2]>1 ) status( test_stencil2<Interpolant, YVolField, ZSurfYField>(opdb, srcYVol, refZSurfY, testZSurfY), "Interpolant YVolField -> ZSurfYField (2)" );
      if( npts[1]>1 & npts[2]>1 ) status( test_stencil2<Gradient,    YVolField, ZSurfYField>(opdb, srcYVol, refZSurfY, testZSurfY), "Gradient    YVolField -> ZSurfYField (2)" );

      if( npts[2]>1 & npts[0]>1 ) status( test_stencil2<Interpolant, ZVolField, XSurfZField>(opdb, srcZVol, refXSurfZ, testXSurfZ), "Interpolant ZVolField -> XSurfZField (2)" );
      if( npts[2]>1 & npts[0]>1 ) status( test_stencil2<Gradient,    ZVolField, XSurfZField>(opdb, srcZVol, refXSurfZ, testXSurfZ), "Gradient    ZVolField -> XSurfZField (2)" );

      if( npts[2]>1 & npts[1]>1 ) status( test_stencil2<Interpolant, ZVolField, YSurfZField>(opdb, srcZVol, refYSurfZ, testYSurfZ), "Interpolant ZVolField -> YSurfZField (2)" );
      if( npts[2]>1 & npts[1]>1 ) status( test_stencil2<Gradient,    ZVolField, YSurfZField>(opdb, srcZVol, refYSurfZ, testYSurfZ), "Gradient    ZVolField -> YSurfZField (2)" );

      if( npts[0]>1 ) status( test_stencil2<Interpolant, SVolField, XVolField>(opdb, srcSVol, refXVol, testXVol), "Interpolant SVolField -> XVolField (2)" );
      if( npts[0]>1 ) status( test_stencil2<Gradient,    SVolField, XVolField>(opdb, srcSVol, refXVol, testXVol), "Gradient    SVolField -> XVolField (2)" );

      if( npts[1]>1 ) status( test_stencil2<Interpolant, SVolField, YVolField>(opdb, srcSVol, refYVol, testYVol), "Interpolant SVolField -> YVolField (2)" );
      if( npts[1]>1 ) status( test_stencil2<Gradient,    SVolField, YVolField>(opdb, srcSVol, refYVol, testYVol), "Gradient    SVolField -> YVolField (2)" );

      if( npts[2]>1 ) status( test_stencil2<Interpolant, SVolField, ZVolField>(opdb, srcSVol, refZVol, testZVol), "Interpolant SVolField -> ZVolField (2)" );
      if( npts[2]>1 ) status( test_stencil2<Gradient,    SVolField, ZVolField>(opdb, srcSVol, refZVol, testZVol), "Gradient    SVolField -> ZVolField (2)" );

      if( npts[0]>1 ) status( test_stencil2<Interpolant, XVolField, SVolField>(opdb, srcXVol, refSVol, testSVol), "Interpolant XVolField -> SVolField (2)" );
      if( npts[0]>1 ) status( test_stencil2<Gradient,    XVolField, SVolField>(opdb, srcXVol, refSVol, testSVol), "Gradient    XVolField -> SVolField (2)" );

      if( npts[1]>1 ) status( test_stencil2<Interpolant, YVolField, SVolField>(opdb, srcYVol, refSVol, testSVol), "Interpolant YVolField -> SVolField (2)" );
      if( npts[1]>1 ) status( test_stencil2<Gradient,    YVolField, SVolField>(opdb, srcYVol, refSVol, testSVol), "Gradient    YVolField -> SVolField (2)" );

      if( npts[2]>1 ) status( test_stencil2<Interpolant, ZVolField, SVolField>(opdb, srcZVol, refSVol, testSVol), "Interpolant ZVolField -> SVolField (2)" );
      if( npts[2]>1 ) status( test_stencil2<Gradient,    ZVolField, SVolField>(opdb, srcZVol, refSVol, testSVol), "Gradient    ZVolField -> SVolField (2)" );

      if( npts[0]>1 ) status( test_stencil2<Interpolant, XSurfXField, XVolField>(opdb, srcXSurfX, refXVol, testXVol), "Interpolant XSurfXField -> XVolField (2)" );
      if( npts[1]>1 ) status( test_stencil2<Interpolant, XSurfYField, XVolField>(opdb, srcXSurfY, refXVol, testXVol), "Interpolant XSurfYField -> XVolField (2)" );
      if( npts[2]>1 ) status( test_stencil2<Interpolant, XSurfZField, XVolField>(opdb, srcXSurfZ, refXVol, testXVol), "Interpolant XSurfZField -> XVolField (2)" );

      if( npts[0]>1 ) status( test_stencil2<Interpolant, YSurfXField, YVolField>(opdb, srcYSurfX, refYVol, testYVol), "Interpolant YSurfXField -> YVolField (2)" );
      if( npts[1]>1 ) status( test_stencil2<Interpolant, YSurfYField, YVolField>(opdb, srcYSurfY, refYVol, testYVol), "Interpolant YSurfYField -> YVolField (2)" );
      if( npts[2]>1 ) status( test_stencil2<Interpolant, YSurfZField, YVolField>(opdb, srcYSurfZ, refYVol, testYVol), "Interpolant YSurfZField -> YVolField (2)" );

      if( npts[0]>1 ) status( test_stencil2<Interpolant, ZSurfXField, ZVolField>(opdb, srcZSurfX, refZVol, testZVol), "Interpolant ZSurfXField -> ZVolField (2)" );
      if( npts[1]>1 ) status( test_stencil2<Interpolant, ZSurfYField, ZVolField>(opdb, srcZSurfY, refZVol, testZVol), "Interpolant ZSurfYField -> ZVolField (2)" );
      if( npts[2]>1 ) status( test_stencil2<Interpolant, ZSurfZField, ZVolField>(opdb, srcZSurfZ, refZVol, testZVol), "Interpolant ZSurfZField -> ZVolField (2)" );

      //NullStencil tests:
      //Not yet added

      //Stencil4 tests:

      if( npts[0]>1 & npts[1]>1 ) status( test_stencil4<Interpolant, SVolField, XSurfYField>(opdb, srcSVol, refXSurfY, testXSurfY), "Interpolant SVolField -> XSurfYField (4)" );
      if( npts[0]>1 & npts[2]>1 ) status( test_stencil4<Interpolant, SVolField, XSurfZField>(opdb, srcSVol, refXSurfZ, testXSurfZ), "Interpolant SVolField -> XSurfZField (4)" );

      if( npts[1]>1 & npts[0]>1 ) status( test_stencil4<Interpolant, SVolField, YSurfXField>(opdb, srcSVol, refYSurfX, testYSurfX), "Interpolant SVolField -> YSurfXField (4)" );
      if( npts[1]>1 & npts[2]>1 ) status( test_stencil4<Interpolant, SVolField, YSurfZField>(opdb, srcSVol, refYSurfZ, testYSurfZ), "Interpolant SVolField -> YSurfZField (4)" );

      if( npts[2]>1 & npts[0]>1 ) status( test_stencil4<Interpolant, SVolField, ZSurfXField>(opdb, srcSVol, refZSurfX, testZSurfX), "Interpolant SVolField -> ZSurfXField (4)" );
      if( npts[2]>1 & npts[1]>1 ) status( test_stencil4<Interpolant, SVolField, ZSurfYField>(opdb, srcSVol, refZSurfY, testZSurfY), "Interpolant SVolField -> ZSurfYField (4)" );

      if( npts[0]>1 & npts[1]>1 ) status( test_stencil4<Interpolant, XSurfYField, SVolField>(opdb, srcXSurfY, refSVol, testSVol), "Interpolant XSurfYField -> SVolField (4)" );
      if( npts[0]>1 & npts[2]>1 ) status( test_stencil4<Interpolant, XSurfZField, SVolField>(opdb, srcXSurfZ, refSVol, testSVol), "Interpolant XSurfZField -> SVolField (4)" );

      if( npts[1]>1 & npts[0]>1 ) status( test_stencil4<Interpolant, YSurfXField, SVolField>(opdb, srcYSurfX, refSVol, testSVol), "Interpolant YSurfXField -> SVolField (4)" );
      if( npts[1]>1 & npts[2]>1 ) status( test_stencil4<Interpolant, YSurfZField, SVolField>(opdb, srcYSurfZ, refSVol, testSVol), "Interpolant YSurfZField -> SVolField (4)" );

      if( npts[2]>1 & npts[0]>1 ) status( test_stencil4<Interpolant, ZSurfXField, SVolField>(opdb, srcZSurfX, refSVol, testSVol), "Interpolant ZSurfXField -> SVolField (4)" );
      if( npts[2]>1 & npts[1]>1 ) status( test_stencil4<Interpolant, ZSurfYField, SVolField>(opdb, srcZSurfY, refSVol, testSVol), "Interpolant ZSurfYField -> SVolField (4)" );

      if( npts[0]>1 & npts[1]>1 ) status( test_stencil4<Interpolant, XVolField, YVolField>(opdb, srcXVol, refYVol, testYVol), "Interpolant XVolField -> YVolField (4)" );
      if( npts[0]>1 & npts[2]>1 ) status( test_stencil4<Interpolant, XVolField, ZVolField>(opdb, srcXVol, refZVol, testZVol), "Interpolant XVolField -> ZVolField (4)" );

      if( npts[1]>1 & npts[0]>1 ) status( test_stencil4<Interpolant, YVolField, XVolField>(opdb, srcYVol, refXVol, testXVol), "Interpolant YVolField -> XVolField (4)" );
      if( npts[1]>1 & npts[2]>1 ) status( test_stencil4<Interpolant, YVolField, ZVolField>(opdb, srcYVol, refZVol, testZVol), "Interpolant YVolField -> ZVolField (4)" );

      if( npts[2]>1 & npts[0]>1 ) status( test_stencil4<Interpolant, ZVolField, XVolField>(opdb, srcZVol, refXVol, testXVol), "Interpolant ZVolField -> XVolField (4)" );
      if( npts[2]>1 & npts[1]>1 ) status( test_stencil4<Interpolant, ZVolField, YVolField>(opdb, srcZVol, refYVol, testYVol), "Interpolant ZVolField -> YVolField (4)" );

      //Box filter tests:
      //Not yet added

      //Finite Difference (FDStencil2) tests:

      if( npts[0]>1 ) status( test_fd_stencil2<InterpolantX, SVolField>(opdb, srcSVol, refSVol, testSVol), "InterpolantX SVolField -> SVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<InterpolantY, SVolField>(opdb, srcSVol, refSVol, testSVol), "InterpolantY SVolField -> SVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<InterpolantZ, SVolField>(opdb, srcSVol, refSVol, testSVol), "InterpolantZ SVolField -> SVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<GradientX, SVolField>(opdb, srcSVol, refSVol, testSVol), "GradientX    SVolField -> SVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<GradientY, SVolField>(opdb, srcSVol, refSVol, testSVol), "GradientY    SVolField -> SVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<GradientZ, SVolField>(opdb, srcSVol, refSVol, testSVol), "GradientZ    SVolField -> SVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<InterpolantX, XVolField>(opdb, srcXVol, refXVol, testXVol), "InterpolantX XVolField -> XVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<InterpolantY, XVolField>(opdb, srcXVol, refXVol, testXVol), "InterpolantY XVolField -> XVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<InterpolantZ, XVolField>(opdb, srcXVol, refXVol, testXVol), "InterpolantZ XVolField -> XVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<GradientX, XVolField>(opdb, srcXVol, refXVol, testXVol), "GradientX    XVolField -> XVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<GradientY, XVolField>(opdb, srcXVol, refXVol, testXVol), "GradientY    XVolField -> XVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<GradientZ, XVolField>(opdb, srcXVol, refXVol, testXVol), "GradientZ    XVolField -> XVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<InterpolantX, YVolField>(opdb, srcYVol, refYVol, testYVol), "InterpolantX YVolField -> YVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<InterpolantY, YVolField>(opdb, srcYVol, refYVol, testYVol), "InterpolantY YVolField -> YVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<InterpolantZ, YVolField>(opdb, srcYVol, refYVol, testYVol), "InterpolantZ YVolField -> YVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<GradientX, YVolField>(opdb, srcYVol, refYVol, testYVol), "GradientX    YVolField -> YVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<GradientY, YVolField>(opdb, srcYVol, refYVol, testYVol), "GradientY    YVolField -> YVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<GradientZ, YVolField>(opdb, srcYVol, refYVol, testYVol), "GradientZ    YVolField -> YVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<InterpolantX, ZVolField>(opdb, srcZVol, refZVol, testZVol), "InterpolantX ZVolField -> ZVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<InterpolantY, ZVolField>(opdb, srcZVol, refZVol, testZVol), "InterpolantY ZVolField -> ZVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<InterpolantZ, ZVolField>(opdb, srcZVol, refZVol, testZVol), "InterpolantZ ZVolField -> ZVolField (FD 2)" );

      if( npts[0]>1 ) status( test_fd_stencil2<GradientX, ZVolField>(opdb, srcZVol, refZVol, testZVol), "GradientX    ZVolField -> ZVolField (FD 2)" );
      if( npts[1]>1 ) status( test_fd_stencil2<GradientY, ZVolField>(opdb, srcZVol, refZVol, testZVol), "GradientY    ZVolField -> ZVolField (FD 2)" );
      if( npts[2]>1 ) status( test_fd_stencil2<GradientZ, ZVolField>(opdb, srcZVol, refZVol, testZVol), "GradientZ    ZVolField -> ZVolField (FD 2)" );

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

