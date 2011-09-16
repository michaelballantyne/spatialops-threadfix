#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>

#include <test/TestHelper.h>

#include <string>

using std::string;

#define TEST_EXTRAPTS_WITH_BC( FIELDT, XDIR_RESULT, YDIR_RESULT, ZDIR_RESULT ) \
  status( ExtraPoints<FIELDT, XDIR>::value(true)== XDIR_RESULT, string(#FIELDT)+",  XDIR, with BC" ); \
  status( ExtraPoints<FIELDT, YDIR>::value(true)== YDIR_RESULT, string(#FIELDT)+",  YDIR, with BC" ); \
  status( ExtraPoints<FIELDT, ZDIR>::value(true)== ZDIR_RESULT, string(#FIELDT)+",  ZDIR, with BC" );

#define TEST_EXTRAPTS_WITHOUT_BC( FIELDT, XDIR_RESULT, YDIR_RESULT, ZDIR_RESULT ) \
  status( ExtraPoints<FIELDT, XDIR>::value(false)== XDIR_RESULT, string(#FIELDT)+",  XDIR, without BC" ); \
  status( ExtraPoints<FIELDT, YDIR>::value(false)== YDIR_RESULT, string(#FIELDT)+",  YDIR, without BC" ); \
  status( ExtraPoints<FIELDT, ZDIR>::value(false)== ZDIR_RESULT, string(#FIELDT)+",  ZDIR, without BC" );

#define TEST_STAGGER_IX( FIELDT,  XDIR_RESULT, YDIR_RESULT, ZDIR_RESULT ) \
  status( IndexStagger<FIELDT, XDIR>::value== XDIR_RESULT, "Stagger shift for " + string(#FIELDT) ); \
  status( IndexStagger<FIELDT, YDIR>::value== YDIR_RESULT, "Stagger shift for " + string(#FIELDT) ); \
  status( IndexStagger<FIELDT, ZDIR>::value== ZDIR_RESULT, "Stagger shift for " + string(#FIELDT) );

int main()
{
  TestHelper status(true);

  using namespace SpatialOps;
  using namespace structured;

  TEST_EXTRAPTS_WITH_BC   ( SVolField, 0, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( SVolField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( SSurfXField, 1, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( SSurfXField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( SSurfYField, 0, 1, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( SSurfYField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( SSurfZField, 0, 0, 1 );
  TEST_EXTRAPTS_WITHOUT_BC( SSurfZField, 0, 0, 0 );


  TEST_EXTRAPTS_WITH_BC   ( XVolField, 1, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( XVolField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( XSurfXField, 0, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( XSurfXField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( XSurfYField, 0, 1, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( XSurfYField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( XSurfZField, 0, 0, 1 );
  TEST_EXTRAPTS_WITHOUT_BC( XSurfZField, 0, 0, 0 );


  TEST_EXTRAPTS_WITH_BC   ( YVolField, 0, 1, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( YVolField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( YSurfXField, 1, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( YSurfXField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( YSurfYField, 0, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( YSurfYField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( YSurfZField, 0, 0, 1 );
  TEST_EXTRAPTS_WITHOUT_BC( YSurfZField, 0, 0, 0 );


  TEST_EXTRAPTS_WITH_BC   ( ZVolField, 0, 0, 1 );
  TEST_EXTRAPTS_WITHOUT_BC( ZVolField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( ZSurfXField, 1, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( ZSurfXField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( ZSurfYField, 0, 1, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( ZSurfYField, 0, 0, 0 );
  TEST_EXTRAPTS_WITH_BC   ( ZSurfZField, 0, 0, 0 );
  TEST_EXTRAPTS_WITHOUT_BC( ZSurfZField, 0, 0, 0 );


  TEST_STAGGER_IX( SVolField  , 0, 0, 0 );
  TEST_STAGGER_IX( SSurfXField,-1, 0, 0 );
  TEST_STAGGER_IX( SSurfYField, 0,-1, 0 );
  TEST_STAGGER_IX( SSurfZField, 0, 0,-1 );

  TEST_STAGGER_IX( XVolField  ,-1, 0, 0 );
  TEST_STAGGER_IX( XSurfXField, 0, 0, 0 );
  TEST_STAGGER_IX( XSurfYField, 0,-1, 0 );
  TEST_STAGGER_IX( XSurfZField, 0, 0,-1 );

  TEST_STAGGER_IX( YVolField  , 0,-1, 0 );
  TEST_STAGGER_IX( YSurfXField,-1, 0, 0 );
  TEST_STAGGER_IX( YSurfYField, 0, 0, 0 );
  TEST_STAGGER_IX( YSurfZField, 0, 0,-1 );

  TEST_STAGGER_IX( ZVolField  , 0, 0,-1 );
  TEST_STAGGER_IX( ZSurfXField,-1, 0, 0 );
  TEST_STAGGER_IX( ZSurfYField, 0,-1, 0 );
  TEST_STAGGER_IX( ZSurfZField, 0, 0, 0 );

  if( status.ok() ) return 0;
  return -1;
}
