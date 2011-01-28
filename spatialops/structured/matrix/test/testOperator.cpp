#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/matrix/FVStaggeredOperatorTypes.h>

#include <spatialops/structured/matrix/FVStaggeredInterpolant.h>
#include <spatialops/structured/matrix/FVStaggeredGradient.h>
#include <spatialops/structured/matrix/FVStaggeredDivergence.h>
#include <spatialops/structured/matrix/FVStaggeredScratch.h>
#include <spatialops/structured/matrix/FVTopHatFilter.h>
#include <spatialops/structured/matrix/FVRestrictOp.h>

#include <spatialops/structured/FVTools.h>

#include <test/TestHelper.h>
#include "buildOps.h"

#include <iostream>

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;
using std::vector;

template< typename FieldLoc >
struct FieldRelate{};

template<> struct FieldRelate<SVol>{
  typedef SVolField   Vol;
  typedef SSurfXField XSurf;
  typedef SSurfYField YSurf;
  typedef SSurfZField ZSurf;
};


template<> struct FieldRelate<XVol>{
  typedef XVolField   Vol;
  typedef XSurfXField XSurf;
  typedef XSurfYField YSurf;
  typedef XSurfZField ZSurf;
};


template<> struct FieldRelate<YVol>{
  typedef YVolField   Vol;
  typedef YSurfXField XSurf;
  typedef YSurfYField YSurf;
  typedef YSurfZField ZSurf;
};


template<> struct FieldRelate<ZVol>{
  typedef ZVolField   Vol;
  typedef ZSurfXField XSurf;
  typedef ZSurfYField YSurf;
  typedef ZSurfZField ZSurf;
};

template<typename VolT>
bool test_ops( const IntVec dim,
               const vector<double> spc,
               const vector<bool> bcPlus )
{
  OperatorDatabase opDB;
  build_ops( dim, spc, bcPlus, opDB );

  typedef FieldRelate<VolT> Fields;
  typedef typename Fields::Vol   Vol;
  typedef typename Fields::XSurf XSurf;
  typedef typename Fields::YSurf YSurf;
  typedef typename Fields::ZSurf ZSurf;

  try{
    Vol   v( get_window_with_ghost<Vol  >( dim,bcPlus[0],bcPlus[1],bcPlus[2]), NULL );
    XSurf x( get_window_with_ghost<XSurf>( dim,bcPlus[0],bcPlus[1],bcPlus[2]), NULL );
    YSurf y( get_window_with_ghost<YSurf>( dim,bcPlus[0],bcPlus[1],bcPlus[2]), NULL );
    ZSurf z( get_window_with_ghost<ZSurf>( dim,bcPlus[0],bcPlus[1],bcPlus[2]), NULL );

    v = 1.0;

    typedef SpatialOperator< LinAlg, Gradient, Vol, XSurf > GradX;
    typedef SpatialOperator< LinAlg, Gradient, Vol, YSurf > GradY;
    typedef SpatialOperator< LinAlg, Gradient, Vol, ZSurf > GradZ;

    const GradX* xgrad = opDB.retrieve_operator<GradX>();  xgrad->apply_to_field(v,x);
    const GradY* ygrad = opDB.retrieve_operator<GradY>();  ygrad->apply_to_field(v,y);
    const GradZ* zgrad = opDB.retrieve_operator<GradZ>();  zgrad->apply_to_field(v,z);

    typedef SpatialOperator< LinAlg, Divergence, XSurf, Vol > DivX;
    typedef SpatialOperator< LinAlg, Divergence, YSurf, Vol > DivY;
    typedef SpatialOperator< LinAlg, Divergence, ZSurf, Vol > DivZ;

    const DivX* xdiv = opDB.retrieve_operator<DivX>();  xdiv->apply_to_field(x,v);
    const DivY* ydiv = opDB.retrieve_operator<DivY>();  ydiv->apply_to_field(y,v);
    const DivZ* zdiv = opDB.retrieve_operator<DivZ>();  zdiv->apply_to_field(z,v);

    typedef SpatialOperator< LinAlg, Interpolant, Vol, XSurf > InterpVX;
    typedef SpatialOperator< LinAlg, Interpolant, Vol, YSurf > InterpVY;
    typedef SpatialOperator< LinAlg, Interpolant, Vol, ZSurf > InterpVZ;

    const InterpVX* vxinterp = opDB.retrieve_operator<InterpVX>();  vxinterp->apply_to_field(v,x);
    const InterpVY* vyinterp = opDB.retrieve_operator<InterpVY>();  vyinterp->apply_to_field(v,y);
    const InterpVZ* vzinterp = opDB.retrieve_operator<InterpVZ>();  vzinterp->apply_to_field(v,z);
  }
  catch( std::exception& e ){
    cout << e.what() << endl;
    return false;
  }
  return true;
}



int main()
{
  const IntVec dim(4,4,4);
  const vector<double> spc(3,1.0);
  vector<bool> bcPlus(3,true);

  TestHelper status(true);

  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [T,T,T]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [T,T,T]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [T,T,T]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [T,T,T]" );

  bcPlus[0] = bcPlus[1] = bcPlus[2] = false;
  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [F,F,F]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [F,F,F]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [F,F,F]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [F,F,F]" );

  bcPlus[1] = bcPlus[2] = true;
  bcPlus[0]=false;
  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [F,T,T]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [F,T,T]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [F,T,T]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [F,T,T]" );

  bcPlus[1]=false;
  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [F,F,T]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [F,F,T]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [F,F,T]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [F,F,T]" );

  bcPlus[2]=false;
  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [F,F,F]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [F,F,F]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [F,F,F]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [F,F,F]" );

  bcPlus[0]=true;
  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [T,F,F]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [T,F,F]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [T,F,F]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [T,F,F]" );

  bcPlus[1]=true;
  status( test_ops<SVol>( dim, spc, bcPlus ), "SVol [T,T,F]" );
  status( test_ops<XVol>( dim, spc, bcPlus ), "XVol [T,T,F]" );
  status( test_ops<YVol>( dim, spc, bcPlus ), "YVol [T,T,F]" );
  status( test_ops<ZVol>( dim, spc, bcPlus ), "ZVol [T,T,F]" );

  if( status.ok() ) return 0;
  return -1;
}
