#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
using namespace std;


#include <spatialops/OperatorDatabase.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/SpatialFieldStore.h>

#include <spatialops/structured/stencil/FVStaggeredBCOp.h>
#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/stencil/StencilBuilder.h>

#include <spatialops/structured/Grid.h>
#include <test/TestHelper.h>

using namespace SpatialOps;
using namespace structured;

//--------------------------------------------------------------------

template<typename OpT>
bool test_bc_helper( const OperatorDatabase& opDB,
                     const IntVec&dim,
                     const std::vector<bool>& bcFlag,
                     const IntVec ijk,
                     const double bcVal,
                     const BCSide side )
{
  using namespace SpatialOps;
  using namespace structured;

  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const OpT& op = *opDB.retrieve_operator<OpT>();
  //  op.write_matlab("op");

  SpatFldPtr<SrcFieldT > f  = SpatialFieldStore<SrcFieldT >::self().get( get_window_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]) );
  SpatFldPtr<DestFieldT> df = SpatialFieldStore<DestFieldT>::self().get( get_window_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) );

  int icnt=0;
  for( typename SrcFieldT::iterator ifld=f->begin(); ifld!=f->end(); ++ifld,++icnt ) *ifld = icnt;

  // assign the BC.
  BoundaryConditionOp<OpT,ConstValEval> bc( dim,
                                            bcFlag[0], bcFlag[1], bcFlag[2],
                                            ijk,
                                            side,
                                            ConstValEval(bcVal),
                                            opDB );

  // evaluate the BC and set it in the field.
  bc(*f);

  // calculate the dest field
  op.apply_to_field( *f, *df );

  // verify that the BC was set properly.
  const int ix = df->window_without_ghost().flat_index(ijk);

  const double abserr = abs( (*df)[ix] - bcVal );
  const double relerr = abserr/abs(bcVal);

  bool isOkay = (abserr<1.0e-9 && relerr<1.0e-9);
//   if( !isOkay ) cout << "  " << abserr << "  " << relerr << "     "
//                   << df[ix] << "  " << bcVal << endl;
  return isOkay;
}

//--------------------------------------------------------------------

template<typename OpT>
bool
test_bc_loop( const OperatorDatabase& opDB,
              const std::string opName,
              const IntVec& dim,
              const std::vector<bool>& bcFlag,
              const BCSide side,
              const int bcFaceIndex,
              const double bcVal )
{
  TestHelper status(false);

  switch(side){

  case X_MINUS_SIDE:
  case X_PLUS_SIDE: {
    TestHelper tmp(false);
    const int i=bcFaceIndex;
    for( int j=0; j<dim[1]; ++j ){
      for( int k=0; k<dim[2]; ++k ){
        tmp( test_bc_helper<OpT>( opDB, dim, bcFlag, IntVec(i,j,k), bcVal, side ) );
      }
    }
    status( tmp.ok(), "X BC " + opName );
    break;
  }

  case Y_MINUS_SIDE:
  case Y_PLUS_SIDE: {
    TestHelper tmp(false);
    const int j=bcFaceIndex;
    for( int i=0; i<dim[0]; ++i ){
      for( int k=0; k<dim[2]; ++k ){
        tmp( test_bc_helper<OpT>( opDB, dim, bcFlag, IntVec(i,j,k), bcVal, side ) );
      }
    }
    status( tmp.ok(), "Y BC " + opName );
    break;
  }

  case Z_MINUS_SIDE:
  case Z_PLUS_SIDE: {
    TestHelper tmp(false);
    const int k=bcFaceIndex;
    for( int i=0; i<dim[0]; ++i ){
      for( int j=0; j<dim[1]; ++j ){
        tmp( test_bc_helper<OpT>( opDB, dim, bcFlag, IntVec(i,j,k), bcVal, side ) );
      }
    }
    status( tmp.ok(), "Z BC " + opName );
    break;
  }

  case NO_SHIFT:
    assert(1);

  } // switch(side)

  return status.ok();
}

//--------------------------------------------------------------------

bool test_bc( const OperatorDatabase& opDB,
              const Grid& g,
              const std::vector<bool>& bcFlag )
{
  using namespace SpatialOps;
  using namespace structured;

  typedef BasicOpTypes<SVolField> SOps;
  typedef BasicOpTypes<XVolField> XOps;
  typedef BasicOpTypes<YVolField> YOps;
  typedef BasicOpTypes<ZVolField> ZOps;

  const IntVec& dim = g.extent();

  TestHelper status(true);

  if( dim[0]>1 ){
    // X BCs - Left side
    int i=0;

    status( test_bc_loop<SOps::GradX     >( opDB, "GradSVolSSurfX",   dim, bcFlag, X_MINUS_SIDE, i, 1.2345 ), "-X GradSVolSSurfX" );
    status( test_bc_loop<SOps::InterpC2FX>( opDB, "InterpSVolSSurfX", dim, bcFlag, X_MINUS_SIDE, i, 123.45 ), "-X InterpSVolSSurfX" );

    status( test_bc_loop<XOps::GradX     >( opDB, "GradXVolXSurfX",   dim, bcFlag, X_MINUS_SIDE, i, 5.4321  ), "-X GradXVolXSurfX" );
    status( test_bc_loop<XOps::InterpC2FX>( opDB, "InterpXVolXSurfX", dim, bcFlag, X_MINUS_SIDE, i, 54.321  ), "-X InterpXVolXSurfX" );

    if( dim[1]>1 ){
      status( test_bc_loop<YOps::GradX     >( opDB, "GradYVolYSurfX",   dim, bcFlag, X_MINUS_SIDE, i, 5.4321 ), "-X GradYVolYSurfX" );
      status( test_bc_loop<YOps::InterpC2FX>( opDB, "InterpYVolXSurfX", dim, bcFlag, X_MINUS_SIDE, i, 54.321 ), "-X InterpYVolXSurfX" );
    }
    if( dim[2]>1 ){
      status( test_bc_loop<ZOps::GradX     >( opDB, "GradZVolZSurfX",   dim, bcFlag, X_MINUS_SIDE, i, 5.4321 ), "-X GradZVolZSurfX" );
      status( test_bc_loop<ZOps::InterpC2FX>( opDB, "InterpZVolZSurfX", dim, bcFlag, X_MINUS_SIDE, i, 54.321 ), "-X InterpZVolZSurfX" );
    }

    // X BCs - Right side
    cout << endl;
    i=dim[0]-1;

    status( test_bc_loop<SOps::GradX     >( opDB, "GradSVolSSurfX",   dim, bcFlag, X_PLUS_SIDE, i, 1.2345 ), "+X GradSVolSSurfX" );
    status( test_bc_loop<SOps::InterpC2FX>( opDB, "InterpSVolSSurfX", dim, bcFlag, X_PLUS_SIDE, i, 123.45 ), "+X InterpSVolSSurfX" );

    status( test_bc_loop<XOps::GradX     >( opDB, "GradXVolXSurfX",   dim, bcFlag, X_PLUS_SIDE, i, 5.4321 ), "+X GradXVolXSurfX" );
    status( test_bc_loop<XOps::InterpC2FX>( opDB, "InterpXVolXSurfX", dim, bcFlag, X_PLUS_SIDE, i, 54.321 ), "+X InterpXVolXSurfX" );
    if( dim[1]>1 ){
      status( test_bc_loop<YOps::GradX     >( opDB, "GradYVolYSurfX",   dim, bcFlag, X_PLUS_SIDE, i, 5.4321 ), "+X GradYVolYSurfX" );
      status( test_bc_loop<YOps::InterpC2FX>( opDB, "InterpYVolXSurfX", dim, bcFlag, X_PLUS_SIDE, i, 54.321 ), "+X InterpYVolXSurfX" );
    }
    if( dim[2]>1 ){
      status( test_bc_loop<ZOps::GradX     >( opDB, "GradZVolZSurfX",   dim, bcFlag, X_PLUS_SIDE, i, 5.4321 ), "+X GradZVolZSurfX" );
      status( test_bc_loop<ZOps::InterpC2FX>( opDB, "InterpZVolZSurfX", dim, bcFlag, X_PLUS_SIDE, i, 54.321 ), "+X InterpZVolZSurfX" );
    }
    cout << endl;
  }

  if( dim[1]>1 ){
    // Y BCs - Left side
    int j=0;
    status( test_bc_loop<SOps::GradY     >( opDB, "GradSVolSSurfY",   dim, bcFlag, Y_MINUS_SIDE, j, 1.23456 ), "-Y GradSVolSSurfY"   );
    status( test_bc_loop<SOps::InterpC2FY>( opDB, "InterpSVolSSurfY", dim, bcFlag, Y_MINUS_SIDE, j, 123.456 ), "-Y InterpSVolSSurfY" );

    if( dim[0]>1 ){
      status( test_bc_loop<XOps::GradY     >( opDB, "GradXVolXSurfY",   dim, bcFlag, Y_MINUS_SIDE, j, 1.23456 ), "-Y GradXVolXSurfY"   );
      status( test_bc_loop<XOps::InterpC2FY>( opDB, "InterpXVolXSurfY", dim, bcFlag, Y_MINUS_SIDE, j, 123.456 ), "-Y InterpXVolXSurfY" );
    }

    status( test_bc_loop<YOps::GradY     >( opDB, "GradYVolYSurfY",   dim, bcFlag, Y_MINUS_SIDE, j, 1.23456 ), "-Y GradYVolYSurfY"   );
    status( test_bc_loop<YOps::InterpC2FY>( opDB, "InterpYVolYSurfY", dim, bcFlag, Y_MINUS_SIDE, j, 123.456 ), "-Y InterpYVolYSurfY" );

    if( dim[2]>1 ){
      status( test_bc_loop<ZOps::GradY     >( opDB, "GradZVolZSurfY",   dim, bcFlag, Y_MINUS_SIDE, j, 1.23456 ), "-Y GradZVolZSurfY"   );
      status( test_bc_loop<ZOps::InterpC2FY>( opDB, "InterpZVolZSurfY", dim, bcFlag, Y_MINUS_SIDE, j, 123.456 ), "-Y InterpZVolZSurfY" );
    }

    cout << endl;

    // Y BCs - Right side
    j=dim[1]-1;
    status( test_bc_loop<SOps::GradY     >( opDB, "GradSVolSSurfY",   dim, bcFlag, Y_PLUS_SIDE, j, 6.54321 ), "+Y GradSVolSSurfY"   );
    status( test_bc_loop<SOps::InterpC2FY>( opDB, "InterpSVolSSurfY", dim, bcFlag, Y_PLUS_SIDE, j, 123.456 ), "+Y InterpSVolSSurfY" );

    if( dim[0]>1 ){
      status( test_bc_loop<XOps::GradY     >( opDB, "GradXVolXSurfY",   dim, bcFlag, Y_PLUS_SIDE, j, 1.23456 ), "+Y GradXVolXSurfY" );
      status( test_bc_loop<YOps::InterpC2FY>( opDB, "InterpXVolXSurfY", dim, bcFlag, Y_PLUS_SIDE, j, 123.456 ), "+Y InterpXVolXSurfY" );
    }

    status( test_bc_loop<YOps::GradY     >( opDB, "GradYVolYSurfY",   dim, bcFlag, Y_PLUS_SIDE, j, 1.23456 ), "+Y GradYVolYSurfY" );
    status( test_bc_loop<YOps::InterpC2FY>( opDB, "InterpYVolYSurfY", dim, bcFlag, Y_PLUS_SIDE, j, 123.456 ), "+Y InterpYVolYSurfY" );

    if( dim[2]>1 ){
      status( test_bc_loop<ZOps::GradY     >( opDB, "GradZVolZSurfY",   dim, bcFlag, Y_PLUS_SIDE, j, 1.23456 ), "+Y GradZVolZSurfY" );
      status( test_bc_loop<ZOps::InterpC2FY>( opDB, "InterpZVolZSurfY", dim, bcFlag, Y_PLUS_SIDE, j, 123.456 ), "+Y InterpZVolZSurfY" );
    }

    cout << endl;
  }

  if( dim[2]>1 ){
    // Z BCs - Left side
    int k=0;
    status( test_bc_loop<SOps::GradZ     >( opDB, "GradSVolSSurfZ",   dim, bcFlag, Z_MINUS_SIDE, k, 1.2345 ), "-Z GradSVolSSurfZ" );
    status( test_bc_loop<SOps::InterpC2FZ>( opDB, "InterpSVolSSurfZ", dim, bcFlag, Z_MINUS_SIDE, k, 123.456 ), "-Z InterpSVolSSurfZ" );

    if( dim[0]>1 ){
      status( test_bc_loop<XOps::GradZ     >( opDB, "GradXVolXSurfZ",   dim, bcFlag, Z_MINUS_SIDE, k, 1.2345 ), "-Z GradXVolXSurfZ" );
      status( test_bc_loop<XOps::InterpC2FZ>( opDB, "InterpXVolXSurfZ", dim, bcFlag, Z_MINUS_SIDE, k, 123.456 ), "-Z InterpXVolXSurfZ" );
    }                                      
    if( dim[1]>1 ){                        
      status( test_bc_loop<YOps::GradZ     >( opDB, "GradYVolYSurfZ",   dim, bcFlag, Z_MINUS_SIDE, k, 1.2345 ), "-Z GradYVolYSurfZ" );
      status( test_bc_loop<YOps::InterpC2FZ>( opDB, "InterpYVolYSurfZ", dim, bcFlag, Z_MINUS_SIDE, k, 123.456 ), "-Z InterpYVolYSurfZ" );
    }
    status( test_bc_loop<ZOps::GradZ     >( opDB, "GradZVolZSurfZ",   dim, bcFlag, Z_MINUS_SIDE, k, 1.2345 ), "-Z GradZVolZSurfZ" );
    status( test_bc_loop<ZOps::InterpC2FZ>( opDB, "InterpZVolZSurfZ", dim, bcFlag, Z_MINUS_SIDE, k, 123.456 ), "-Z InterpZVolZSurfZ" );

    cout << endl;

    // Z BCs - Right side
    k=dim[2]-1;
    status( test_bc_loop<SOps::GradZ     >( opDB, "GradSVolSSurfZ",   dim, bcFlag, Z_PLUS_SIDE, k, 1.2345 ), "+Z GradSVolSSurfZ" );
    status( test_bc_loop<SOps::InterpC2FZ>( opDB, "InterpSVolSSurfZ", dim, bcFlag, Z_PLUS_SIDE, k, 123.456 ), "+Z InterpSVolSSurfZ" );

    if( dim[0]>1 ){
      status( test_bc_loop<XOps::GradZ     >( opDB, "GradXVolXSurfZ",   dim, bcFlag, Z_PLUS_SIDE, k, 1.2345 ), "+Z GradXVolXSurfZ" );
      status( test_bc_loop<XOps::InterpC2FZ>( opDB, "InterpXVolXSurfZ", dim, bcFlag, Z_PLUS_SIDE, k, 123.456 ), "+Z InterpXVolXSurfZ" );
    }                                      
    if( dim[1]>1 ){                        
      status( test_bc_loop<YOps::GradZ     >( opDB, "GradYVolYSurfZ",   dim, bcFlag, Z_PLUS_SIDE, k, 1.2345 ), "+Z GradYVolYSurfZ" );
      status( test_bc_loop<YOps::InterpC2FZ>( opDB, "InterpYVolYSurfZ", dim, bcFlag, Z_PLUS_SIDE, k, 123.456 ), "+Z InterpYVolYSurfZ" );
    }
    status( test_bc_loop<ZOps::GradZ     >( opDB, "GradZVolZSurfZ",   dim, bcFlag, Z_PLUS_SIDE, k, 1.2345 ), "+Z GradZVolZSurfZ" );
    status( test_bc_loop<ZOps::InterpC2FZ>( opDB, "InterpZVolZSurfZ", dim, bcFlag, Z_PLUS_SIDE, k, 123.456 ), "+Z InterpZVolZSurfZ" );
    cout << endl;
  }
 
  return status.ok();
}

//--------------------------------------------------------------------

bool test_driver( const IntVec& dim )
{
  std::vector<double> length(3,1);
  std::vector<bool> bcFlag(3,true);

  OperatorDatabase opDB;
  build_stencils( dim[0], dim[1], dim[2], length[0], length[1], length[2], opDB );
  const Grid grid( dim, length );

  return test_bc( opDB, grid, bcFlag );
}

//--------------------------------------------------------------------

int main()
{
  TestHelper status(true);

  status( test_driver( IntVec(10,1 ,1 ) ), "Mesh: (10,1,1)\n" );
  status( test_driver( IntVec(1 ,10,1 ) ), "Mesh: (1,10,1)\n" );
  status( test_driver( IntVec(1 ,1 ,10) ), "Mesh: (1,1,10)\n" );

  status( test_driver( IntVec(10,10,1 ) ), "Mesh: (10,10,1)\n" );
  status( test_driver( IntVec(10,1 ,10) ), "Mesh: (10,1,10)\n" );
  status( test_driver( IntVec(1 ,10,10) ), "Mesh: (1,10,10)\n" );

# ifdef NDEBUG  // too slow in debug mode - only run in release
  status( test_driver( IntVec(10,10,10) ), "Mesh: (10,10,10)\n" );
# endif

    cout << endl << "----------" << endl
         << "BC Op Test: ";
  if( status.ok() ){
    cout << "PASS" << endl << "----------" << endl;
    return 0;
  }
  else{
    cout << "FAIL" << endl << "----------" << endl;
  }
  return -1;
}
