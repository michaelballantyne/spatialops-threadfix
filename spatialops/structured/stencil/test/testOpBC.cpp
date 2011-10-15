#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

using std::cout;
using std::endl;

#include <spatialops/OperatorDatabase.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/SpatialFieldStore.h>

#include <spatialops/structured/stencil/FVStaggeredBCOp.h>
#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/FVStaggeredBCTools.h>
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
                     const BCSide side,
                     const double bcVal )
{
  using namespace SpatialOps;
  using namespace structured;

  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const OpT& op = *opDB.retrieve_operator<OpT>();

  SpatFldPtr<SrcFieldT > f  = SpatialFieldStore<SrcFieldT >::self().get( get_window_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]) );
  SpatFldPtr<DestFieldT> df = SpatialFieldStore<DestFieldT>::self().get( get_window_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) );

  int icnt=0;
  for( typename SrcFieldT::iterator ifld=f->begin(); ifld!=f->end(); ++ifld,++icnt ) *ifld = icnt;
  *df = 0.0;

  // assign the BC.
  BoundaryConditionOp<OpT,ConstValEval> bc( ijk, side, ConstValEval(bcVal), opDB );

  // evaluate the BC and set it in the field.
  bc(*f);

  // calculate the dest field
  op.apply_to_field( *f, *df );

  // verify that the BC was set properly.
  const int ix = df->window_without_ghost().flat_index(ijk);

  const double abserr = std::abs( (*df)[ix] - bcVal );
  const double relerr = abserr/std::abs(bcVal);

  const bool isOkay = (abserr<2.0e-12 && relerr<2.0e-12);
  if( !isOkay )
    cout << "expected: " << bcVal << endl
         << "found: " << (*df)[ix] << endl
         << "abserr: " << abserr << ", relerr: " << relerr << endl
         << endl;
  return isOkay;
}

//--------------------------------------------------------------------

template<typename OpT>
bool
test_bc_loop( const OperatorDatabase& opDB,
              const std::string opName,
              const IntVec& dim,
              const std::vector<bool>& bcFlag,
              const int bcFaceIndex,
              const BCSide side,
              const double bcVal )
{
  typedef typename OpT:: SrcFieldType::Location::Offset SO;
  typedef typename OpT::DestFieldType::Location::Offset DO;
  typedef typename GetNonzeroDir< typename Subtract<SO,DO>::result >::DirT ActiveDir;

  TestHelper status(false);

  if( IsSameType<ActiveDir,XDIR>::result ){
    TestHelper tmp(false);
    const int i=bcFaceIndex;
    for( int j=0; j<dim[1]; ++j ){
      for( int k=0; k<dim[2]; ++k ){
        tmp( test_bc_helper<OpT>( opDB, dim, bcFlag, IntVec(i,j,k), side, bcVal ) );
      }
    }
    status( tmp.ok(), "X BC " + opName );
  }
  else if( IsSameType<ActiveDir,YDIR>::result ){
    TestHelper tmp(false);
    const int j=bcFaceIndex;
    for( int i=0; i<dim[0]; ++i ){
      for( int k=0; k<dim[2]; ++k ){
        tmp( test_bc_helper<OpT>( opDB, dim, bcFlag, IntVec(i,j,k), side, bcVal ) );
      }
    }
    status( tmp.ok(), "Y BC " + opName );
  }
  else{ // Z-Dir.
    TestHelper tmp(false);
    const int k=bcFaceIndex;
    for( int i=0; i<dim[0]; ++i ){
      for( int j=0; j<dim[1]; ++j ){
        tmp( test_bc_helper<OpT>( opDB, dim, bcFlag, IntVec(i,j,k), side, bcVal ) );
      }
    }
    status( tmp.ok(), "Z BC " + opName );
  }

  return status.ok();
}

//--------------------------------------------------------------------

bool test_bc( const OperatorDatabase& opDB,
              const Grid& g,
              const std::vector<bool>& bcFlag )
{
  using namespace SpatialOps;
  using namespace structured;

  typedef BasicOpTypes<SVolField> Ops;

  const IntVec& dim = g.extent();

  TestHelper status(true);

  if( dim[0]>1 ){
    // X BCs - Left side
    int i=0;
    status( test_bc_loop<Ops::GradX     >( opDB, "GradSVolSSurfX",   dim, bcFlag, i,   MINUS_SIDE, 1.2345 ), "-X GradSVolSSurfX"   );
    status( test_bc_loop<Ops::InterpC2FX>( opDB, "InterpSVolSSurfX", dim, bcFlag, i,   MINUS_SIDE, 123.45 ), "-X InterpSVolSSurfX" );
    status( test_bc_loop<Ops::DivX      >( opDB, "DivSSurfXSVol",    dim, bcFlag, i,   MINUS_SIDE, 123.45 ), "-X DivSSurfXSVol"    );
    cout << endl;

    // X BCs - Right side
    i=dim[0];
    status( test_bc_loop<Ops::GradX     >( opDB, "GradSVolSSurfX",   dim, bcFlag, i,   PLUS_SIDE, 1.2345 ), "+X GradSVolSSurfX"   );
    status( test_bc_loop<Ops::InterpC2FX>( opDB, "InterpSVolSSurfX", dim, bcFlag, i,   PLUS_SIDE, 123.45 ), "+X InterpSVolSSurfX" );
    status( test_bc_loop<Ops::DivX      >( opDB, "DivSSurfXSVol",    dim, bcFlag, i-1, PLUS_SIDE, 123.45 ), "+X DivSSurfXSVol"    );
    cout << endl;
  }

  if( dim[1]>1 ){
    // Y BCs - Left side
    int j=0;
    status( test_bc_loop<Ops::GradY     >( opDB, "GradSVolSSurfY",   dim, bcFlag, j,   MINUS_SIDE, 1.23456 ), "-Y GradSVolSSurfY"   );
    status( test_bc_loop<Ops::InterpC2FY>( opDB, "InterpSVolSSurfY", dim, bcFlag, j,   MINUS_SIDE, 123.456 ), "-Y InterpSVolSSurfY" );
    status( test_bc_loop<Ops::DivY      >( opDB, "DivSSurfYSVol",    dim, bcFlag, j,   MINUS_SIDE, 123.45  ), "-Y DivSSurfYSVol"    );
    cout << endl;

    // Y BCs - Right side
    j=dim[1];
    status( test_bc_loop<Ops::GradY     >( opDB, "GradSVolSSurfY",   dim, bcFlag, j,   PLUS_SIDE, 6.54321 ), "+Y GradSVolSSurfY"   );
    status( test_bc_loop<Ops::InterpC2FY>( opDB, "InterpSVolSSurfY", dim, bcFlag, j,   PLUS_SIDE, 123.456 ), "+Y InterpSVolSSurfY" );
    status( test_bc_loop<Ops::DivY      >( opDB, "DivSSurfYSVol",    dim, bcFlag, j-1, PLUS_SIDE, 123.45  ), "+Y DivSSurfYSVol"    );
    cout << endl;
  }

  if( dim[2]>1 ){
    // Z BCs - Left side
    int k=0;
    status( test_bc_loop<Ops::GradZ     >( opDB, "GradSVolSSurfZ",   dim, bcFlag, k,   MINUS_SIDE, 1.2345  ), "-Z GradSVolSSurfZ"   );
    status( test_bc_loop<Ops::InterpC2FZ>( opDB, "InterpSVolSSurfZ", dim, bcFlag, k,   MINUS_SIDE, 123.456 ), "-Z InterpSVolSSurfZ" );
    status( test_bc_loop<Ops::DivZ      >( opDB, "DivSSurfZSVol",    dim, bcFlag, k,   MINUS_SIDE, 123.45  ), "-Z DivSSurfZSVol"    );
    cout << endl;

    // Z BCs - Right side
    k=dim[2];
    status( test_bc_loop<Ops::GradZ     >( opDB, "GradSVolSSurfZ",   dim, bcFlag, k,   PLUS_SIDE, 1.2345  ), "+Z GradSVolSSurfZ"   );
    status( test_bc_loop<Ops::InterpC2FZ>( opDB, "InterpSVolSSurfZ", dim, bcFlag, k,   PLUS_SIDE, 123.456 ), "+Z InterpSVolSSurfZ" );
    status( test_bc_loop<Ops::DivZ      >( opDB, "DivSSurfZSVol",    dim, bcFlag, k-1, PLUS_SIDE, 123.45  ), "+Z DivSSurfZSVol"    );
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

  status( test_driver( IntVec(3,1 ,1 ) ), "Mesh: (10,1,1)\n" );
  status( test_driver( IntVec(1 ,10,1 ) ), "Mesh: (1,10,1)\n" );
  status( test_driver( IntVec(1 ,1 ,10) ), "Mesh: (1,1,10)\n" );

  status( test_driver( IntVec(10,10,1 ) ), "Mesh: (10,10,1)\n" );
  status( test_driver( IntVec(10,1 ,10) ), "Mesh: (10,1,10)\n" );
  status( test_driver( IntVec(1 ,10,10) ), "Mesh: (1,10,10)\n" );

  status( test_driver( IntVec(10,10,10) ), "Mesh: (10,10,10)\n" );

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
