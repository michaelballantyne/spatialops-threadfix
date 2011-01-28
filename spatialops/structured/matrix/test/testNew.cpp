#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
using namespace std;

#include <spatialops/OperatorDatabase.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>

#include <spatialops/structured/matrix/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/matrix/FVStaggeredInterpolant.h>
#include <spatialops/structured/matrix/FVStaggeredGradient.h>
#include <spatialops/structured/matrix/FVStaggeredDivergence.h>
#include <spatialops/structured/matrix/FVStaggeredScratch.h>
#include <spatialops/structured/matrix/FVTopHatFilter.h>
#include <spatialops/structured/matrix/FVRestrictOp.h>

#ifdef LINALG_TRILINOS
# include <spatialops/LinearSystem.h>
#endif

#include "Grid.h"
#include "buildOps.h"
#include <test/Functions.h>

using namespace SpatialOps;
using namespace structured;

OperatorDatabase opDB;  // jcs hack global variable

//--------------------------------------------------------------------

template<typename FieldT>
void compare( const FieldT& f1, const FieldT& f2,
              double& maxAbsErr, double& maxRelErr,
              double& avgAbsErr, double& avgRelErr )
{
  static const double SMALL = numeric_limits<double>::epsilon() * 10.0;

  maxAbsErr = 0.0;
  maxRelErr = 0.0;
  avgAbsErr = 0.0;
  avgRelErr = 0.0;

  int npts = 0;

  typename FieldT::const_interior_iterator if1 = f1.interior_begin();
  typename FieldT::const_interior_iterator if2 = f2.interior_begin();
  const typename FieldT::const_interior_iterator if1e= f1.interior_end();

  for( ; if1!=if1e; ++if1, ++if2 ){
    const double absErr = std::abs( *if1 - *if2 );
    const double relErr = absErr / (std::abs(*if1)+SMALL);
    maxAbsErr = std::max( maxAbsErr, absErr );
    maxRelErr = std::max( maxRelErr, relErr );
    avgRelErr += relErr;
    avgAbsErr += absErr;
    ++npts;
  }
  avgRelErr /= npts;
  avgAbsErr /= npts;
}

//--------------------------------------------------------------------

template<typename FieldT>
void report_errors( const FieldT& phiExact, const FieldT& phi )
{
  double absErr, relErr, avgAbsErr, avgRelErr;
  compare( phiExact, phi, absErr, relErr, avgAbsErr, avgRelErr );

  cout << scientific << setprecision(5)
       << setw(12) << absErr    << " |"
       << setw(12) << relErr    << " |"
       << setw(12) << avgAbsErr << " |"
       << setw(12) << avgRelErr << " |"
       << endl;
}

//--------------------------------------------------------------------

template<typename OpT, typename FuncType1, typename FuncType2>
void test_grad_op( const Grid& grid,
                   const FuncType1& funcPhi,
                   const FuncType2& funcFPhi,
                   const std::vector<bool>& bcFlag,
                   const int dir )
{
  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const IntVec& dim = grid.extent();

  if( get_ntot_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]) == 1 ||
      get_ntot_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) == 1 ) return;

  const MemoryWindow srcWindow( get_window_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]));
  const MemoryWindow dstWindow( get_window_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]));

  SrcFieldT  phi      ( srcWindow, NULL );
  DestFieldT fphi     ( dstWindow, NULL );
  DestFieldT fphiExact( dstWindow, NULL );

  funcPhi.evaluate( phi );

  switch( dir )
  {
    case XDIR::value : funcFPhi.dx( fphiExact );  break;
    case YDIR::value : funcFPhi.dy( fphiExact );  break;
    case ZDIR::value : funcFPhi.dz( fphiExact );  break;
    default: assert(0);
  }

  const OpT* const op = opDB.retrieve_operator<OpT>();

  op->apply_to_field( phi, fphi );

  report_errors( fphiExact, fphi );
}
//--------------------------------------------------------------------

template<typename OpT, typename FuncType1, typename FuncType2>
void test_interp_op( const Grid& grid,
                     const FuncType1& funcPhi,
                     const FuncType2& funcFPhi,
                     const std::vector<bool>& bcFlag )
{
  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const IntVec& dim = grid.extent();

  if( get_ntot_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]) == 1 ||
      get_ntot_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) == 1 ) return;

  const MemoryWindow srcWindow( get_window_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]));
  const MemoryWindow dstWindow( get_window_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]));

  SrcFieldT  phi      ( srcWindow, NULL );
  DestFieldT fphi     ( dstWindow, NULL );
  DestFieldT fphiExact( dstWindow, NULL );

  funcPhi.evaluate( phi );
  funcFPhi.evaluate( fphiExact );

  const OpT* const op = opDB.retrieve_operator<OpT>();
  op->apply_to_field( phi, fphi );

  report_errors( fphiExact, fphi );
}

//--------------------------------------------------------------------

template<typename OpT, typename FuncType1, typename FuncType2>
void test_div_op( const Grid& grid,
                  const FuncType1& funcPhi,
                  const FuncType2& funcFPhi,
                  const std::vector<bool>& bcFlag )
{
  typedef typename OpT::SrcFieldType  SrcFieldT;
  typedef typename OpT::DestFieldType DestFieldT;

  const IntVec& dim = grid.extent();

  if( get_ntot_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]) == 1 ||
      get_ntot_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) == 1 ) return;

  const MemoryWindow srcWindow( get_window_with_ghost<SrcFieldT >(dim,bcFlag[0],bcFlag[1],bcFlag[2]));
  const MemoryWindow dstWindow( get_window_with_ghost<DestFieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]));

  SrcFieldT  phi      ( srcWindow, NULL );
  DestFieldT fphi     ( dstWindow, NULL );
  DestFieldT fphiExact( dstWindow, NULL );

  switch( FuncType1::FieldType::Location::FaceDir::value ){
  case XDIR::value :
    funcPhi.dx( phi );
    funcFPhi.d2x( fphiExact );
    break;
  case YDIR::value :
    funcPhi.dy( phi );
    funcFPhi.d2y( fphiExact );
    break;
  case ZDIR::value :
    funcPhi.dz( phi );
    funcFPhi.d2z( fphiExact );
    break;
  default:
    assert(0);
  }

  const OpT* const op = opDB.retrieve_operator<OpT>();

  op->apply_to_field( phi, fphi );

  report_errors( fphiExact, fphi );

}

//--------------------------------------------------------------------

void test_ops()
{
  ScratchSVolSSurfX& Sx = *opDB.retrieve_operator<ScratchSVolSSurfX>();
  GradSVolSSurfX&    Gx = *opDB.retrieve_operator<GradSVolSSurfX>();
  Sx.reset_entries(1.0);

  Sx += Gx;
  Sx -= Gx;

  ScratchSVolSSurfY& Sy = *opDB.retrieve_operator<ScratchSVolSSurfY>();
  GradSVolSSurfY&    Gy = *opDB.retrieve_operator<GradSVolSSurfY>();
  Sy += Gy;
}

//--------------------------------------------------------------------

void test_poisson( const Grid& grid,
                   const IntVec& dim,
                   const vector<bool>& bcFlag )
{
# ifdef LINALG_TRILINOS
  //
  // here we use a solution of the form
  //   phi = ax^2 + by^2 + cz^2
  // because this should be solved exactly for a second order scheme.
  //
  // Laplacian(phi) = 2a + 2b + 2c
  //
  const double a=2.0, b=3.0, c=4.0;

  cout << "Setting up Poisson equation test...";

  ScratchSVol& Lx = *opDB.retrieve_operator<ScratchSVol>(1);
  ScratchSVol& Ly = *opDB.retrieve_operator<ScratchSVol>(2);
  ScratchSVol& Lz = *opDB.retrieve_operator<ScratchSVol>(3);

  LinSysInfo lsi( dim, bcFlag[0], bcFlag[1], bcFlag[2] );
  LinearSystem& linsys = LinSysFactory::self().get_linsys( lsi );
  RHS& rhs = linsys.get_rhs();
  LHS& lhs = linsys.get_lhs();
  lhs.reset();

  const GradSVolSSurfX& Gx = *opDB.retrieve_operator<GradSVolSSurfX>();
  const GradSVolSSurfY& Gy = *opDB.retrieve_operator<GradSVolSSurfY>();
  const GradSVolSSurfZ& Gz = *opDB.retrieve_operator<GradSVolSSurfZ>();

  const DivSSurfXSVol& Dx = *opDB.retrieve_operator<DivSSurfXSVol>();  
  const DivSSurfYSVol& Dy = *opDB.retrieve_operator<DivSSurfYSVol>();  
  const DivSSurfZSVol& Dz = *opDB.retrieve_operator<DivSSurfZSVol>();  

  //
  // set up the Laplacian operator in each direction and assemble the
  // linear system to be solved.
  //
  if( dim[0]>1 ){
    Dx.apply_to_op( Gx, Lx );
    lhs.add_op_contribution( Lx );
  }
  if( dim[1]>1 ){
    Dy.apply_to_op( Gy, Ly );
    lhs.add_op_contribution( Ly );
  }
  if( dim[2]>1 ){
    Dz.apply_to_op( Gz, Lz );
    lhs.add_op_contribution( Lz );
  }

  const SVolField& x = grid.xcoord_svol();
  const SVolField& y = grid.ycoord_svol();
  const SVolField& z = grid.zcoord_svol();

  //
  // set the RHS field
  //
  double q=0;
  if( dim[0]>1 ) q+=2*a;
  if( dim[1]>1 ) q+=2*b;
  if( dim[2]>1 ) q+=2*c;
  rhs.reset( q );

  //
  // set the boundary conditions - dirichlet at the cell centers.  In
  // practice we won't know this usually.
  //
  const int ighost = dim[0]>1 ? SVolField::Ghost::NGHOST : 0;
  const int jghost = dim[1]>1 ? SVolField::Ghost::NGHOST : 0;
  const int kghost = dim[2]>1 ? SVolField::Ghost::NGHOST : 0;

  // set bcs: x faces
  if( dim[0]>1 ){
    for( int ix=0; ix<2; ++ix ){
      int i=0;
      if( ix!=0 ) i=dim[0]-1;
      for( int j=0; j<dim[1]; ++j ){
        for( int k=0; k<dim[2]; ++k ){
          const IntVec ijk( i+ighost, j+jghost, k+kghost );
          const int ii = ijk2flat<SVolField>::value(dim,ijk,bcFlag[0],bcFlag[1],bcFlag[2]);
          double bcval = x[ii]*x[ii]*a;
          if( dim[1]>1 ) bcval += y[ii]*y[ii]*b;
          if( dim[2]>1 ) bcval += z[ii]*z[ii]*c;
          const int irow = i + j*dim[0] + k*dim[0]*dim[1];
          linsys.set_dirichlet_condition( irow, bcval );
        }
      }
    }
  }

  // set bcs: y faces
  if( dim[1]>1 ){
    for( int iy=0; iy<2; ++iy ){
      int j=0;
      if( iy!=0 ) j=dim[1]-1;
      for( int i=0; i<dim[0]; ++i ){
        for( int k=0; k<dim[2]; ++k ){
          const IntVec ijk( i+ighost, j+jghost, k+kghost );
          const int ii = ijk2flat<SVolField>::value(dim,ijk,bcFlag[0],bcFlag[1],bcFlag[2]);
          double bcval = y[ii]*y[ii]*b;
          if( dim[0]>1 ) bcval += x[ii]*x[ii]*a;
          if( dim[2]>1 ) bcval += z[ii]*z[ii]*c;
          const int irow = i + j*dim[0] + k*dim[0]*dim[1];
          linsys.set_dirichlet_condition( irow, bcval );
        }
      }
    }
  }

  // set bcs: z faces
  if( dim[2]>1 ){
    for( int iz=0; iz<2; ++iz ){
      int k=0;
      if( iz!=0 ) k=dim[2]-1;
      for( int i=0; i<dim[0]; ++i ){
        for( int j=0; j<dim[1]; ++j ){
          const IntVec ijk( i+ighost, j+jghost, k+kghost );
          const int ii = ijk2flat<SVolField>::value(dim,ijk,bcFlag[0],bcFlag[1],bcFlag[2]);
          double bcval = z[ii]*z[ii]*c;
          if( dim[0]>1 ) bcval += x[ii]*x[ii]*a;
          if( dim[1]>1 ) bcval += y[ii]*y[ii]*b;
          const int irow = i + j*dim[0] + k*dim[0]*dim[1];
          linsys.set_dirichlet_condition( irow, bcval );
        }
      }
    }
  }

//   lhs.Print(cout);


  //
  // Solve the linear system for the solution.
  //
  linsys.solve();

  //
  // examine the solution to determine error
  //
  const SOLN& soln = linsys.get_soln_field();

  SVolField::const_interior_iterator ix = x.interior_begin();
  SVolField::const_interior_iterator iy = y.interior_begin();
  SVolField::const_interior_iterator iz = z.interior_begin();

  SVolField::const_interior_iterator ixe = x.interior_end();

  SOLN::const_interior_iterator isoln = soln.interior_begin();
  RHS::const_interior_iterator   irhs = rhs.interior_begin();

  double maxAbsErr=0.0, maxRelErr=0.0;
  double avgAbsErr=0.0, avgRelErr=0.0;
  int nrel=0, nabs=0;
  for( ; ix!=ixe; ++ix, ++iy, ++iz, ++isoln, ++irhs ){
    const double x=*ix, y=*iy, z=*iz;
    double phi =0;
    if( dim[0]>1 ) phi += a*x*x;
    if( dim[1]>1 ) phi += b*y*y;
    if( dim[2]>1 ) phi += c*z*z;
    const double err = abs(phi-*isoln);
    avgAbsErr += err;
    maxAbsErr = max(err,maxAbsErr);
    ++nabs;
    if( abs(phi)>1e-10 ){
      const double relErr = abs(err / phi);
      avgRelErr += relErr;
      maxRelErr = max(relErr,maxRelErr);
      ++nrel;
    }
  }
  avgRelErr /= double(nrel);
  avgAbsErr /= double(nabs);

  if( maxRelErr>1.0e-10 || maxAbsErr>1.0e-10 ){
    cout << "FAIL" << endl
         << "  max abs error: " << setw(12) << maxAbsErr << "  max rel err: " << setw(12) << maxRelErr << endl
         << "  avg abs error: " << setw(12) << avgAbsErr << "  avg rel err: " << setw(12) << avgRelErr << endl << endl;
  }
  else{
    cout << "PASS." << endl;
  }
# endif // LINALG_TRILINOS
}

//--------------------------------------------------------------------

int main()
{
  IntVec dim(1,1,1);

  cout << "interior nx = "; cin >> dim[0];
  cout << "interior ny = "; cin >> dim[1];
  cout << "interior nz = "; cin >> dim[2];
  cout << endl;

  std::vector<double> length(3,1);
  std::vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    if( dim[i]>1 ) spacing[i] = length[i]/dim[i];
  }

  std::vector<bool> bcFlag(3,true);

  build_ops( dim, spacing, bcFlag, opDB );
  const Grid grid( dim, spacing, bcFlag, opDB );
  //grid.write();

  const MemoryWindow svolWindow( get_window_with_ghost<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) );

  // Scalar-Volume to scalar face gradients and laplacians
  {
    // sin function
    const SinFun<SVolField  > fun     ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );
    const SinFun<SSurfXField> gradFunX( grid.xcoord_sxsurf(), grid.ycoord_sxsurf(), grid.zcoord_sxsurf() );
    const SinFun<SSurfYField> gradFunY( grid.xcoord_sysurf(), grid.ycoord_sysurf(), grid.zcoord_sysurf() );
    const SinFun<SSurfZField> gradFunZ( grid.xcoord_szsurf(), grid.ycoord_szsurf(), grid.zcoord_szsurf() );
    const SinFun<SVolField  > divFun  ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );
    const SinFun<SSurfXField> interpX ( grid.xcoord_sxsurf(), grid.ycoord_sxsurf(), grid.zcoord_sxsurf()  );
    const SinFun<SSurfYField> interpY ( grid.xcoord_sysurf(), grid.ycoord_sysurf(), grid.zcoord_sysurf()  );
    const SinFun<SSurfZField> interpZ ( grid.xcoord_szsurf(), grid.ycoord_szsurf(), grid.zcoord_szsurf()  );

    cout << "=====================================================" << endl
         << "Interpolant scalar volume -> scalar surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_interp_op<InterpSVolSSurfX>( grid, fun, interpX, bcFlag );
    if( dim[1]>1 ) test_interp_op<InterpSVolSSurfY>( grid, fun, interpY, bcFlag );
    if( dim[2]>1 ) test_interp_op<InterpSVolSSurfZ>( grid, fun, interpZ, bcFlag );
    cout << "=====================================================" << endl << endl;


    cout << "=====================================================" << endl
         << "Gradient scalar volume -> scalar surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_grad_op<GradSVolSSurfX>( grid, fun, gradFunX, bcFlag, XDIR::value );
    if( dim[1]>1 ) test_grad_op<GradSVolSSurfY>( grid, fun, gradFunY, bcFlag, YDIR::value );
    if( dim[2]>1 ) test_grad_op<GradSVolSSurfZ>( grid, fun, gradFunZ, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;


    cout << "=====================================================" << endl
         << "Divergence scalar surfaces -> scalar volume" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_div_op<DivSSurfXSVol>( grid, gradFunX, divFun, bcFlag );
    if( dim[1]>1 ) test_div_op<DivSSurfYSVol>( grid, gradFunY, divFun, bcFlag );
    if( dim[2]>1 ) test_div_op<DivSSurfZSVol>( grid, gradFunZ, divFun, bcFlag );
    cout << "=====================================================" << endl << endl;
  }

  {
    const SinFun<SVolField> svolfun( grid.xcoord_svol(), grid.ycoord_svol(), grid.zcoord_svol() );
    const SinFun<XVolField> xvolfun( grid.xcoord_xvol(), grid.ycoord_xvol(), grid.zcoord_xvol() );
    const SinFun<YVolField> yvolfun( grid.xcoord_yvol(), grid.ycoord_yvol(), grid.zcoord_yvol() );
    const SinFun<ZVolField> zvolfun( grid.xcoord_zvol(), grid.ycoord_zvol(), grid.zcoord_zvol() );

    cout << "=====================================================" << endl
         << "Interpolate scalar volume to staggered volumes" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_interp_op<InterpSVolXVol>( grid, svolfun, xvolfun, bcFlag );
    if( dim[1]>1 ) test_interp_op<InterpSVolYVol>( grid, svolfun, yvolfun, bcFlag );
    if( dim[2]>1 ) test_interp_op<InterpSVolZVol>( grid, svolfun, zvolfun, bcFlag );

    cout << endl
         << "Gradient scalar volume to staggered volumes" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_grad_op<GradSVolXVol>( grid, svolfun, xvolfun, bcFlag, XDIR::value );
    if( dim[1]>1 ) test_grad_op<GradSVolYVol>( grid, svolfun, yvolfun, bcFlag, YDIR::value );
    if( dim[2]>1 ) test_grad_op<GradSVolZVol>( grid, svolfun, zvolfun, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;
  }

  {
    const SinFun<SVolField>  svolfun( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );

    const SinFun<XSurfXField> xsurfx( grid.xcoord_xxsurf(), grid.ycoord_xxsurf(), grid.zcoord_xxsurf() );
    const SinFun<XSurfYField> xsurfy( grid.xcoord_xysurf(), grid.ycoord_xysurf(), grid.zcoord_xysurf() );
    const SinFun<XSurfZField> xsurfz( grid.xcoord_xzsurf(), grid.ycoord_xzsurf(), grid.zcoord_xzsurf() );

    const SinFun<YSurfXField> ysurfx( grid.xcoord_yxsurf(), grid.ycoord_yxsurf(), grid.zcoord_yxsurf() );
    const SinFun<YSurfYField> ysurfy( grid.xcoord_yysurf(), grid.ycoord_yysurf(), grid.zcoord_yysurf() );
    const SinFun<YSurfZField> ysurfz( grid.xcoord_yzsurf(), grid.ycoord_yzsurf(), grid.zcoord_yzsurf() );

    const SinFun<ZSurfXField> zsurfx( grid.xcoord_zxsurf(), grid.ycoord_zxsurf(), grid.zcoord_zxsurf() );
    const SinFun<ZSurfYField> zsurfy( grid.xcoord_zysurf(), grid.ycoord_zysurf(), grid.zcoord_zysurf() );
    const SinFun<ZSurfZField> zsurfz( grid.xcoord_zzsurf(), grid.ycoord_zzsurf(), grid.zcoord_zzsurf() );

    cout << "=====================================================" << endl
         << "Interpolate scalar volume to staggered surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ){
      if( dim[0]>1 ) test_interp_op<InterpSVolXSurfX>( grid, svolfun, xsurfx, bcFlag );
      if( dim[1]>1 ) test_interp_op<InterpSVolXSurfY>( grid, svolfun, xsurfy, bcFlag );
      if( dim[2]>1 ) test_interp_op<InterpSVolXSurfZ>( grid, svolfun, xsurfz, bcFlag );
    }

    if( dim[2]>1 ){
      if( dim[0]>1 ) test_interp_op<InterpSVolYSurfX>( grid, svolfun, ysurfx, bcFlag );
      if( dim[1]>1 ) test_interp_op<InterpSVolYSurfY>( grid, svolfun, ysurfy, bcFlag );
      if( dim[2]>1 ) test_interp_op<InterpSVolYSurfZ>( grid, svolfun, ysurfz, bcFlag );
    }
    if( dim[2]>1 ){
      if( dim[0]>1 ) test_interp_op<InterpSVolZSurfX>( grid, svolfun, zsurfx, bcFlag );
      if( dim[1]>1 ) test_interp_op<InterpSVolZSurfY>( grid, svolfun, zsurfy, bcFlag );
      if( dim[2]>1 ) test_interp_op<InterpSVolZSurfZ>( grid, svolfun, zsurfz, bcFlag );
    }
    cout << "=====================================================" << endl << endl;
  }

  {
    const SinFun<SSurfXField> fsx( grid.xcoord_sxsurf(), grid.ycoord_sxsurf(), grid.zcoord_sxsurf() );
    const SinFun<SSurfYField> fsy( grid.xcoord_sysurf(), grid.ycoord_sysurf(), grid.zcoord_sysurf() );
    const SinFun<SSurfZField> fsz( grid.xcoord_szsurf(), grid.ycoord_szsurf(), grid.zcoord_szsurf() );
    const SinFun<XVolField  > fxv( grid.xcoord_xvol(),   grid.ycoord_xvol(),   grid.zcoord_xvol()   );
    const SinFun<YVolField  > fyv( grid.xcoord_yvol(),   grid.ycoord_yvol(),   grid.zcoord_yvol()   );
    const SinFun<ZVolField  > fzv( grid.xcoord_zvol(),   grid.ycoord_zvol(),   grid.zcoord_zvol()   );
    cout << "=====================================================" << endl
         << "Interpolate scalar surface to staggered volumes" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_interp_op<InterpSSurfXXVol>( grid, fsx, fxv, bcFlag );
    if( dim[1]>1 ) test_interp_op<InterpSSurfYYVol>( grid, fsy, fyv, bcFlag );
    if( dim[2]>1 ) test_interp_op<InterpSSurfZZVol>( grid, fsz, fzv, bcFlag );
    cout << "=====================================================" << endl << endl;
  }

  {
    const SinFun<SVolField> svolf( grid.xcoord_svol(), grid.ycoord_svol(), grid.zcoord_svol() );
    const SinFun<XVolField> xvolf( grid.xcoord_xvol(), grid.ycoord_xvol(), grid.zcoord_xvol() );
    const SinFun<YVolField> yvolf( grid.xcoord_yvol(), grid.ycoord_yvol(), grid.zcoord_yvol() );
    const SinFun<ZVolField> zvolf( grid.xcoord_zvol(), grid.ycoord_zvol(), grid.zcoord_zvol() );
    cout << "=====================================================" << endl
         << "Gradient staggered volumes to scalar volume (dilatation)" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_grad_op<GradXVolSVol>( grid, xvolf, svolf, bcFlag, XDIR::value );
    if( dim[1]>1 ) test_grad_op<GradYVolSVol>( grid, yvolf, svolf, bcFlag, YDIR::value );
    if( dim[2]>1 ) test_grad_op<GradZVolSVol>( grid, zvolf, svolf, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl;
  }

  // x-Volume to y-volume x-surface and z-volume x-surface
  if( dim[0]>1 && dim[1]>1 || dim[2]>1 ){
    const SinFun<XVolField  >   xvolfun( grid.xcoord_xvol(),   grid.ycoord_xvol(),   grid.zcoord_xvol()   );
    const SinFun<YSurfXField> ysurfxfun( grid.xcoord_yxsurf(), grid.ycoord_yxsurf(), grid.zcoord_yxsurf() );
    const SinFun<ZSurfXField> zsurfxfun( grid.xcoord_zxsurf(), grid.ycoord_zxsurf(), grid.zcoord_zxsurf() );
    cout << "=====================================================" << endl
         << "Interpolate x volume to y and z volume x-surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpXVolYSurfX>( grid, xvolfun, ysurfxfun, bcFlag );
    test_interp_op<InterpXVolZSurfX>( grid, xvolfun, zsurfxfun, bcFlag );

    cout << endl
         << "Gradient xvol->yvol x surf and xvol->zvol x surf" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradXVolYSurfX>( grid, xvolfun, ysurfxfun, bcFlag, YDIR::value );
    test_grad_op<GradXVolZSurfX>( grid, xvolfun, zsurfxfun, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;
  }

  // y-volume to x-volume y-surface z-volume y-surface
  if( dim[1]>1 && dim[0]>1 || dim[2]>1 ){
    const SinFun<YVolField  >   yvolfun( grid.xcoord_yvol(),   grid.ycoord_yvol(),   grid.zcoord_yvol()   );
    const SinFun<XSurfYField> xsurfyfun( grid.xcoord_xysurf(), grid.ycoord_xysurf(), grid.zcoord_xysurf() );
    const SinFun<ZSurfYField> zsurfyfun( grid.xcoord_zysurf(), grid.ycoord_zysurf(), grid.zcoord_zysurf() );
    cout << "=====================================================" << endl
         << "Interpolate y volume to x and z volume y-surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpYVolXSurfY>( grid, yvolfun, xsurfyfun, bcFlag );
    test_interp_op<InterpYVolZSurfY>( grid, yvolfun, zsurfyfun, bcFlag );

    cout << endl
         << "Gradient yvol->xvol y surf and yvol->zvol y surf" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradYVolXSurfY>( grid, yvolfun, xsurfyfun, bcFlag, XDIR::value );
    test_grad_op<GradYVolZSurfY>( grid, yvolfun, zsurfyfun, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;
  }

  // z-volume to x-volume z-surface y-volume z-surface
  if( dim[2]>1 && dim[0]>1 || dim[1]>1 ){
    const SinFun<ZVolField  >   zvolfun( grid.xcoord_zvol(),   grid.ycoord_zvol(),   grid.zcoord_zvol()   );
    const SinFun<XSurfZField> xsurfzfun( grid.xcoord_xzsurf(), grid.ycoord_xzsurf(), grid.zcoord_xzsurf() );
    const SinFun<YSurfZField> ysurfzfun( grid.xcoord_yzsurf(), grid.ycoord_yzsurf(), grid.zcoord_yzsurf() );
    cout << "=====================================================" << endl
         << "Interpolate z volume to x and y volume z-surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    test_interp_op<InterpZVolXSurfZ>( grid, zvolfun, xsurfzfun, bcFlag );
    test_interp_op<InterpZVolYSurfZ>( grid, zvolfun, ysurfzfun, bcFlag );

    cout << endl
         << "Gradient zvol->xvol z surf and zvol->yvol z surf" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    test_grad_op<GradZVolXSurfZ>( grid, zvolfun, xsurfzfun, bcFlag, XDIR::value );
    test_grad_op<GradZVolYSurfZ>( grid, zvolfun, ysurfzfun, bcFlag, YDIR::value );
    cout << "=====================================================" << endl << endl;
  }

  // X-volume to x-surface face component gradients
  if( dim[0]>1 ){

    // sin function
    const SinFun<XVolField  > sinfun ( grid.xcoord_xvol(),   grid.ycoord_xvol(),   grid.zcoord_xvol()   );
    const SinFun<XSurfXField> gradX  ( grid.xcoord_xxsurf(), grid.ycoord_xxsurf(), grid.zcoord_xxsurf() );
    const SinFun<XSurfYField> gradY  ( grid.xcoord_xysurf(), grid.ycoord_xysurf(), grid.zcoord_xysurf() );
    const SinFun<XSurfZField> gradZ  ( grid.xcoord_xzsurf(), grid.ycoord_xzsurf(), grid.zcoord_xzsurf() );
    const SinFun<XVolField  > divFun ( grid.xcoord_xvol(),   grid.ycoord_xvol(),   grid.zcoord_xvol()   );
    const SinFun<XSurfXField> interpX( grid.xcoord_xxsurf(), grid.ycoord_xxsurf(), grid.zcoord_xxsurf() );
    const SinFun<XSurfYField> interpY( grid.xcoord_xysurf(), grid.ycoord_xysurf(), grid.zcoord_xysurf() );
    const SinFun<XSurfZField> interpZ( grid.xcoord_xzsurf(), grid.ycoord_xzsurf(), grid.zcoord_xzsurf() );
//     const SinFun<XSurfField > interp ( grid.xcoord_xsurf(),  grid.ycoord_xsurf(),  grid.zcoord_xsurf()  );

    cout << "=====================================================" << endl
         << "Gradient x-volume -> x-surface" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_grad_op<GradXVolXSurfX>( grid, sinfun, gradX, bcFlag, XDIR::value );
    if( dim[1]>1 ) test_grad_op<GradXVolXSurfY>( grid, sinfun, gradY, bcFlag, YDIR::value );
    if( dim[2]>1 ) test_grad_op<GradXVolXSurfZ>( grid, sinfun, gradZ, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
         << "Divergence x-surface -> x-volume" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_div_op<DivXSurfXXVol>( grid, gradX, divFun, bcFlag );
    if( dim[1]>1 ) test_div_op<DivXSurfYXVol>( grid, gradY, divFun, bcFlag );
    if( dim[2]>1 ) test_div_op<DivXSurfZXVol>( grid, gradZ, divFun, bcFlag );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
         << "Interpolate x-volume -> x-surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_interp_op<InterpXVolXSurfX>( grid, sinfun, interpX, bcFlag );
    if( dim[1]>1 ) test_interp_op<InterpXVolXSurfY>( grid, sinfun, interpY, bcFlag );
    if( dim[2]>1 ) test_interp_op<InterpXVolXSurfZ>( grid, sinfun, interpZ, bcFlag );
    cout << "=====================================================" << endl << endl;
  }

  // Y-volume to y-surface face component gradients
  if( dim[1]>1 ){

    // sin function
    const SinFun<YVolField  > sinfun ( grid.xcoord_yvol(),   grid.ycoord_yvol(),   grid.zcoord_yvol()   );
    const SinFun<YSurfXField> gradX  ( grid.xcoord_yxsurf(), grid.ycoord_yxsurf(), grid.zcoord_yxsurf() );
    const SinFun<YSurfYField> gradY  ( grid.xcoord_yysurf(), grid.ycoord_yysurf(), grid.zcoord_yysurf() );
    const SinFun<YSurfZField> gradZ  ( grid.xcoord_yzsurf(), grid.ycoord_yzsurf(), grid.zcoord_yzsurf() );
    const SinFun<YVolField  > divFun ( grid.xcoord_yvol(),   grid.ycoord_yvol(),   grid.zcoord_yvol()   );
    const SinFun<YSurfXField> interpX( grid.xcoord_yxsurf(), grid.ycoord_yxsurf(), grid.zcoord_yxsurf() );
    const SinFun<YSurfYField> interpY( grid.xcoord_yysurf(), grid.ycoord_yysurf(), grid.zcoord_yysurf() );
    const SinFun<YSurfZField> interpZ( grid.xcoord_yzsurf(), grid.ycoord_yzsurf(), grid.zcoord_yzsurf() );

    cout << "=====================================================" << endl
         << "Gradient y-volume -> y-surface" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_grad_op<GradYVolYSurfX>( grid, sinfun, gradX, bcFlag, XDIR::value );
    if( dim[1]>1 ) test_grad_op<GradYVolYSurfY>( grid, sinfun, gradY, bcFlag, YDIR::value );
    if( dim[2]>1 ) test_grad_op<GradYVolYSurfZ>( grid, sinfun, gradZ, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
         << "Divergence y-surface -> y-volume" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_div_op<DivYSurfXYVol>( grid, gradX, divFun, bcFlag );
    if( dim[1]>1 ) test_div_op<DivYSurfYYVol>( grid, gradY, divFun, bcFlag );
    if( dim[2]>1 ) test_div_op<DivYSurfZYVol>( grid, gradZ, divFun, bcFlag );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
         << "Interpolate y-volume -> y-surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_interp_op<InterpYVolYSurfX>( grid, sinfun, interpX, bcFlag );
    if( dim[1]>1 ) test_interp_op<InterpYVolYSurfY>( grid, sinfun, interpY, bcFlag );
    if( dim[2]>1 ) test_interp_op<InterpYVolYSurfZ>( grid, sinfun, interpZ, bcFlag );
    cout << "=====================================================" << endl << endl;
  }

  // Z-volume to z-surface face component gradients
  if( dim[2]>1 ){

    // sin function
    const SinFun<ZVolField  > sinfun ( grid.xcoord_zvol(),   grid.ycoord_zvol(),   grid.zcoord_zvol()   );
    const SinFun<ZSurfXField> gradX  ( grid.xcoord_zxsurf(), grid.ycoord_zxsurf(), grid.zcoord_zxsurf() );
    const SinFun<ZSurfYField> gradY  ( grid.xcoord_zysurf(), grid.ycoord_zysurf(), grid.zcoord_zysurf() );
    const SinFun<ZSurfZField> gradZ  ( grid.xcoord_zzsurf(), grid.ycoord_zzsurf(), grid.zcoord_zzsurf() );
    const SinFun<ZVolField  > divFun ( grid.xcoord_zvol(),   grid.ycoord_zvol(),   grid.zcoord_zvol()   );
    const SinFun<ZSurfXField> interpX( grid.xcoord_zxsurf(), grid.ycoord_zxsurf(), grid.zcoord_zxsurf() );
    const SinFun<ZSurfYField> interpY( grid.xcoord_zysurf(), grid.ycoord_zysurf(), grid.zcoord_zysurf() );
    const SinFun<ZSurfZField> interpZ( grid.xcoord_zzsurf(), grid.ycoord_zzsurf(), grid.zcoord_zzsurf() );

    cout << "=====================================================" << endl
         << "Gradient z-volume -> z-surface" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_grad_op<GradZVolZSurfX>( grid, sinfun, gradX, bcFlag, XDIR::value );
    if( dim[1]>1 ) test_grad_op<GradZVolZSurfY>( grid, sinfun, gradY, bcFlag, YDIR::value );
    if( dim[2]>1 ) test_grad_op<GradZVolZSurfZ>( grid, sinfun, gradZ, bcFlag, ZDIR::value );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
         << "Divergence z-surface -> z-volume" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_div_op<DivZSurfXZVol>( grid, gradX, divFun, bcFlag );
    if( dim[1]>1 ) test_div_op<DivZSurfYZVol>( grid, gradY, divFun, bcFlag );
    if( dim[2]>1 ) test_div_op<DivZSurfZZVol>( grid, gradZ, divFun, bcFlag );
    cout << "=====================================================" << endl << endl;

    cout << "=====================================================" << endl
         << "Interpolate z-volume -> z-surfaces" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    if( dim[0]>1 ) test_interp_op<InterpZVolZSurfX>( grid, sinfun, interpX, bcFlag );
    if( dim[1]>1 ) test_interp_op<InterpZVolZSurfY>( grid, sinfun, interpY, bcFlag );
    if( dim[2]>1 ) test_interp_op<InterpZVolZSurfZ>( grid, sinfun, interpZ, bcFlag );
    cout << "=====================================================" << endl << endl;
  }

  {
    cout << endl << "Scalar volume scratch operators" << endl
         << " max abs err | max rel err | avg abs err | avg rel err |" << endl
         << "-------------|-------------|-------------|-------------|" << endl;
    ScratchSVol* const Ssx = opDB.retrieve_operator<ScratchSVol>(1);
    ScratchSVol* const Ssy = opDB.retrieve_operator<ScratchSVol>(2);
    ScratchSVol* const Ssz = opDB.retrieve_operator<ScratchSVol>(3);

    const GradSVolSSurfX* const Gsx = opDB.retrieve_operator<GradSVolSSurfX>();
    const GradSVolSSurfY* const Gsy = opDB.retrieve_operator<GradSVolSSurfY>();
    const GradSVolSSurfZ* const Gsz = opDB.retrieve_operator<GradSVolSSurfZ>();

    const DivSSurfXSVol * const Dsx = opDB.retrieve_operator<DivSSurfXSVol>();
    const DivSSurfYSVol * const Dsy = opDB.retrieve_operator<DivSSurfYSVol>();
    const DivSSurfZSVol * const Dsz = opDB.retrieve_operator<DivSSurfZSVol>();

    if( dim[0]>1 ) Dsx->apply_to_op( *Gsx, *Ssx );
    if( dim[1]>1 ) Dsy->apply_to_op( *Gsy, *Ssy );
    if( dim[2]>1 ) Dsz->apply_to_op( *Gsz, *Ssz );

//     Dsx->write_matlab("Dx");
//     Gsx->write_matlab("Gx");

    const SinFun<SVolField  > fun     ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );
    const SinFun<SVolField  > divFun  ( grid.xcoord_svol(),   grid.ycoord_svol(),   grid.zcoord_svol()   );

    SVolField phi       ( svolWindow, NULL );
    SVolField d2phi     ( svolWindow, NULL );
    SVolField d2phiExact( svolWindow, NULL );

//     write_matlab(phi,"phi");

    fun.evaluate( phi );

    if( dim[0]>1 ){
      Ssx->apply_to_field( phi, d2phi );
      divFun.d2x( d2phiExact );
      report_errors( d2phiExact, d2phi );
    }
    if( dim[1]>1 ){
      Ssy->apply_to_field( phi, d2phi );
      divFun.d2y( d2phiExact );
      report_errors( d2phiExact, d2phi );
    }
    if( dim[2]>1 ){
      Ssz->apply_to_field( phi, d2phi );
      divFun.d2z( d2phiExact );
      report_errors( d2phiExact, d2phi );
    }
  }

  {
    cout << endl << "x-volume scratch operators ... ";
    ScratchXVol* const Sxx = opDB.retrieve_operator<ScratchXVol>( 1 );
    ScratchXVol* const Sxy = opDB.retrieve_operator<ScratchXVol>( 2 );
    ScratchXVol* const Sxz = opDB.retrieve_operator<ScratchXVol>( 3 );

    const GradXVolXSurfX* const Gsx = opDB.retrieve_operator<GradXVolXSurfX>();
    const GradXVolXSurfY* const Gsy = opDB.retrieve_operator<GradXVolXSurfY>();
    const GradXVolXSurfZ* const Gsz = opDB.retrieve_operator<GradXVolXSurfZ>();

    const DivXSurfXXVol * const Dsx = opDB.retrieve_operator<DivXSurfXXVol>();
    const DivXSurfYXVol * const Dsy = opDB.retrieve_operator<DivXSurfYXVol>();
    const DivXSurfZXVol * const Dsz = opDB.retrieve_operator<DivXSurfZXVol>();

    if( dim[0]>1 ) Dsx->apply_to_op( *Gsx, *Sxx );
    if( dim[1]>1 ) Dsy->apply_to_op( *Gsy, *Sxy );
    if( dim[2]>1 ) Dsz->apply_to_op( *Gsz, *Sxz );

    cout << "done" << endl;
  }


  test_ops();
  test_poisson( grid, dim, bcFlag );
}
