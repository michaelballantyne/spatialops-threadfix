#include <spatialops/structured/FVOneDimensional.h>

#include "Functions.h"

#include <iostream>
#include <iomanip>
#include <vector>

//====================================================================

using namespace std;
using namespace SpatialOps;
using namespace structured;

//====================================================================

template<typename FieldT>
void compare( const FieldT& f1, const FieldT& f2,
              double& maxAbsErr, double& maxRelErr,
              double& avgAbsErr, double& avgRelErr )
{
  static const double SMALL = numeric_limits<double>::epsilon() * 10.0;

  maxAbsErr = 0.0;  maxRelErr = 0.0;
  avgAbsErr = 0.0;  avgRelErr = 0.0;

  int npts = 0;

  typename FieldT::const_interior_iterator if1 = f1.interior_begin();
  typename FieldT::const_interior_iterator if2 = f2.interior_begin();
  const typename FieldT::const_interior_iterator if1e= f1.interior_end();

  for( ; if1!=if1e; ++if1, ++if2 ){
    const double absErr = abs( *if1 - *if2 );
    const double relErr = absErr / (abs(*if1)+SMALL);
    maxAbsErr = max( maxAbsErr, absErr );
    maxRelErr = max( maxRelErr, relErr );
    avgRelErr += relErr;
    avgAbsErr += absErr;
    ++npts;
  }
  avgRelErr /= npts;
  avgAbsErr /= npts;
}

//====================================================================

void
setup_mesh( const int npts,
            vector<double>& xv,
            vector<double>& xs,
            const double scaleFac=1.0 )
{
  xv.resize(npts,0.0);
  xs.resize(npts+1,0.0);

  const double xend = 3.1415;

  if( scaleFac==1.0 ){
    const double dx = xend / npts;
    vector<double>::iterator ixv=xv.begin(), ixs=xs.begin();
    for( int i=0; i<npts; ++i, ++ixs, ++ixv ){
      *ixs = i*dx;
      *ixv = *ixs + 0.5*dx;
      //     cout << *ixs << "  " << *ixv << endl;
    }
    *ixs = npts*dx;
    //   cout << *ixs << endl;
  }
  else{
    // stretched mesh
    // first create a uniform mesh
    const double dxtmp = 1.0/npts;
    vector<double>::iterator ixs=xs.begin();
    for( int i=0; i<=npts; ++i, ++ixs ){
      *ixs = i*dxtmp;
    }
    const double xcrit = 0.23;
    const double xo = 1.0/(2.0*scaleFac) * log( (1.+(exp( scaleFac)-1.)*xcrit) /
                                                (1.+(exp(-scaleFac)-1.)*xcrit) );
    const double A = sinh(scaleFac*xo);
    xs[0]=0.0;
    for( ixs=xs.begin()+1; ixs!=xs.end()-1; ++ixs ){
      *ixs = xcrit/A * (sinh(scaleFac*(*ixs-xo)) + A);
      *ixs *= xend;
      assert( *ixs >= 0.0  && *ixs <= xend );
    }
    *ixs = xend;
    ixs=xs.begin();
    for( vector<double>::iterator ixv=xv.begin(); ixv!=xv.end(); ++ixv, ++ixs ){
      *ixv = 0.5*(*ixs+*(ixs+1));
    }
  }
}

//====================================================================

void calculate_fields( const int npts,
                       const int order,
                       const double scaleFac,
                       double& interpErr,
                       double& gradErr,
                       double& divErr,
                       double& fcinterpErr,
                       double& fcgradErr,
                       double& ccgradErr )
{
  vector<double> xv, xs;
  setup_mesh( npts, xv, xs, scaleFac );

  //
  // build operators
  //
  typedef InterpSVolSSurfX            Interp;
  typedef GradSVolSSurfX              Grad;
  typedef DivSSurfXSVol               Div;
  typedef OneDimInterpolantAssembler  InterpAssembler;
  typedef OneDimGradientAssembler     GradAssembler;
  typedef OneDimDivergenceAssembler   DivAssembler;

  const int nGhost = SVolField::Ghost::NGHOST;

  InterpAssembler Rxa(order,nGhost,xv,xs);
  GradAssembler   Gxa(order,nGhost,xv,xs);
  DivAssembler    Dxa(nGhost,xv,xs);

  const Interp Rx(Rxa);
  const Grad   Gx(Gxa);
  const Div    Dx(Dxa);

  //
  // build fields
  // 
  IntVec dim(npts,1,1);

  const MemoryWindow svolWindow( get_window_with_ghost<SVolField  >(dim,true,true,true) );
  const MemoryWindow ssxWindow ( get_window_with_ghost<SSurfXField>(dim,true,true,true) );

  SVolField          f( svolWindow, NULL );
  SSurfXField  fxexact( ssxWindow,  NULL );
  SSurfXField fxinterp( ssxWindow,  NULL );
  SSurfXField    fgrad( ssxWindow,  NULL );
  SSurfXField  fgexact( ssxWindow,  NULL );
  SVolField   d2fexact( svolWindow, NULL );
  SVolField        d2f( svolWindow, NULL );

  const SVolField    xvol( svolWindow, &get_x_src (nGhost,xv,xs)[0] );
  const SSurfXField xsurf( ssxWindow,  &get_x_dest(nGhost,xv,xs)[0] );
  SVolField    yvol( svolWindow, NULL );  yvol=0.0;
  SVolField    zvol( svolWindow, NULL );  zvol=0.0;
  SSurfXField ysurf( ssxWindow,  NULL );  ysurf=0.0;
  SSurfXField zsurf( ssxWindow,  NULL );  zsurf=0.0;

  //
  // build mesh functions for testing and get exact solutions.
  //
  SinFun<SVolField  >  sinvol( xvol, yvol, zvol );
  SinFun<SSurfXField> sinsurf( xsurf, ysurf, zsurf );
  sinvol.evaluate( f );
  sinsurf.evaluate( fxexact );
  sinsurf.dx( fgexact  );
  sinvol.d2x( d2fexact );

  //
  // Get the field values from the operators
  //
  Rx.apply_to_field( f, fxinterp );
  Gx.apply_to_field( f, fgrad );
  Dx.apply_to_field( fgrad, d2f );

  double avgAbsErr, avgRelErr, maxRelErr;
  compare( fxexact, fxinterp, interpErr, maxRelErr, avgAbsErr, avgRelErr );
  compare( fgexact, fgrad,    gradErr,   maxRelErr, avgAbsErr, avgRelErr );
  compare( d2fexact, d2f,     divErr,    maxRelErr, avgAbsErr, avgRelErr );

  /*
  //
  // Dump results
  //
  Rx.write_matlab("Rx");
  Gx.write_matlab("Gx");
  Dx.write_matlab("Dx");
//   cout << endl
//        << "Interpolant operator has " << Rxa.num_nonzeros() << " nonzero columns" << endl
//        << "Gradient    operator has " << Gxa.num_nonzeros() << " nonzero columns" << endl
//        << "Divergence  operator has " << Dxa.num_nonzeros() << " nonzero columns" << endl
//        << endl;
//   xvol.write_matlab("x");
//   xsurf.write_matlab("xs");
  f.write_matlab("f");
  fxexact.write_matlab("fxexact");
  fxinterp.write_matlab("fxinterp");
  fgexact.write_matlab("dfexact");
  fgrad.write_matlab("df");
  d2fexact.write_matlab("d2fexact");
  d2f.write_matlab("d2f");
  */

  InterpAssembler Rxsva(order,nGhost,xs,xv);
  typedef SpatialOperator< LinAlg, Interpolant, SSurfXField, SVolField > InterpFace2Vol;
  const InterpFace2Vol Rxsv(Rxsva);
  Rxsv.apply_to_field( fxexact, d2f );
  compare( f, d2f, fcinterpErr, maxRelErr, avgAbsErr, avgRelErr );

  /*
  f.write_matlab("ficell");
  d2f.write_matlab("ficellexact");
  */

  GradAssembler Gxsva(order,nGhost,xs,xv);
  typedef SpatialOperator< LinAlg, Gradient, SSurfXField, SVolField > GradFace2Vol;
  const GradFace2Vol Gxsv(Gxsva);
  Gxsv.apply_to_field( fxexact, f );
  sinvol.dx( d2fexact );
  compare( f, d2fexact, fcgradErr, maxRelErr, avgAbsErr, avgRelErr );
  /*
  f.write_matlab( "dfcell" );
  d2fexact.write_matlab( "dfcellexact" );
  xvol.write_matlab( "xv" );
  xsurf.write_matlab("xs");
  Gxsv.write_matlab("Gxsv");
  Rxsv.write_matlab("Rxsv");
  */

  GradAssembler Gxvxva(order,nGhost,xv,xv);
  typedef SpatialOperator<LinAlg,Gradient,SVolField,SVolField> GradC2C;
  const GradC2C Gvv( Gxvxva );
  sinvol.evaluate(f);
  Gvv.apply_to_field( f, d2f );
  compare( d2f, d2fexact, ccgradErr, maxRelErr, avgAbsErr, avgRelErr );
//   Gvv.write_matlab( "Gvv" );
//   d2f.write_matlab("dfdx"); d2fexact.write_matlab("dfdxexact");
//   xvol.write_matlab("xv");
}

//====================================================================

void driver( const double scalefac,
             const int order,
             ostream& os )
{
  const int n=8;
  const int npts [] = { 20, 40, 80, 160, 320, 640, 1280, 2560 };

  os << "# scalefac = " << scalefac << endl
     << "# order = " << order << endl
     << "# npts   Interp     Grad       Div    " << endl;

  double idealErr = 0.1;
  double interpErr, gradErr, divErr, fciErr, fcgErr, ccgErr;
  for( int i=0; i<n; ++i ){
    idealErr /= order;
    calculate_fields( npts[i], order, scalefac, interpErr, gradErr, divErr, fciErr, fcgErr, ccgErr );
    os << setw(6) << npts[i]
       << scientific << setprecision(2)
       << setw(11) << interpErr
       << setw(11) << gradErr
       << setw(11) << divErr
       << setw(11) << idealErr
       << endl;
  }
  os << endl;
}

//====================================================================

void check_err( ofstream& fout, const int npts, const int order, const double scalefac )
{
  double ierr, gerr, derr, fcierr, fcgerr, ccgerr;
  cout << "running test for " << npts << " points and polynomial order " << order
       << " with stretch factor " << scalefac << endl;
  calculate_fields( npts, order, scalefac, ierr, gerr, derr, fcierr, fcgerr, ccgerr );
  fout << "npts="  << npts << ", poly order=" << order << ", stretch factor=" << scalefac << endl;
  fout << "interpolation    : " << ierr << endl
       << "gradient         : " << gerr << endl
       << "divergence       : " << derr << endl
       << "face-cell interp : " << fcierr << endl
       << "face-cell grad   : " << fcgerr << endl
       << "cell-cell grad   : " << ccgerr << endl
       << endl;
}

//====================================================================

int main()
{
  /*
  // jcs could do this once we get a test harness built:
  ofstream os("OneDimOpsConvergence.out",ios::trunc);
  driver( 1.0, 2, os );
  driver( 4.0, 2, os );
  driver( 1.0, 4, os );
  driver( 4.2, 4, os );

  const int val = system( "diff -q OneDimOpsConvergence.out OneDimOpsConvergence.gold" );
  cout << "convergence test on 1D operators on uniform and nonuniform meshes ... ";
  if( val==0 )
    cout << "PASS" << endl;
  else
    cout << "FAIL" << endl;
  */

  ofstream fout( "results_test1d.out" );
  check_err( fout, 20, 2, 1.0 );
  check_err( fout, 40, 2, 1.0 );
  check_err( fout, 80, 2, 1.0 );

  check_err( fout, 20, 2, 4.0 );
  check_err( fout, 40, 2, 4.0 );
  check_err( fout, 80, 2, 4.0 );

  check_err( fout, 20, 4, 1.0 );
  check_err( fout, 40, 4, 1.0 );
  check_err( fout, 80, 4, 1.0 );

  check_err( fout, 20, 4, 4.0 );
  check_err( fout, 40, 4, 4.0 );
  check_err( fout, 80, 4, 4.0 );

  const int val = system( "diff -q results_test1d.out results_test1d.gold" );
  cout << "convergence test on 1D operators on uniform and nonuniform meshes ... ";

  int returncode = 0;
  if( val==0 ){
    cout << "PASS" << endl;
    returncode = 0;
  }
  else{
    cout << "FAIL" << endl;
    returncode = -1;
  }
  return returncode;
}

//====================================================================
