#include <FVOneDimensional.h>

#include <Functions.h>

#include <iostream>
#include <iomanip>
#include <vector>

//====================================================================

using namespace std;
using namespace SpatialOps;
using namespace FVStaggered;

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
  //const double xend = 1.0;

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
                       double& divErr )
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

  InterpAssembler Rxa(order,xv,xs);
  GradAssembler   Gxa(order,xv,xs);
  DivAssembler    Dxa(xv,xs);

  const Interp Rx(Rxa);
  const Grad   Gx(Gxa);
  const Div    Dx(Dxa);

  //
  // build fields
  // 
  vector<int> dim(3,1); dim[0]=npts;
  const int ntotSVol  = get_n_tot<SVolField  >(dim,true,true,true);
  const int ntotSSurf = get_n_tot<SSurfXField>(dim,true,true,true);
  const set<int> ghostSVol  = get_ghost_set<SVolField  >(dim,true,true,true);
  const set<int> ghostSSurf = get_ghost_set<SSurfXField>(dim,true,true,true);

  SVolField          f( ntotSVol,  ghostSVol,  NULL );
  SSurfXField  fxexact( ntotSSurf, ghostSSurf, NULL );
  SSurfXField fxinterp( ntotSSurf, ghostSSurf, NULL );
  SSurfXField    fgrad( ntotSSurf, ghostSSurf, NULL );
  SSurfXField  fgexact( ntotSSurf, ghostSSurf, NULL );
  SVolField   d2fexact( ntotSVol,  ghostSVol,  NULL );
  SVolField        d2f( ntotSVol,  ghostSVol,  NULL );

  const SVolField    xvol( ntotSVol,  ghostSVol,  &Rxa.get_x (xv,xs)[0] );
  const SSurfXField xsurf( ntotSSurf, ghostSSurf, &Rxa.get_xs(xv,xs)[0] );
  SVolField    yvol( ntotSVol,  ghostSVol,  NULL );  yvol=0.0;
  SVolField    zvol( ntotSVol,  ghostSVol,  NULL );  zvol=0.0;
  SSurfXField ysurf( ntotSSurf, ghostSSurf, NULL );  ysurf=0.0;
  SSurfXField zsurf( ntotSSurf, ghostSSurf, NULL );  zsurf=0.0;

  //
  // build mesh functions for testing and get exact solutions.
  //
  SinFun<SVolField  >  sinvol( xvol, yvol, zvol );
  SinFun<SSurfXField> sinsurf( xsurf, ysurf, zsurf );
  sinvol.evaluate( f );
  sinsurf.evaluate( fxexact );
  sinsurf.dx( fgexact );
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
  cout << endl
       << "Interpolant operator has " << Rxa.num_nonzeros() << " nonzero columns" << endl
       << "Gradient    operator has " << Gxa.num_nonzeros() << " nonzero columns" << endl
       << "Divergence  operator has " << Dxa.num_nonzeros() << " nonzero columns" << endl
       << endl;
  xvol.write_matlab("x");
  xsurf.write_matlab("xs");
  f.write_matlab("f");
  fxexact.write_matlab("fxexact");
  fxinterp.write_matlab("fxinterp");
  fgexact.write_matlab("dfexact");
  fgrad.write_matlab("df");
  d2fexact.write_matlab("d2fexact");
  d2f.write_matlab("d2f");
  */
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
  double interpErr, gradErr, divErr;
  for( int i=0; i<n; ++i ){
    idealErr /= order;
    calculate_fields( npts[i], order, scalefac, interpErr, gradErr, divErr );
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

bool check_err( const int npts, const int order, const double scalefac,
                const double ie, const double ge, const double de )
{
  bool isFailed = false;
  double ierr, gerr, derr;
  calculate_fields( npts, order, scalefac, ierr, gerr, derr );
  if( abs(ie-ierr)/ie > 1e-8 ){
    isFailed=true;
    cout << setprecision(10) << "interpolation failed: " << ie << ", " << ierr << endl;
  }
  if( abs(ge-gerr)/ge > 1e-8 ){
    isFailed=true;
    cout << setprecision(10) << "gradient failed: " << ge << ", " << gerr << endl;
  }
  if( abs(de-derr)/de > 1e-8 ){
    isFailed=true;
    cout << setprecision(10) << "divergence failed: " << de << ", " << derr << endl;
  }
  return isFailed;
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

  bool isFailed = false;

  isFailed = isFailed || check_err( 20, 2, 1.0, 3.0258405e-2, 3.17782874e-2, 1.98599559e-1 );
  isFailed = isFailed || check_err( 80, 2, 4.0, 3.06998703e-3, 1.99207515e-2, 8.06760267e-2 );
  isFailed = isFailed || check_err( 25, 4, 4.0, 0.0223312291, 0.1020888904, 0.4327000552 );

  if( isFailed ) cout << "FAIL" << endl;
  else cout << "PASS" << endl;

}

//====================================================================
