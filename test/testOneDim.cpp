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

  const int nGhost = SVolField::Ghost::NM;

  InterpAssembler Rxa(order,nGhost,xv,xs);
  GradAssembler   Gxa(order,nGhost,xv,xs);
  DivAssembler    Dxa(nGhost,xv,xs);

  const Interp Rx(Rxa);
  const Grad   Gx(Gxa);
  const Div    Dx(Dxa);

  //
  // build fields
  // 
  vector<int> dim(3,1); dim[0]=npts;
  const int ntotSVol  = get_n_tot<SVolField  >(dim,true,true,true);
  const int ntotSSurf = get_n_tot<SSurfXField>(dim,true,true,true);
  const set<size_t> ghostSVol  = get_ghost_set<SVolField  >(dim,true,true,true);
  const set<size_t> ghostSSurf = get_ghost_set<SSurfXField>(dim,true,true,true);

  SVolField          f( ntotSVol,  ghostSVol,  NULL );
  SSurfXField  fxexact( ntotSSurf, ghostSSurf, NULL );
  SSurfXField fxinterp( ntotSSurf, ghostSSurf, NULL );
  SSurfXField    fgrad( ntotSSurf, ghostSSurf, NULL );
  SSurfXField  fgexact( ntotSSurf, ghostSSurf, NULL );
  SVolField   d2fexact( ntotSVol,  ghostSVol,  NULL );
  SVolField        d2f( ntotSVol,  ghostSVol,  NULL );

  const SVolField    xvol( ntotSVol,  ghostSVol,  &get_x_src(nGhost,xv,xs)[0] );
  const SSurfXField xsurf( ntotSSurf, ghostSSurf, &get_x_dest(nGhost,xv,xs)[0] );
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

bool check_err( const int npts, const int order, const double scalefac,
                const double ie, const double ge, const double de,
                const double fcie, const double fcge, const double ccge )
{
  bool isFailed = false;
  double ierr, gerr, derr, fcierr, fcgerr, ccgerr;
  cout << "running test for " << npts << " points and polynomial order " << order
       << " with stretch factor " << scalefac << " ... " << flush;
  calculate_fields( npts, order, scalefac, ierr, gerr, derr, fcierr, fcgerr, ccgerr );
  if( abs(ie-ierr)/ie > 1e-8 ){
    isFailed=true;
    cout << "FAIL" << endl
         << setprecision(10) << "  interpolation failed: " << ie << ", " << ierr << endl;
  }
  if( abs(ge-gerr)/ge > 1e-8 ){
    isFailed=true;
    cout << "FAIL" << endl
         << setprecision(10) << "  gradient failed: " << ge << ", " << gerr << endl;
  }
  if( abs(de-derr)/de > 1e-8 ){
    isFailed=true;
    cout << "FAIL" << endl
         << setprecision(10) << "  divergence failed: " << de << ", " << derr << endl;
  }
  if( abs(fcie-fcierr)/fcie > 1e-8 ){
    isFailed = true;
    cout << "FAIL" << endl
         << setprecision(10) << "  face->cell interp failed: " << fcie << ", " << fcierr << endl;
  }
  if( abs(fcge-fcgerr)/fcge > 1e-8 ){
    isFailed = true;
    cout << "FAIL" << endl
         << setprecision(10) << "  face->cell grad failed: " << fcge << ", " << fcgerr << endl;
  }
  if( abs(ccge-ccgerr)/ccge>1e-8 ){
    isFailed = true;
    cout << "FAIL" << endl
         << setprecision(10) << "  cell->cell grad failed: " << ccge << ", " << ccgerr << endl;
  }
  if( !isFailed ) cout << "PASS" << endl;
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

  isFailed |= check_err( 20, 2, 1.0, 7.3502968e-3, 3.17782874e-2, 1.98599559e-1, 7.410158999e-3, 3.170923728e-2, 0.1256835881   );
  isFailed |= check_err( 40, 2, 1.0, 9.3492781e-4, 7.96273079e-3, 4.99407572e-2, 9.352071377e-4, 7.962443016e-3, 0.03177713888  );
  isFailed |= check_err( 80, 2, 1.0, 1.1722637e-4, 1.99181917e-3, 1.25084456e-2, 1.172462558e-4, 1.991801181e-3, 0.007962658853 );

  isFailed |= check_err( 20, 2, 4.0, 1.44722118e-1, 2.245428298e-1, 1.2583794e-0, 1.037205182e-1, 2.407903611e-1, 0.682537828  );
  isFailed |= check_err( 40, 2, 4.0, 2.27238703e-2, 7.54298446e-2, 3.33347057e-1, 2.234508564e-2, 6.945968971e-2, 0.1936764955 );
  isFailed |= check_err( 80, 2, 4.0, 3.06998703e-3, 1.99207515e-2, 8.06760267e-2, 3.144696195e-3, 1.777579115e-2, 0.0659789562 );

  isFailed |= check_err( 20, 4, 1.0, 7.6162859e-4, 6.591328897e-3, 1.024802205e-1, 3.328861949e-4, 8.587804659e-4, 0.0088222634 );
  isFailed |= check_err( 40, 4, 1.0, 2.45591292e-5, 4.205666662e-4, 2.517237273e-2, 1.061089627e-5, 5.437630461e-5, 0.0005737300688 );
  isFailed |= check_err( 80, 4, 1.0, 7.7041488e-7, 2.658622206e-5, 6.266903321e-3, 3.342428068e-7, 3.407858931e-6, 3.582233267e-05 );
  
  isFailed |= check_err( 20, 4, 4.0, 0.03979495377, 0.1674880963 , 0.6706331688  , 0.02826635825 , 0.04905237428 , 0.1374255994 );
  isFailed |= check_err( 40, 4, 4.0, 3.96734470e-3, 0.02517044905, 0.1514529967  , 0.001911244633, 0.004195780012, 0.05340286527 );
  isFailed |= check_err( 80, 4, 4.0, 1.77251488e-4, 2.06917947e-3, 3.545939321e-2, 7.762276848e-5, 2.846857322e-4, 0.01794738432 );
//   isFailed |= check_err(160, 4, 4.0, 6.18109727e-6, 1.39831006e-4, 8.837083831e-3, 2.645839646e-6, 1.811638886e-5, 0.008605711392 );

  return (int) isFailed;
}

//====================================================================
