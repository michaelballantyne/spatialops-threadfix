#include <spatialops/FieldExpressionsExtended.h>
#include <spatialops/FieldReductions.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/Stencil2.h>

#include <spatialops/WriteMatlab.h>

#include "Grid.h"
#include <test/TestHelper.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;

#include <stdexcept>
using std::cout;
using std::endl;

//====================================================================

template< typename FieldT >
void function( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& f )
{
  f <<= sin(x) + cos(y) + sin(z);
}

template< typename CoordT, typename FieldT >
struct FunctionDer{ void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df ); };

template< typename FieldT >
struct FunctionDer<XDIR,FieldT>
{
  static void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df )
  {
    df <<= cos(x);
  }
};
template< typename FieldT >
struct FunctionDer<YDIR,FieldT>
{
  static void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df )
  {
    df <<= -1.0*sin(y);
  }
};
template< typename FieldT >
struct FunctionDer<ZDIR,FieldT>
{
  static void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df )
  {
    df <<= cos(z);
  }
};

template< typename CoordT, typename FieldT > struct Function2Der;

template< typename FieldT > struct Function2Der<XDIR,FieldT>
{
  static void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df )
  {
    df <<= -1.0*sin(x);
  }
};
template< typename FieldT > struct Function2Der<YDIR,FieldT>
{
  static void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df )
  {
    df <<= -1.0*cos(y);
  }
};
template< typename FieldT > struct Function2Der<ZDIR,FieldT>
{
  static void value( const FieldT& x, const FieldT& y, const FieldT& z, FieldT& df )
  {
    df <<= -1.0*sin(z);
  }
};

//====================================================================

template< typename FieldT >
double interior_norm( const FieldT& f1, const FieldT& f2 )
{
  const FieldT f1interior( f1.window_without_ghost(), &f1[0], ExternalStorage );
  const FieldT f2interior( f2.window_without_ghost(), &f2[0], ExternalStorage );
  const double l2 = field_norm( f1interior-f2interior ) / field_norm(f2interior);
//   const double linf = field_max( abs(f1interior-f2interior) );
//   const double l1 = field_sum( abs(f1interior-f2interior) ) / f1.window_without_ghost().npts();
  return l2;
}

//====================================================================

bool check_convergence( const std::vector<double>& spacings,
                        const std::vector<double>& norms,
                        const unsigned int order )
{
  const size_t n = spacings.size();
  std::vector<double> errRatio(n,0), calcOrder(n-1,0);
  for( size_t i=0; i<n; ++i ){
    const double ideal = norms[0] * std::pow( spacings[i] / spacings[0], double(order) );
    errRatio[i] = norms[i] / norms[0];
  }
  for( size_t i=0; i<n-1; ++i ){
    const double num = log10( errRatio[i+1]/errRatio[i] );
    const double den = log10( spacings[i]/spacings[i+1] );
    calcOrder[i] = std::abs(num/den);
  }
  const double maxOrd = *std::max_element( calcOrder.begin(), calcOrder.end() );

  if( maxOrd < 0.99*order ){
    cout << "Order of accuracy: " << maxOrd << endl;
    return false;
  }
  return true;
}

//--------------------------------------------------------------------

template< typename VolT, typename FaceT >
void
apply_stencil( const IntVec& npts, const bool* bcPlus,
               double& interpNorm, double& gradNorm, double& divNorm )
{
  const MemoryWindow vmw = get_window_with_ghost<VolT >( npts, bcPlus[0], bcPlus[1], bcPlus[2] );
  const MemoryWindow fmw = get_window_with_ghost<FaceT>( npts, bcPlus[0], bcPlus[1], bcPlus[2] );

  std::vector<double> length(3,10.0);
  const Grid grid( npts, length );

  VolT   vol( vmw, NULL ), xvol(vmw,NULL), yvol(vmw,NULL), zvol(vmw,NULL);
  FaceT face( fmw, NULL ), xface(fmw,NULL), yface(fmw,NULL), zface(fmw,NULL);

  grid.set_coord<XDIR>( xvol );
  grid.set_coord<YDIR>( yvol );
  grid.set_coord<ZDIR>( zvol );

  grid.set_coord<XDIR>( xface );
  grid.set_coord<YDIR>( yface );
  grid.set_coord<ZDIR>( zface );

  // set the function value
  function( xvol, yvol, zvol, vol );

  FaceT faceExact( fmw, NULL );

  // interpolant
  {
    const Stencil2< Interpolant, VolT, FaceT > interpOp( 0.5, 0.5 );
    interpOp.apply_to_field( vol, face );
    function( xface, yface, zface, faceExact );

    interpNorm = interior_norm( face, faceExact );
//     write_matlab( xvol, "xvol" );
//     write_matlab( vol,  "phiVol" );
//     write_matlab( xface, "xface" );
//     write_matlab( face, "phiFace" );
//     write_matlab( faceExact, "phiFaceExact" );
  }

  // gradient
  {
    const double spc = grid.spacing<typename FaceT::Location::FaceDir>();
    const Stencil2< Gradient, VolT, FaceT > gradOp( -1.0/spc, 1.0/spc );
    gradOp.apply_to_field( vol, face );

    // set the exact gradient value
    FunctionDer<typename FaceT::Location::FaceDir,FaceT>::value( xface, yface, zface, faceExact );

    gradNorm = interior_norm( face, faceExact );
  }

  // divergence
  {
    FunctionDer<typename FaceT::Location::FaceDir,FaceT>::value( xface, yface, zface, faceExact );

    const double spc = grid.spacing<typename FaceT::Location::FaceDir>();
    Stencil2< Divergence, FaceT, VolT > divOp( -1.0/spc, 1.0/spc );
    divOp.apply_to_field( faceExact, vol );

    // set the exact divergence value
    VolT volExact( vmw, NULL );
    Function2Der<typename FaceT::Location::FaceDir,VolT>::value( xvol, yvol, zvol, volExact );
    divNorm = interior_norm( vol, volExact );

//     write_matlab( vol, "div");
//     write_matlab( volExact, "divx" );
  }
}

//--------------------------------------------------------------------

template< typename Vol >
bool
run_variants( const IntVec npts,
              const bool* bcPlus )
{
  TestHelper status(false);

  typedef typename FaceTypes<Vol>::XFace  XFace;
  typedef typename FaceTypes<Vol>::YFace  YFace;
  typedef typename FaceTypes<Vol>::ZFace  ZFace;

  const int nrefine = 5;
  std::vector<double> interpNorms(nrefine,0.0), gradNorms(nrefine,0.0), divNorms(nrefine,0.0);
  std::vector<double> spacings(nrefine,0.0);

  // x-convergence
  if( npts[0]>1 ){
    IntVec n(npts);
    const unsigned int ix = 0;
    for( unsigned int icount=0; icount<nrefine; ++icount ){
      spacings[icount] = 1.0/n[ix];
      apply_stencil<Vol,XFace>( n, bcPlus, interpNorms[icount], gradNorms[icount], divNorms[icount] );
      n[ix] *= 2;
    }
    status( check_convergence( spacings, interpNorms, 2 ), "x-Interpolant" );
    status( check_convergence( spacings, gradNorms,   2 ), "x-Gradient" );
    status( check_convergence( spacings, divNorms,    2 ), "x-Divergence" );
  }

  // y-convergence
  if( npts[1]>1 ){
    IntVec n(npts);
    const unsigned int ix = 1;
    for( unsigned int icount=0; icount<nrefine; ++icount ){
      spacings[icount] = 1.0/n[ix];
      apply_stencil<Vol,YFace>( n, bcPlus, interpNorms[icount], gradNorms[icount], divNorms[icount] );
      n[ix] *= 2;
    }
    status( check_convergence( spacings, interpNorms, 2 ), "y-Interpolant" );
    status( check_convergence( spacings, gradNorms,   2 ), "y-Gradient" );
    status( check_convergence( spacings, divNorms,    2 ), "y-Divergence" );
  }

  // z-convergence
  if( npts[2]>1 ){
    IntVec n(npts);
    const unsigned int ix = 2;
    for( unsigned int icount=0; icount<nrefine; ++icount ){
      spacings[icount] = 1.0/n[ix];
      apply_stencil<Vol,ZFace>( n, bcPlus, interpNorms[icount], gradNorms[icount], divNorms[icount] );
      n[ix] *= 2;
    }
    status( check_convergence( spacings, interpNorms, 2 ), "z-Interpolant" );
    status( check_convergence( spacings, gradNorms,   2 ), "z-Gradient" );
    status( check_convergence( spacings, divNorms,    2 ), "z-Divergence" );
  }

  return status.ok();
}

//--------------------------------------------------------------------

template< typename SrcT, typename DestT >
bool interp_test( const unsigned int dir,
                  IntVec npts,
                  const bool* bcPlus )
{
  const int nrefine = 5;
  std::vector<double> interpNorms(nrefine,0.0);
  std::vector<double> spacings(nrefine,0.0);

  for( unsigned int icount=0; icount<nrefine; ++icount ){

    const MemoryWindow vmw = get_window_with_ghost<SrcT >( npts, bcPlus[0], bcPlus[1], bcPlus[2] );
    const MemoryWindow fmw = get_window_with_ghost<DestT>( npts, bcPlus[0], bcPlus[1], bcPlus[2] );

    std::vector<double> length(3,10.0);
    const Grid grid( npts, length );

    SrcT   src( vmw, NULL ), xsrc(vmw,NULL), ysrc(vmw,NULL), zsrc(vmw,NULL);
    DestT dest( fmw, NULL ), xdest(fmw,NULL), ydest(fmw,NULL), zdest(fmw,NULL);

    grid.set_coord<XDIR>( xsrc );
    grid.set_coord<YDIR>( ysrc );
    grid.set_coord<ZDIR>( zsrc );

    grid.set_coord<XDIR>( xdest );
    grid.set_coord<YDIR>( ydest );
    grid.set_coord<ZDIR>( zdest );

    // set the function value
    function( xsrc, ysrc, zsrc, src );

    DestT destExact( fmw, NULL );

    // interpolant
    const Stencil2< Interpolant, SrcT, DestT > interpOp( 0.5, 0.5 );
    interpOp.apply_to_field( src, dest );
    function( xdest, ydest, zdest, destExact );

    interpNorms[icount] = interior_norm( dest, destExact );
    spacings[icount] = 1.0/npts[dir];
    npts[dir] *= 2;

    write_matlab( xdest, "x" );
    write_matlab( destExact, "phix");
    write_matlab( dest, "phi");
  }
  return check_convergence( spacings, interpNorms, 2 );
}

//--------------------------------------------------------------------

int main( int iarg, char* carg[] )
{
  int nx, ny, nz;
  bool bcplus[] = { false, false, false };

  {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message\n" )
      ( "nx",   po::value<int>(&nx)->default_value(5), "number of points in x-dir for base mesh" )
      ( "ny",   po::value<int>(&ny)->default_value(5), "number of points in y-dir for base mesh" )
      ( "nz",   po::value<int>(&nz)->default_value(5), "number of points in z-dir for base mesh" )
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
    std::string bcx = bcplus[0] ? "ON" : "OFF";
    std::string bcy = bcplus[1] ? "ON" : "OFF";
    std::string bcz = bcplus[2] ? "ON" : "OFF";
    cout << "Run information: " << endl
         << "  bcx    : " << bcx << endl
         << "  bcy    : " << bcy << endl
         << "  bcz    : " << bcz << endl
         << "  domain : " << npts << endl
         << endl;
  }

  try{
    status( run_variants< SVolField >( npts, bcplus ), "SVol operators" );
    if( npts[0] > 1 ) status( run_variants< XVolField >( npts, bcplus ), "XVol operators" );
    if( npts[1] > 1 ) status( run_variants< YVolField >( npts, bcplus ), "YVol operators" );
    if( npts[2] > 1 ) status( run_variants< ZVolField >( npts, bcplus ), "ZVol operators" );

    status( interp_test< XVolField, YSurfXField >( 1, npts, bcplus ), "InterpXVolYSurfX" );
    status( interp_test< XVolField, ZSurfXField >( 2, npts, bcplus ), "InterpXVolZSurfX" );

    status( interp_test< YVolField, XSurfYField >( 0, npts, bcplus ), "InterpYVolXSurfY" );
    status( interp_test< YVolField, ZSurfYField >( 2, npts, bcplus ), "InterpYVolZSurfY" );

    status( interp_test< ZVolField, XSurfZField >( 0, npts, bcplus ), "InterpZVolXSurfZ" );
    status( interp_test< ZVolField, YSurfZField >( 1, npts, bcplus ), "InterpZVolYSurfZ" );

    status( interp_test< SVolField, XVolField >( 0, npts, bcplus ), "InterpSVolXVol" );
    status( interp_test< SVolField, YVolField >( 1, npts, bcplus ), "InterpSVolYVol" );
    status( interp_test< SVolField, ZVolField >( 2, npts, bcplus ), "InterpSVolZVol" );
  }
  catch( std::runtime_error& e ){
    cout << e.what() << endl;
    return -1;
  }
  if( status.ok() ) return 0;
  return -1;
}

//--------------------------------------------------------------------

