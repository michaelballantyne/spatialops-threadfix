#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/Stencil2.h>

#include <spatialops/FieldOperationDefinitions.h>

#include <spatialops/WriteMatlab.h>

#include "Grid.h"

using namespace SpatialOps;
using namespace structured;

#include <stdexcept>
using std::cout;
using std::endl;

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
    df <<= -1.0*cos(z);
  }
};

//--------------------------------------------------------------------

template< typename VolT, typename FaceT >
void
apply_stencil( const IntVec& npts, const bool* bcPlus )
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

    // interpNorm = norm( face-faceExact );
  }

  // gradient
  {
    const double spc = grid.spacing<typename FaceT::Location::FaceDir>();
    const Stencil2< Gradient, VolT, FaceT > gradOp( -1.0/spc, 1.0/spc );
    gradOp.apply_to_field( vol, face );

    // set the exact gradient value
    FunctionDer<typename FaceT::Location::FaceDir,FaceT>::value( xface, yface, zface, faceExact );

    // gradNorm = norm( face-faceExact );
  }

  // divergence
  {
    const double spc = grid.spacing<typename FaceT::Location::FaceDir>();
    Stencil2< Divergence, FaceT, VolT > divOp( -1.0/spc, 1.0/spc );

    divOp.apply_to_field( faceExact, vol );

    // set the exact gradient value
    VolT volExact( vmw, NULL );
    Function2Der<typename FaceT::Location::FaceDir,VolT>::value( xvol, yvol, zvol, volExact );
    // divNorm = norm( vol-volExact );
  }
}

//--------------------------------------------------------------------

template< typename Vol >
void
run_variants( const IntVec& npts, const bool* bcPlus )
{
  typedef typename FaceTypes<Vol>::XFace  XFace;
  typedef typename FaceTypes<Vol>::YFace  YFace;
  typedef typename FaceTypes<Vol>::ZFace  ZFace;

  apply_stencil< Vol, XFace >( npts, bcPlus );
  apply_stencil< Vol, YFace >( npts, bcPlus );
  apply_stencil< Vol, ZFace >( npts, bcPlus );
}

//--------------------------------------------------------------------

int main()
{
  const IntVec npts( 3, 3, 1 );
  bool bcplus[] = { true, true, true };

  try{
    run_variants< SVolField >( npts, bcplus );
    run_variants< XVolField >( npts, bcplus );
    run_variants< YVolField >( npts, bcplus );
    run_variants< ZVolField >( npts, bcplus );

//   apply_stencil< Interpolant, XVolField, SSurfXField >( StencilNull<Interpolant,XVolField,SSurfXField>(), npts );
  }
  catch( std::runtime_error& e ){
    cout << e.what() << endl;
    return -1;
  }
  return 0;
}

//--------------------------------------------------------------------

