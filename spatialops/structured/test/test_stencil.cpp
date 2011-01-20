#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/Stencil.h>

#include <spatialops/WriteMatlab.h>

using namespace SpatialOps;
using namespace structured;

//--------------------------------------------------------------------

template< typename OpT, typename SrcT, typename DestT >
void
apply_stencil( const Stencil2<OpT,SrcT,DestT>& op, const IntVec& npts )
{
  const bool bcx=true;
  const bool bcy=true;
  const bool bcz=true;
//   const bool bcx=false;
//   const bool bcy=false;
//   const bool bcz=false;

  const MemoryWindow smw = get_window_with_ghost<SrcT >( npts, bcx, bcy, bcz );
  const MemoryWindow dmw = get_window_with_ghost<DestT>( npts, bcx, bcy, bcz );

  SrcT   src( smw, NULL );
  DestT dest( dmw, NULL );

  const double dx = 1.0;
  const double dy = 1.0;
  const double dz = 1.0;

  const int ilo = smw.offset(0);
  const int ihi = ilo + smw.extent(0);

  const int jlo = smw.offset(1);
  const int jhi = jlo + smw.extent(1);

  const int klo = smw.offset(2);
  const int khi = klo + smw.extent(2);

  for( int k=klo; k!=khi; ++k ){
    for( int j=jlo; j!=jhi; ++j ){
      for( int i=ilo; i!=ihi; ++i ){
        const double x = i*dx;
        const double y = j*dy;
        const double z = k*dz;
        src[ smw.flat_index( IntVec(i,j,k) ) ] = x*x + y*y + z*z;
      }
    }
  }

  op.apply_to_field( src, dest );

//   write_matlab( src,  "src"  );
//   write_matlab( dest, "dest" );
}

//--------------------------------------------------------------------

template< typename Vol >
void
run_variants( const IntVec& npts )
{
  typedef typename FaceTypes<Vol>::XFace  XFace;
  typedef typename FaceTypes<Vol>::YFace  YFace;
  typedef typename FaceTypes<Vol>::ZFace  ZFace;

  apply_stencil< Interpolant, Vol, XFace >( Stencil2<Interpolant,Vol,XFace>(0.5,0.5), npts );
  apply_stencil< Interpolant, Vol, YFace >( Stencil2<Interpolant,Vol,YFace>(0.5,0.5), npts );
  apply_stencil< Interpolant, Vol, ZFace >( Stencil2<Interpolant,Vol,ZFace>(0.5,0.5), npts );

  apply_stencil< Gradient, Vol, XFace >( Stencil2<Gradient,Vol,XFace>(-1,1), npts );
  apply_stencil< Gradient, Vol, YFace >( Stencil2<Gradient,Vol,YFace>(-1,1), npts );
  apply_stencil< Gradient, Vol, ZFace >( Stencil2<Gradient,Vol,ZFace>(-1,1), npts );

  apply_stencil< Divergence, XFace, Vol >( Stencil2<Divergence,XFace,Vol>(-1,1), npts );
  apply_stencil< Divergence, YFace, Vol >( Stencil2<Divergence,YFace,Vol>(-1,1), npts );
  apply_stencil< Divergence, ZFace, Vol >( Stencil2<Divergence,ZFace,Vol>(-1,1), npts );
}

//--------------------------------------------------------------------

int main()
{
  const IntVec npts( 80, 100, 70 );

  run_variants< SVolField >( npts );
  run_variants< XVolField >( npts );
  run_variants< YVolField >( npts );
  run_variants< ZVolField >( npts );

//   apply_stencil< Interpolant, XVolField, SSurfXField >( StencilNull<Interpolant,XVolField,SSurfXField>(), npts );
}

//--------------------------------------------------------------------

