#include <vector>
#include <iostream>


#include <spatialops/structured/FVStaggered.h>

#include "Grid.h"
#include "Functions.h"
#include "buildOps.h"

using namespace SpatialOps;
using namespace structured;

//====================================================================

void build_filter( const int filterWidth,
                   const std::vector<int>& dim,
                   const std::vector<bool>& bcPlus,
                   OperatorDatabase& opDB )
{
  FilterSVol::Assembler tfa( filterWidth, dim, bcPlus[0], bcPlus[1], bcPlus[2] );
  opDB.register_new_operator<FilterSVol>( new FilterSVol(tfa) );
}

//====================================================================

void manual_filter( const IntVec& dim,
                    const int npts,
                    const SVolField& f,
                    SVolField& fbar )
{
  const int nside = npts/2;

  const int nx = get_nx_with_ghost<SVolField>( dim[0], true );
  const int ny = get_ny_with_ghost<SVolField>( dim[1], true );
  const int nz = get_nz_with_ghost<SVolField>( dim[2], true );

  fbar = 0.0;
  // jcs this will only work for the trivial memory window.
  SVolField::const_iterator ifval=f.begin(), ifvale=f.end();
  SVolField::iterator ifbar=fbar.begin();
  for( size_t ipt=0; ifval!=ifvale; ++ifval, ++ifbar, ++ipt ){

    const IntVec ijk = flat2ijk<SVolField>::value( dim, ipt, true, true, true );

    std::vector<int> half(3,0);
    half[0] = std::min( std::min( nside, ijk.i ), std::min( nside, nx-ijk.i-1 ) );
    half[1] = std::min( std::min( nside, ijk.j ), std::min( nside, ny-ijk.j-1 ) );
    half[2] = std::min( std::min( nside, ijk.k ), std::min( nside, nz-ijk.k-1 ) );

    double& fb = *ifbar;

    int n=0;
    for( int i=ijk.i-half[0]; i<=ijk.i+half[0]; ++i ){
      for( int j=ijk.j-half[1]; j<=ijk.j+half[1]; ++j ){
        for( int k=ijk.k-half[2]; k<=ijk.k+half[2]; ++k ){
          const int ixf = ijk2flat<SVolField>::value( dim, IntVec( i,j,k ), true, true, true );
          fb += f[ixf];
          ++n;
        }
      }
    }
    fb /= double(n);
  }

}

//====================================================================

struct abs_comparison{
  bool operator()( const double x1, const double x2 )
  {
    return std::fabs(x1) < std::fabs(x2);
  }
};

//====================================================================

int main()
{
  using namespace std;

  IntVec dim(30,15,20);

//   cout << "interior nx = "; cin >> dim[0];
//   cout << "interior ny = "; cin >> dim[1];
//   cout << "interior nz = "; cin >> dim[2];
//   cout << endl;

  const int filterWidth = 5;

  vector<double> length(3,6);
  vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    if( dim[i]>1 ) spacing[i] = length[i]/dim[i];
  }

  vector<bool> bcFlag(3,true);

  OperatorDatabase opDB;
  build_filter( filterWidth, dim, bcFlag, opDB );
  build_ops( dim, spacing, bcFlag, opDB );  // need this to form the mesh - stupid for this test case...
  Grid grid( dim, spacing, bcFlag, opDB );

  const MemoryWindow svolWindow( get_window_with_ghost<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]) );
  SVolField f    ( svolWindow, NULL );
  SVolField fbar ( svolWindow, NULL );
  SVolField fbar2( svolWindow, NULL );

  SinFun<SVolField> fun( grid.xcoord_svol(), grid.ycoord_svol(), grid.zcoord_svol() );

  fun.evaluate( f );

  const FilterSVol* const filter = opDB.retrieve_operator<FilterSVol>();
  filter->apply_to_field( f, fbar );

  manual_filter( dim, filterWidth, f, fbar2 );

//   f.write_matlab( "f", true );
//   fbar.write_matlab( "fbar", true );
//   fbar2.write_matlab( "fbar2", true );
//   grid.xcoord_svol().write_matlab("x", true);
//   grid.ycoord_svol().write_matlab("y", true);
//   grid.zcoord_svol().write_matlab("z", true);

  f = fbar2-fbar;
//   f.write_matlab("err",true);

  const double maxerr = *std::max_element( f.begin(), f.end(), abs_comparison() );

  if( maxerr > 1e-12 )
    return 1;
  else
    return 0;
}

//====================================================================
