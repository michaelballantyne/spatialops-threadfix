#include <vector>
#include <iostream>


#include <spatialops/FVStaggered.h>
#include <spatialops/FVTopHatFilter.h>

#include <Grid.h>
#include <Functions.h>
#include <buildOps.h>

using namespace SpatialOps;
using namespace FVStaggered;

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

void manual_filter( const std::vector<int>& dim,
                    const int npts,
                    const SVolField& f,
                    SVolField& fbar )
{
  const int ng = SVolField::Ghost::NM;
  const int nside = npts/2;

  const int nx = get_nx<SVolField>( dim, true );
  const int ny = get_ny<SVolField>( dim, true );
  const int nz = get_nz<SVolField>( dim, true );

  fbar = 0.0;

  for( int ipt=0; ipt<f.get_ntotal(); ++ipt ){

    const IndexTriplet ijk = flat2ijk<SVolField>::value( dim, ipt, true, true, true );

    std::vector<int> half(3,0);
    half[0] = std::min( std::min( nside, ijk.i ), std::min( nside, nx-ijk.i-1 ) );
    half[1] = std::min( std::min( nside, ijk.j ), std::min( nside, ny-ijk.j-1 ) );
    half[2] = std::min( std::min( nside, ijk.k ), std::min( nside, nz-ijk.k-1 ) );

    double& fb = fbar[ipt];

    int n=0;
    for( int i=ijk.i-half[0]; i<=ijk.i+half[0]; ++i ){
      for( int j=ijk.j-half[1]; j<=ijk.j+half[1]; ++j ){
        for( int k=ijk.k-half[2]; k<=ijk.k+half[2]; ++k ){
          const int ixf = ijk2flat<SVolField>::value( dim, IndexTriplet( i,j,k ), true, true, true );
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

  vector<int> dim(3,1);
  dim[0] = 30;
  dim[1] = 15;
  dim[2] = 20;

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

  SVolField f   ( get_n_tot    <SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                  get_ghost_set<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                  NULL );
  SVolField fbar( get_n_tot    <SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                  get_ghost_set<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                  NULL );

  SVolField fbar2( get_n_tot    <SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                   get_ghost_set<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                   NULL );

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
