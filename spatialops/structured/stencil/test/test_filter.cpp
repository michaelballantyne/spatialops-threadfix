#include <spatialops/structured/stencil/BoxFilter.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/Grid.h>
#include <spatialops/Nebo.h>

#include <test/TestHelper.h>

#include <vector>

using namespace SpatialOps;
using namespace structured;

template<typename DirT>
bool run_test( const IntVec& dim )
{
  TestHelper status(false);

  const bool bc[3] = {false,false,false};
  const std::vector<double> length(3,1.0);

  MemoryWindow mw = get_window_with_ghost<SVolField>( dim, bc[0], bc[1], bc[2] );

  SVolField    f( mw, NULL );
  SVolField fbar( mw, NULL );

  Grid grid(dim,length);
  grid.set_coord<DirT>(f);

  const BoxFilter<SVolField> filter;
  filter.apply_to_field( f, fbar );

  SVolField::const_interior_iterator i1=f.interior_begin();
  SVolField::const_interior_iterator i1e=f.interior_end();
  SVolField::const_interior_iterator i2=fbar.interior_begin();
  for( ; i1!=i1e; ++i1, ++i2 ){
    const double err = std::abs( *i1-*i2 ) / *i1;
    status( err<1e-15 );
//    std::cout << *i1 << " : " << *i2 << ", err=" << err << std::endl;
  }

  return status.ok();
}

int main()
{
  TestHelper status(true);

  status( run_test<XDIR>( IntVec(30, 1, 1) ), "[30, 1, 1], X" );
  status( run_test<XDIR>( IntVec(30,30, 1) ), "[30,30, 1], X" );
  status( run_test<XDIR>( IntVec(30, 1,30) ), "[30, 1,30], X" );
  status( run_test<XDIR>( IntVec(30,30,30) ), "[30,30,30], X" );

  std::cout << std::endl;

  status( run_test<YDIR>( IntVec( 1,30, 1) ), "[ 1,30, 1], Y" );
  status( run_test<YDIR>( IntVec(30,30, 1) ), "[30,30, 1], Y" );
  status( run_test<YDIR>( IntVec( 1,30,30) ), "[ 1,30,30], Y" );
  status( run_test<YDIR>( IntVec(30,30,30) ), "[30,30,30], Y" );

  std::cout << std::endl;

  status( run_test<ZDIR>( IntVec( 1, 1,30) ), "[ 1, 1,30], Z" );
  status( run_test<ZDIR>( IntVec(30, 1,30) ), "[30, 1,30], Z" );
  status( run_test<ZDIR>( IntVec( 1,30,30) ), "[ 1,30,30], Z" );
  status( run_test<ZDIR>( IntVec(30,30,30) ), "[30,30,30], Z" );

  std::cout << std::endl;

  if( status.ok() ){
    std::cout << "PASS" << std::endl;
    return 0;
  }
  std::cout << "FAIL" << std::endl;
  return -1;
}
