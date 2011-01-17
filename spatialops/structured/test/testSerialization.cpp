#include <test/TestHelper.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>


//#define BINARY_IO
#define ASCII_IO

#ifdef BINARY_IO
# include <boost/archive/binary_oarchive.hpp>
# include <boost/archive/binary_iarchive.hpp>
  typedef boost::archive::binary_oarchive OutputArchive;
  typedef boost::archive::binary_iarchive InputArchive;
#elif defined ASCII_IO
# include <boost/archive/text_oarchive.hpp>
# include <boost/archive/text_iarchive.hpp>
  typedef boost::archive::text_oarchive OutputArchive;
  typedef boost::archive::text_iarchive InputArchive;
#else
# error "NO VALID IO SCHEME DEFINED"
#endif

using namespace SpatialOps;
using namespace structured;

int main()
{
  TestHelper status(true);

  const IntVec npts(5,1,1);
  const MemoryWindow window( get_window_with_ghost<SVolField>(npts,true,true,true) );
  SVolField fx( window, NULL );
  SVolField fy( window, NULL );
  SVolField fz( window, NULL );

  const SVolField::iterator ixend=fx.end();
  SVolField::iterator ix=fx.begin();
  SVolField::iterator iy=fy.begin();
  SVolField::iterator iz=fz.begin();
  for( ; ix!=ixend; ++ix, ++iy, ++iz ){
    *ix = ix.i();
    *iy = iy.j();
    *iz = iz.k();
  }

  {
    std::ofstream fout("test.out",std::ios_base::trunc|std::ios_base::out);
    OutputArchive ar(fout);

    ar << npts << window << fx << fy << fz;
  }
  {
    std::ifstream fin("test.out",std::ios_base::in);
    InputArchive ar(fin);

    IntVec v;
    ar >> v;
    status( v==npts, "IntVec" );

    MemoryWindow w( IntVec(0,0,0) );
    ar >> w;
    status( w==window, "MemoryWindow" );

    const MemoryWindow w0(IntVec(0,0,0));
    SVolField xf( w0, NULL );
    SVolField yf( w0, NULL );
    SVolField zf( w0, NULL );
    ar >> xf >> yf >> zf;

    status( fx==xf, "x field" );
    status( fy==yf, "y field" );
    status( fz==zf, "z field" );
  }

  {
    std::stringstream out;
    OutputArchive ar(out);
    ar << npts << window << fx << fy << fz;

    InputArchive arin(out);
    IntVec n2;
    MemoryWindow w2(n2);
    SVolField xf( w2, NULL );
    SVolField yf( w2, NULL );
    SVolField zf( w2, NULL );

    arin >> n2 >> w2 >> xf >> yf >> zf;

    status( npts==n2, "IntVec - stringstream" );
    status( window==w2, "MemoryWindow - stringstream" );
    status( fx==xf, "x field - stringstream" );
    status( fy==yf, "y field - stringstream" );
    status( fz==zf, "z field - stringstream" );
  }

  if( status.ok() ) return 0;
  return -1;
}
