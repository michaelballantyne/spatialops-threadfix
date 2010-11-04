#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <test/TestHelper.h>

#include <sstream>
#include <fstream>

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;


void test_iterator( const IntVec npts,
                    TestHelper& overall )
{
  TestHelper status(false);

  const MemoryWindow window( get_window_with_ghost<SVolField>(npts,true,true,true) );
  SVolField f1( window, NULL );
  SVolField f2( window, NULL );
  f1 = 2.0;
  f2 = 1.0;

  SVolField::iterator if2=f2.begin();
  const SVolField::iterator if2e=f2.end();
  SVolField::const_iterator if1=f1.begin();
  for( ; if2!=if2e; ++if1, ++if2 ){
    *if2 += *if1;
  }

  if2 = f2.begin() + 2;
  status( f2[2] == *if2, "iterator + operator" );
  status( &f2[2] == &(*if2), "iterator + operator address" );

  if2 += 3;
  status( &f2[5] == &(*if2), "iterator += address" );

  const size_t ng = SVolField::Ghost::NGHOST;
  const size_t ihi = npts[0]>1 ? npts[0] + 2*ng : 1;
  const size_t jhi = npts[1]>1 ? npts[1] + 2*ng : 1;
  const size_t khi = npts[2]>1 ? npts[2] + 2*ng : 1;

  if1=f1.begin();
  if2=f2.begin();
  for( size_t k=0; k<khi; ++k ){
    for( size_t j=0; j<jhi; ++j ){
      for( size_t i=0; i<ihi; ++i ){
        {
          std::ostringstream msg;
          msg << "test_iterator 1.1: [" << i << "," << j << "," << k << "],  found: " << f2(i,j,k) << ", expected: 3.0";
          status( f2(i,j,k) == 3.0, msg.str() );
        }
        {
          std::ostringstream msg;
          msg << "test_iterator mem check f1: [" << i << "," << j << "," << k << "]";
          const double& f1pt = f1(i,j,k);
          status( &f1pt == &*if1, msg.str() );
        }
        {
          std::ostringstream msg;
          msg << "test_iterator mem check f2: [" << i << "," << j << "," << k << "]";
          const double& f2pt = f2(i,j,k);
          status( &f2pt == &*if2, msg.str() );
        }
        ++if1; ++if2;
      }
    }
  }

  status( if1 == f1.end(), "iterator end (1)" );
  status( if2 == f2.end(), "iterator end (2)" );

  f1 = 2.0;
  f2 = 1.0;
  f2 += f1;
  for( size_t k=0; k<khi; ++k ){
    for( size_t j=0; j<jhi; ++j ){
      for( size_t i=0; i<ihi; ++i ){
        std::ostringstream msg;
        msg << "test_iterator 2: [" << i << "," << j << "," << k << "],  found: " << f2(i,j,k) << ", expected: 3.0";
        status( f2(i,j,k) == 3.0, msg.str() );
      }
    }
  }

  //  status.report_status(true);
  SVolField::iterator iter=f1.end()-1;
  for( size_t k=khi; k>0; --k ){
    for( size_t j=jhi; j>0; --j ){
      for( size_t i=ihi; i>0; --i ){
        std::ostringstream msg;
        msg << "test_iterator backward iterator mem check f1: [" << i-1 << "," << j-1 << "," << k-1 << "]";
        status( &f1(i-1,j-1,k-1) == &*iter, msg.str() );
        --iter;
      }
    }
  }

  std::ostringstream msg;
  msg << "full iterator test " << npts[0] << "x" << npts[1] << "x" << npts[2];
  overall( status.ok(), msg.str() );
}

//--------------------------------------------------------------------

void test_interior( const IntVec npts,
                    TestHelper& overall )
{
  const MemoryWindow window( get_window_with_ghost<SVolField>(npts,true,true,true) );
  SVolField f1( window, NULL );
  SVolField f2( window, NULL );
  f1 = 2.0;

  const int ng = SVolField::Ghost::NGHOST;
  const size_t ilo = npts[0]>1 ? ng : 0;
  const size_t jlo = npts[1]>1 ? ng : 0;
  const size_t klo = npts[2]>1 ? ng : 0;
  const size_t ihi = npts[0]>1 ? npts[0]+ng : 1;
  const size_t jhi = npts[1]>1 ? npts[1]+ng : 1;
  const size_t khi = npts[2]>1 ? npts[2]+ng : 1;

  f2 = 0.0;
  // set interior values
  for( size_t k=klo; k<khi; ++k ){
    for( size_t j=jlo; j<jhi; ++j ){
      for( size_t i=ilo; i<ihi; ++i ){
        f2(i,j,k) = 1+i+j+k;
      }
    }
  }

  TestHelper status(false);

  SVolField::interior_iterator if2=f2.interior_begin();
  const SVolField::interior_iterator if2e=f2.interior_end();
  SVolField::const_interior_iterator if1=f1.interior_begin();
  const SVolField::interior_iterator if1e=f1.interior_end();
  for( size_t k=klo; k<khi; ++k ){
    for( size_t j=jlo; j<jhi; ++j ){
      for( size_t i=ilo; i<ihi; ++i ){
        {
          const double& f2pt = f2(i,j,k);
          const double* f2pti = &*if2;
          std::ostringstream msg;  msg << "f2 mem loc at " << i<<","<<j<<","<<k;
          status( &f2pt == f2pti, msg.str() );
        }
        {
          const double& f1pt = f1(i,j,k);
          const double* f1pti = &*if1;
          std::ostringstream msg;  msg << "f1 mem loc at " << i<<","<<j<<","<<k;
          status( &f1pt == f1pti, msg.str() );
        }
        ++if2; ++if1;
      }
    }
  }

  if1 = f1.interior_begin();
  if2 = f2.interior_begin();
  for( ; if2!=if2e; ++if1, ++if2 ){
    *if2 += *if1;
  }

  for( size_t k=klo; k<khi; ++k ){
    for( size_t j=jlo; j<jhi; ++j ){
      for( size_t i=ilo; i<ihi; ++i ){
        const double val = 1+i+j+k + 2.0;
        std::ostringstream msg;  msg << i<<","<<j<<","<<k << ",  found " << f2(i,j,k) << ", expected " << val;
        status( f2(i,j,k) == val, msg.str() );
      }
    }
  }
  std::ostringstream msg;
  msg << "interior iterator test " << npts[0] << "x" << npts[1] << "x" << npts[2] << std::endl;
  overall( status.ok(), msg.str() );
}

//--------------------------------------------------------------------

int main()
{
  TestHelper overall(true);

  // test iterators
  test_iterator( IntVec(3,3,3), overall );
  test_iterator( IntVec(3,4,1), overall );
  test_iterator( IntVec(4,3,2), overall );

  std::cout << std::endl << std::endl;

  test_interior( IntVec(3,3,3), overall );

  // test basic layout and operators
  {
    TestHelper status(false);

    const int npts[3] = {10,11,12};
    const MemoryWindow window(npts);
    SVolField svol1( window, NULL, InternalStorage );
    SVolField svol2( window, NULL, InternalStorage );

    for( int k=0; k<npts[2]; ++k ){
      for( int j=0; j<npts[1]; ++j ){
        for( int i=0; i<npts[0]; ++i ){
          svol1(i,j,k) = i + j + k;
        }
      }
    }

    svol2 = 2.0;
    svol1 += svol2;
    svol1 *= svol2;
    svol1 /= svol2;
    svol1 -= svol2;

    for( int k=0; k<npts[2]; ++k ){
      for( int j=0; j<npts[1]; ++j ){
        for( int i=0; i<npts[0]; ++i ){
          const double ans = (i + j + k);
          std::ostringstream msg;  msg << "("<<i<<","<<j<<","<<k<<")" << ",  found " << svol1(i,j,k) << ", expected " << ans;
          status( ans==svol1(i,j,k), msg.str() );
        }
      }
    }

    {
      SpatialOps::SpatFldPtr<SVolField> sv3 = SpatialOps::SpatialFieldStore<SVolField>::self().get( svol1 );
      *sv3 = svol1;
      status( *sv3 == svol1, "spatial field pointer from store" );
    }

    overall( status.ok(), "field operations" );
  }
  if( overall.isfailed() ) return -1;
  return 0;
}
