#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/Nebo.h>
#include <test/TestHelper.h>
#include <test/FieldHelper.h>

#ifdef ENABLE_THREADS
#include <spatialops/ThreadPool.h>
#include <boost/thread/barrier.hpp>
#endif

#include <cmath>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <sstream>
#include <fstream>

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;

void jcs_pause()
{
  std::cout << "Press <ENTER> to continue...";
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

//--------------------------------------------------------------------
#define COMPAREFIELDS( F1MEMTYPE, F1VAL, F2MEMTYPE, F2VAL, ERROR, MESSAGE)  \
{                                                                           \
    FieldT F1(window, NULL, InternalStorage, F1MEMTYPE);                    \
    FieldT F2(window, NULL, InternalStorage, F2MEMTYPE);                    \
    F1 <<= F1VAL;                                                           \
    F2 <<= F2VAL;                                                           \
    status(manual_compare(F1, F2, ERROR, verboseOutput), MESSAGE);          \
}


template<typename FieldT>
bool manual_compare(FieldT& f1,
                    FieldT& f2,
                    const double error,
                    const bool verboseOutput)
{
  //copy the fields to local ram if applicable
  if(f1.memory_device_type() == EXTERNAL_CUDA_GPU) {
    f1.add_consumer(LOCAL_RAM, 0);
  }
  if(f2.memory_device_type() == EXTERNAL_CUDA_GPU) {
    f2.add_consumer(LOCAL_RAM, 0);
  }

  //iterate through fields.
  typename FieldT::const_iterator if1 = const_cast<const FieldT&>(f1).begin();
  typename FieldT::const_iterator if1e = const_cast<const FieldT&>(f1).end();
  typename FieldT::const_iterator if2 = const_cast<const FieldT&>(f2).begin();

  //manually compare equality
  bool man_equal = true;
  for(; if1 != if1e; ++if1, ++if2) {
    if( std::abs(*if1 - *if2)/(double)std::abs(*if1) > error ) {
      man_equal = false;
      break;
    }
  }

  return man_equal == f1.field_equal(f2, error);
}

template<typename FieldT>
bool test_field_equal( const IntVec npts,
                       const MemoryType memType1,
                       const MemoryType memType2,
                       const bool verboseOutput)
{
  TestHelper status(verboseOutput);
  const MemoryWindow window(npts);
  const int total = npts[0] * npts[1] * npts[2];

  FieldT* f1;
  FieldT* f2;

  //local fields
  FieldT lf1(window, NULL, InternalStorage); initialize_field(lf1, 0);
  FieldT lf2(window, NULL, InternalStorage); initialize_field(lf2, 0); 
  FieldT lf3(window, NULL, InternalStorage); initialize_field(lf3, total); 
#ifdef __CUDACC__
  //gpu fields
  FieldT gf1(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  FieldT gf2(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  FieldT gf3(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  
  //move local initialized fields to gpu if necessary
  if(memType1 == EXTERNAL_CUDA_GPU) {
    lf1.add_consumer(EXTERNAL_CUDA_GPU, 0); 
    gf1 <<= lf1; 
    f1 = &gf1;
  } 
  else {
    f1 = &lf1;
  }

  if(memType2 == EXTERNAL_CUDA_GPU) {
    lf2.add_consumer(EXTERNAL_CUDA_GPU, 0); 
    gf2 <<= lf2;
    f2 = &gf2;
  } 
  else {
    f2 = &lf2;
  }
#else
  f1 = &lf1;
  f2 = &lf2;
#endif

  //two duplicate fields exactly equal
  status(manual_compare(*f1, *f2, 0.0, verboseOutput), "Duplicate Fields Equal");

  //change second field and compare not equal
#ifdef __CUDACC__
  if(memType2 == EXTERNAL_CUDA_GPU) {
    lf3.add_consumer(EXTERNAL_CUDA_GPU, 0); 
    gf3 <<= lf3;
    f2 = &gf2;
  } 
  else {
    f2 = &lf3;
  }
#endif
  status(manual_compare(*f1, *f2, 0.0, verboseOutput), "Non-Duplicate Fields Not Equal");


  //relative error test
  COMPAREFIELDS(memType1, 3.0, memType2, 7.49, 1.5, "Off By 150% (Equal)");
  COMPAREFIELDS(memType1, 3.0, memType2, 7.51, 1.5, "Off By 150% (Not Equal)");
  COMPAREFIELDS(memType1, 3.0, memType2, 3.98, 1.5, "Off By 33% (Equal)");
  COMPAREFIELDS(memType1, 3.0, memType2, 4.10, 1.5, "Off By 33% (Not Equal)");
  COMPAREFIELDS(memType1, 3.0, memType2, 2.97, 1.5, "Off By 1% (Equal)");
  COMPAREFIELDS(memType1, 3.0, memType2, 2.96, 1.5, "Off By 1% (Not Equal)");

  return status.ok();
}

//--------------------------------------------------------------------

template<typename FieldT>
bool test_iterator( const IntVec npts,
                    const bool verboseOutput )
{
  TestHelper status(verboseOutput);

  const MemoryWindow window( get_window_with_ghost<FieldT>(npts,true,true,true) );
  FieldT f1( window, NULL );
  FieldT f2( window, NULL );
  f1 <<= 2.0;
  f2 <<= 1.0;

  typename FieldT::iterator if2=f2.begin();
  const typename FieldT::iterator if2e=f2.end();
  typename FieldT::const_iterator if1=f1.begin();
  for( ; if2!=if2e; ++if1, ++if2 ){
    *if2 += *if1;
  }

  {
    typename FieldT::iterator i1, i2;
    i1 = f1.begin();
    i2 = f1.begin();
    status( *i1 == *i2 );
  }
  {
    typename FieldT::const_iterator i1, i2;
    i1 = f1.begin();
    i2 = f1.begin();
    status( *i1 == *i2 );
  }

  if2 = f2.begin() + 2;
  status( f2[2] == *if2, "iterator + operator" );
  status( &f2[2] == &(*if2), "iterator + operator address" );

  if2 += 3;
  status( &f2[5] == &(*if2), "iterator += address" );

  const IntVec& hi = window.glob_dim();

  if1=f1.begin();
  if2=f2.begin();
  for( int k=0; k<hi[2]; ++k ){
    for( int j=0; j<hi[1]; ++j ){
      for( int i=0; i<hi[0]; ++i ){
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

  f1 <<= 2.0;
  f2 <<= 1.0;
  f2 <<= f2+f1;
  for( int k=0; k<hi[2]; ++k ){
    for( int j=0; j<hi[1]; ++j ){
      for( int i=0; i<hi[0]; ++i ){
        std::ostringstream msg;
        msg << "test_iterator 2: [" << i << "," << j << "," << k << "],  found: " << f2(i,j,k) << ", expected: 3.0";
        status( f2(i,j,k) == 3.0, msg.str() );
      }
    }
  }

  typename FieldT::iterator iter=f1.end()-1;
  for( size_t k=hi[2]; k>0; --k ){
    for( size_t j=hi[1]; j>0; --j ){
      for( size_t i=hi[0]; i>0; --i ){
        std::ostringstream msg;
        msg << "test_iterator backward iterator mem check f1: [" << i-1 << "," << j-1 << "," << k-1 << "]";
        status( &f1(i-1,j-1,k-1) == &*iter, msg.str() );
        --iter;
      }
    }
  }

  for( typename FieldT::iterator i=f1.begin(); i!=f1.end(); ++i ){
    *i = 0.0;
  }

  return status.ok();
}

//--------------------------------------------------------------------

template< typename FieldT >
bool test_interior( const IntVec npts,
                    const bool verbose )
{
  const MemoryWindow window( get_window_with_ghost<FieldT>(npts,true,true,true) );
  FieldT f1( window, NULL );
  FieldT f2( window, NULL );
  f1 <<= 2.0;
  const FieldT f3(f1);

  const MemoryWindow& interiorWindow = f1.window_without_ghost();
  const IntVec lo = interiorWindow.offset();
  const IntVec hi = lo + interiorWindow.extent();

  f2 <<= 0.0;
  // set interior values
  for( int k=lo[2]; k<hi[2]; ++k ){
    for( int j=lo[1]; j<hi[1]; ++j ){
      for( int i=lo[0]; i<hi[0]; ++i ){
        f2(i,j,k) = 1+i+j+k;
      }
    }
  }

  TestHelper status(verbose);

  typename FieldT::interior_iterator if2=f2.interior_begin();
  const typename FieldT::interior_iterator if2e=f2.interior_end();
  typename FieldT::const_interior_iterator if1=f1.interior_begin();
  const typename FieldT::interior_iterator if1e=f1.interior_end();
  typename FieldT::const_interior_iterator if3=f3.interior_begin();
  const typename FieldT::const_interior_iterator if3e=f3.interior_end();
  for( int k=lo[2]; k<hi[2]; ++k ){
    for( int j=lo[1]; j<hi[1]; ++j ){
      for( int i=lo[0]; i<hi[0]; ++i ){
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
        {
          const double& f3pt = f3(i,j,k);
          const double* f3pti = &*if3;
          std::ostringstream msg;  msg << "f3 mem loc at " << i<<","<<j<<","<<k;
          status( &f3pt == f3pti, msg.str() );
        }
        ++if2; ++if1; ++if3;
      }
    }
  }

  if1 = f1.interior_begin();
  if2 = f2.interior_begin();
  for( ; if2!=if2e; ++if1, ++if2 ){
    *if2 += *if1;
  }

  for( int k=lo[2]; k<hi[2]; ++k ){
    for( int j=lo[1]; j<hi[1]; ++j ){
      for( int i=lo[0]; i<hi[0]; ++i ){
        const double val = 1+i+j+k + 2.0;
        std::ostringstream msg;  msg << i<<","<<j<<","<<k << ",  found " << f2(i,j,k) << ", expected " << val;
        status( f2(i,j,k) == val, msg.str() );
      }
    }
  }
  return status.ok();
}

//--------------------------------------------------------------------

#ifdef ENABLE_THREADS
template<typename T>
struct ThreadWork{
  void doit(){
    const MemoryWindow ww = get_window_with_ghost<T>( IntVec(24,1,1), true,true,true );
    for( size_t i=0; i<100; ++i ){
      SpatFldPtr<T> f1 = SpatialFieldStore::get_from_window<T>( ww );
      SpatFldPtr<T> f2 = SpatialFieldStore::get_from_window<T>( ww );
      *f1 <<= 0.0;
    }
  }
};
#endif

//--------------------------------------------------------------------

template< typename FT1, typename FT2 >
bool test_store( const IntVec& dim, const IntVec& bc )
{
//  jcs_pause();
  TestHelper status(false);

  const MemoryWindow w1 = get_window_with_ghost<FT1>( dim, bc[0]==1, bc[1]==1, bc[2]==1 );
  const MemoryWindow w2 = get_window_with_ghost<FT2>( dim, bc[0]==1, bc[1]==1, bc[2]==1 );

# ifdef ENABLE_THREADS
  set_hard_thread_count( NTHREADS );
  ThreadWork<FT1> tw;
  for( int i=0; i<20; ++i ){
#   ifdef USE_FIFO
    ThreadPoolFIFO::self().schedule( boost::bind(&ThreadWork<FT1>::doit,tw) );
#   else
    ThreadPool::self().schedule( boost::threadpool::prio_task_func(1,boost::bind(&ThreadWork<FT1>::doit,tw)) );
#   endif
  }
  while( !ThreadPool::self().empty()
      || ThreadPool::self().active()>0
#     ifdef USE_FIFO
      || !ThreadPoolFIFO::self().empty()
      || ThreadPoolFIFO::self().active()>0
#     endif
      ){
    // force master thread to wait until queue is empty.
  }
# endif

  SpatFldPtr<FT1> f1 = SpatialFieldStore::get_from_window<FT1>( w1 );
  SpatFldPtr<FT2> f2 = SpatialFieldStore::get_from_window<FT2>( w2 );
  SpatFldPtr<FT1> f1a= SpatialFieldStore::get<FT1>( *f1 );
  SpatFldPtr<FT2> f2a= SpatialFieldStore::get<FT2>( *f2 );
  SpatFldPtr<FT2> f2b= SpatialFieldStore::get<FT2>( *f1 );

  status( f1->window_with_ghost() == f1a->window_with_ghost(), "f1==f1a" );
  status( f2->window_with_ghost() == f2a->window_with_ghost(), "f2==f2a" );
  status( f2->window_with_ghost() == f2b->window_with_ghost(), "f2==f2b" );

  FT1 f4( f1->window_with_ghost(), f1->field_values(), ExternalStorage );
  SpatFldPtr<FT1> f4a = SpatialFieldStore::get<FT1>(f4);
  status( f4a->window_with_ghost() == f1a->window_with_ghost() );

  return status.ok();
}

//--------------------------------------------------------------------

int main()
{
  TestHelper overall(true);

  bool verbose = false;

  boost::posix_time::ptime start, stop;
  start = boost::posix_time::microsec_clock::universal_time();
  try{
    overall( test_store<SVolField,  SVolField  >( IntVec(30,40,50), IntVec(0,0,0) ), "SVol,SVol(bc) store" );
    overall( test_store<SVolField,  SVolField  >( IntVec(30,40,50), IntVec(1,1,1) ), "SVol,SVol     store" );
    overall( test_store<SVolField,  SSurfXField>( IntVec(30,40,50), IntVec(0,0,0) ), "SVol,SSX      store" );
    overall( test_store<SSurfXField,SVolField  >( IntVec(30,40,50), IntVec(0,0,0) ), "SSX ,SVol     store" );
    overall( test_store<SVolField,  SSurfYField>( IntVec(30,40,50), IntVec(1,1,1) ), "SVol,SSY      store" );
    overall( test_store<SVolField,  SSurfYField>( IntVec(30,40,50), IntVec(0,0,0) ), "SVol,SSY      store" );
    overall( test_store<SVolField,  SSurfXField>( IntVec(30,40,50), IntVec(1,1,1) ), "SVol,SSX (bc) store" );
    overall( test_store<SSurfXField,SVolField  >( IntVec(30,40,50), IntVec(1,1,1) ), "SSX ,SVol(bc) store" );
    overall( test_store<SVolField,  SSurfZField>( IntVec(30,40,50), IntVec(1,1,1) ), "SVol,SSZ (bc) store" );
    overall( test_store<SVolField,  SSurfZField>( IntVec(30,40,50), IntVec(0,0,0) ), "SVol,SSZ      store" );
    overall( test_store<XSurfXField,SVolField  >( IntVec(30,40,50), IntVec(0,0,0) ), "XSX ,SVol     store" );
    overall( test_store<XSurfXField,ZVolField  >( IntVec(30,40,50), IntVec(0,0,0) ), "XSX ,ZVol     store" );
    overall( test_store<YSurfZField,ZSurfXField>( IntVec(30,40,50), IntVec(0,0,0) ), "XSX ,ZVol     store" );
  }
  catch(...){
    overall(false);
    std::cout << "exception thrown while running test_store" << std::endl;
    return -1;
  }
  stop = boost::posix_time::microsec_clock::universal_time();
  std::cout << "elapsed time (s): " << (stop-start).total_microseconds()*1e-6 << std::endl;

  overall( test_iterator<SVolField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) SVolField" );
  overall( test_iterator<SVolField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) SVolField" );
  overall( test_iterator<SVolField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) SVolField" );

  std::cout << std::endl;

  overall( test_iterator<SSurfXField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) SSurfXField" );
  overall( test_iterator<SSurfXField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) SSurfXField" );
  overall( test_iterator<SSurfXField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) SSurfXField" );

  std::cout << std::endl;

  overall( test_iterator<SSurfYField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) SSurfYField" );
  overall( test_iterator<SSurfYField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) SSurfYField" );
  overall( test_iterator<SSurfYField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) SSurfYField" );

  std::cout << std::endl;

  overall( test_iterator<SSurfZField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) SSurfZField" );
  overall( test_iterator<SSurfZField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) SSurfZField" );
  overall( test_iterator<SSurfZField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) SSurfZField" );

  std::cout << std::endl << std::endl;

  overall( test_iterator<XVolField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) XVolField" );
  overall( test_iterator<XVolField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) XVolField" );
  overall( test_iterator<XVolField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) XVolField" );

  std::cout << std::endl;

  overall( test_iterator<XSurfXField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) XSurfXField" );
  overall( test_iterator<XSurfXField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) XSurfXField" );
  overall( test_iterator<XSurfXField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) XSurfXField" );

  std::cout << std::endl;

  overall( test_iterator<XSurfYField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) XSurfYField" );
  overall( test_iterator<XSurfYField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) XSurfYField" );
  overall( test_iterator<XSurfYField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) XSurfYField" );

  std::cout << std::endl;

  overall( test_iterator<XSurfZField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) XSurfZField" );
  overall( test_iterator<XSurfZField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) XSurfZField" );
  overall( test_iterator<XSurfZField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) XSurfZField" );

  std::cout << std::endl << std::endl;

  overall( test_iterator<YVolField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) YVolField" );
  overall( test_iterator<YVolField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) YVolField" );
  overall( test_iterator<YVolField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) YVolField" );

  std::cout << std::endl;

  overall( test_iterator<YSurfXField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) YSurfXField" );
  overall( test_iterator<YSurfXField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) YSurfXField" );
  overall( test_iterator<YSurfXField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) YSurfXField" );

  std::cout << std::endl;

  overall( test_iterator<YSurfYField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) YSurfYField" );
  overall( test_iterator<YSurfYField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) YSurfYField" );
  overall( test_iterator<YSurfYField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) YSurfYField" );

  std::cout << std::endl;

  overall( test_iterator<YSurfZField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) YSurfZField" );
  overall( test_iterator<YSurfZField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) YSurfZField" );
  overall( test_iterator<YSurfZField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) YSurfZField" );

  std::cout << std::endl << std::endl;

  overall( test_iterator<ZVolField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) ZVolField" );
  overall( test_iterator<ZVolField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) ZVolField" );
  overall( test_iterator<ZVolField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) ZVolField" );

  std::cout << std::endl;

  overall( test_iterator<ZSurfXField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) ZSurfXField" );
  overall( test_iterator<ZSurfXField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) ZSurfXField" );
  overall( test_iterator<ZSurfXField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) ZSurfXField" );

  std::cout << std::endl;

  overall( test_iterator<ZSurfYField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) ZSurfYField" );
  overall( test_iterator<ZSurfYField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) ZSurfYField" );
  overall( test_iterator<ZSurfYField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) ZSurfYField" );

  std::cout << std::endl;

  overall( test_iterator<ZSurfZField>( IntVec(3,3,3), verbose ), "test_iterator (3,3,3) ZSurfZField" );
  overall( test_iterator<ZSurfZField>( IntVec(3,4,1), verbose ), "test_iterator (3,4,1) ZSurfZField" );
  overall( test_iterator<ZSurfZField>( IntVec(4,3,2), verbose ), "test_iterator (4,3,2) ZSurfZField" );

  std::cout << std::endl << std::endl;

  // -------------------------------------------------------------- //
  // --------------------  Test interior Iterators ---------------- //
  // -------------------------------------------------------------- //
  verbose=false;
  overall( test_interior<SVolField  >( IntVec(3,3,3), verbose ), "test_interior (3,3,3) SVolField" );
  overall( test_interior<SSurfXField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) SSurfXField" );
  overall( test_interior<SSurfYField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) SSurfYField" );
  overall( test_interior<SSurfZField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) SSurfZField" );

  overall( test_interior<XVolField  >( IntVec(3,3,3), verbose ), "test_interior (3,3,3) XVolField" );
  overall( test_interior<XSurfXField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) XSurfXField" );
  overall( test_interior<XSurfYField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) XSurfYField" );
  overall( test_interior<XSurfZField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) XSurfZField" );

  overall( test_interior<YVolField  >( IntVec(3,3,3), verbose ), "test_interior (3,3,3) YVolField" );
  overall( test_interior<YSurfXField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) YSurfXField" );
  overall( test_interior<YSurfYField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) YSurfYField" );
  overall( test_interior<YSurfZField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) YSurfZField" );

  overall( test_interior<ZVolField  >( IntVec(3,3,3), verbose ), "test_interior (3,3,3) ZVolField" );
  overall( test_interior<ZSurfXField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) ZSurfXField" );
  overall( test_interior<ZSurfYField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) ZSurfYField" );
  overall( test_interior<ZSurfZField>( IntVec(3,3,3), verbose ), "test_interior (3,3,3) ZSurfZField" );

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

    svol2 <<= 2.0;
    svol1 <<= (svol1 + svol2) * svol2 / svol2 - svol2;

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
      SpatFldPtr<SVolField> sv3 = SpatialFieldStore::get<SVolField>( svol1 );
      *sv3 = svol1;
      status( sv3->field_equal(svol1, 0.0), "spatial field pointer from store" );
    }

    overall( status.ok(), "field operations" );
  }

  //test field_equal
  overall(test_field_equal<SVolField>(IntVec(10, 11, 12), LOCAL_RAM, LOCAL_RAM, false), "LOCAL_RAM x LOCAL_RAM Equal Test");
#ifdef __CUDACC__
  overall(test_field_equal<SVolField>(IntVec(10, 11, 12), LOCAL_RAM, EXTERNAL_CUDA_GPU, false), "LOCAL_RAM x EXTERNAL_CUDA_GPU Equal Test");
  overall(test_field_equal<SVolField>(IntVec(10, 11, 12), EXTERNAL_CUDA_GPU, LOCAL_RAM, false), "EXTERNAL_CUDA_GPU x LOCAL_RAM Equal Test");
  overall(test_field_equal<SVolField>(IntVec(10, 11, 12), EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, false), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Equal Test");
#endif

  if( overall.isfailed() ){
    std::cout << "FAIL!" << std::endl;
    return -1;
  }
  std::cout << "PASS" << std::endl;

  return 0;
}
