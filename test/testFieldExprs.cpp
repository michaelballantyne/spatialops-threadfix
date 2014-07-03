#include <spatialops/Nebo.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FieldComparisons.h>
#include "TestHelper.h"

#include <numeric>

namespace SS = SpatialOps::structured;

using namespace SpatialOps;

template< typename FieldT >
bool test( const SS::IntVec dim )
{
  TestHelper status(false);
  const SS::GhostData ghost(1);
  const SS::BoundaryCellInfo bc = SS::BoundaryCellInfo::build<FieldT>(true,true,true);
  const SS::MemoryWindow w( SS::get_window_with_ghost(dim,ghost,bc) );
  FieldT f1(w,bc,ghost,NULL),
         f2(w,bc,ghost,NULL),
         f3(w,bc,ghost,NULL),
         f4(w,bc,ghost,NULL),
          x(w,bc,ghost,NULL);

  const SS::IntVec& globDim = w.glob_dim();
  const double dx = 3.1415 / (globDim[0]);
  for( int k=0; k<globDim[2]; ++k ){
    for( int j=0; j<globDim[1]; ++j ){
      for( int i=0; i<globDim[0]; ++i ){
        const size_t ix = w.flat_index( SS::IntVec(i,j,k) );
        x[ ix ] = dx*i;
      }
    }
  }

  typedef typename FieldT::const_iterator constiter;
  typedef typename FieldT::iterator       iter;

  {
    f1 <<= sin( x ) + 3.0;
    constiter i1=f1.begin(), ix=x.begin();
    for( iter i2=f2.begin(); i2!=f2.end(); ++i2, ++ix, ++i1 ){
      *i2 = sin( *ix ) + 3.0;
    }
    status( field_equal_ulp(f2, f1, 1), "sin(x)+3" );
  }

  {
    constiter ix=x.begin();
    iter i2=f2.begin();
    for( iter i1=f1.begin(); i1!=f1.end(); ++i1, ++i2, ++ix ){
      *i1 = std::sin( *ix ) + 4.0;
      *i2 = std::cos( *ix ) + 3.0;
    }
  }
  f3 <<= cos(x) + 3.0;
  status( field_equal_ulp(f2, f3, 1), "cos(x)" );

  f3 <<= f1+(f2*f1)-f2/f1;

  double l2norm = 0;
  constiter i1=f1.begin(), i2=f2.begin(), i3=f3.begin();
  for( iter i4=f4.begin(); i4!=f4.end(); ++i4, ++i3, ++i2, ++i1 ){
    *i4 = *i1 + (*i2 * *i1) - *i2 / *i1;
    l2norm += (*i4 * *i4);
  }
  l2norm = sqrt(l2norm);
  status( field_equal_ulp(f4, f3, 1), "a+(a*b)-b/a" );

  status( l2norm == field_norm(f4), "norm" );
  status( *std::max_element(f4.begin(),f4.end()) == field_max(f4), "max" );
  status( *std::min_element(f4.begin(),f4.end()) == field_min(f4), "min" );
  status( std::accumulate(f4.begin(),f4.end(),0.0) == field_sum(f4), "sum" );

  return status.ok();
}


bool drive_test( const SS::IntVec& dim )
{
  TestHelper status( true );

  status( test<SS::SVolField  >( dim ), "SVolField" );
  status( test<SS::SSurfXField>( dim ), "SSurfXField" );
  status( test<SS::SSurfYField>( dim ), "SSurfYField" );
  status( test<SS::SSurfZField>( dim ), "SSurfZField" );

  status( test<SS::XVolField  >( dim ), "XVolField" );
  status( test<SS::XSurfXField>( dim ), "XSurfXField" );
  status( test<SS::XSurfYField>( dim ), "XSurfYField" );
  status( test<SS::XSurfZField>( dim ), "XSurfZField" );

  status( test<SS::YVolField  >( dim ), "YVolField" );
  status( test<SS::YSurfXField>( dim ), "YSurfXField" );
  status( test<SS::YSurfYField>( dim ), "YSurfYField" );
  status( test<SS::YSurfZField>( dim ), "YSurfZField" );

  status( test<SS::ZVolField  >( dim ), "ZVolField" );
  status( test<SS::ZSurfXField>( dim ), "ZSurfXField" );
  status( test<SS::ZSurfYField>( dim ), "ZSurfYField" );
  status( test<SS::ZSurfZField>( dim ), "ZSurfZField" );

  return status.ok();
}


int main()
{
  TestHelper status( true );

  status( drive_test( SS::IntVec(10,1,1) ), "Dimension: (10,1,1) tests");
  status( drive_test( SS::IntVec(1,10,1) ), "Dimension: (1,10,1) tests");
  status( drive_test( SS::IntVec(1,1,10) ), "Dimension: (1,1,10) tests");

  status( drive_test( SS::IntVec(10,10,1) ), "Dimension: (10,10,1) tests");
  status( drive_test( SS::IntVec(10,1,10) ), "Dimension: (10,1,10) tests");
  status( drive_test( SS::IntVec(1,10,10) ), "Dimension: (1,10,10) tests");

  status( drive_test( SS::IntVec(10,10,10) ), "Dimension: (10,10,10) tests");

  if( status.ok() ) return 0;
  return -1;
}
