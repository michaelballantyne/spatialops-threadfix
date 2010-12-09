#include <spatialops/FieldOperationDefinitions.h>
#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>

#include "TestHelper.h"

namespace SS = SpatialOps::structured;


template< typename FieldT >
bool test( const SS::IntVec dim )
{
  TestHelper status(true);
  const SS::MemoryWindow w( SS::get_window_with_ghost<FieldT>(dim,true,true,true) );
  FieldT f1(w,NULL), f2(w,NULL), f3(w,NULL), f4(w,NULL);
  FieldT x(w,NULL);

  const SS::IntVec& globDim = w.glob_dim();
  const double dx = 3.1415 / (globDim[0]-1);
  for( size_t k=0; k<globDim[2]; ++k ){
    for( size_t j=0; j<globDim[1]; ++j ){
      for( size_t i=0; i<globDim[0]; ++i ){
        const size_t ix = w.flat_index( SS::IntVec(i,j,k) );
        x[ ix ] = dx*i;
      }
    }
  }

  typedef typename FieldT::const_iterator constiter;
  typedef typename FieldT::iterator       iter;

  {
    TestHelper tmp(false);
    f1 <<= sin( x );
    constiter i1=f1.begin(), ix=x.begin();
    for( iter i2=f2.begin(); i2!=f2.end(); ++i2, ++ix, ++i1 ){
      *i2 = sin( *ix );
      tmp( *i1 == *i2 );
    }
    status( tmp.ok(), "sin(x)" );
  }

  {
    TestHelper tmp(false);
    {
      constiter ix=x.begin();
      iter i2=f2.begin();
      for( iter i1=f1.begin(); i1!=f1.end(); ++i1, ++i2, ++ix ){
        *i1 = std::sin( *ix );
        *i2 = std::cos( *ix );
      }
    }

    f3 <<= f1+(f2*f1)-f2;

    constiter i1=f1.begin(), i2=f2.begin(), i3=f3.begin();
    for( iter i4=f4.begin(); i4!=f4.end(); ++i4, ++i3, ++i2, ++i1 ){
      *i4 = *i1 + (*i2 * *i1) - *i2;
      tmp( *i4 == *i3 );
    }
    status( tmp.ok(), "a+(a*b)-b" );
  }

  return true;
}


int main()
{
  TestHelper status( true );
  status( test<SS::SVolField>( SS::IntVec(10,1,1) ), "SVolField" );

  if( status.ok() ) return 0;
  return -1;
}
