#include <spatialops/pointfield/PointFieldTypes.h>
#include <spatialops/pointfield/PointOperators.h>

#include <spatialops/structured/FVStaggeredBCTools.h>
#include <spatialops/structured/GhostData.h>

#include <test/TestHelper.h>

#include <iostream>

typedef SpatialOps::Point::PointField FieldT;
namespace SS=SpatialOps::structured;

int main()
{
  const size_t npts = 10;
  const SS::GhostDataRT ghost(0);
  const SS::BoundaryCellInfo bc = SS::BoundaryCellInfo::build<FieldT>();
  const SS::MemoryWindow mw( SS::IntVec( npts, 1, 1 ) );

  FieldT f( mw, bc, ghost, NULL );

  double x=0.1;
  for( FieldT::iterator ifld=f.begin(); ifld!=f.end(); ++ifld, x+=1.0 ){
    *ifld = x;
  }

  TestHelper status(true);

  FieldT::interior_iterator i2=f.interior_begin();
  for( FieldT::iterator i=f.begin(); i!=f.end(); ++i, ++i2 ){
    status( *i==*i2, "value" );
    status( &*i == &*i2, "address" );
  }

  {
    typedef SS::ConstValEval BCVal;
    typedef SS::BoundaryCondition<FieldT,BCVal> BC;

    BC bc1( SS::IntVec(2,1,1), BCVal(1.234) );
    BC bc2( SS::IntVec(4,1,1), BCVal(3.456) );
  
    bc1(f);
    bc2(f);
    status( f[2] == 1.234, "point BC 1" );
    status( f[4] == 3.456, "point BC 2" );
  }

  {
    std::vector<size_t> ix;
    ix.push_back(4);
    ix.push_back(2);
    SpatialOps::Point::FieldToPoint<FieldT> ftp(ix);
    SpatialOps::Point::PointToField<FieldT> ptf(ix);

    const SS::MemoryWindow mw2( SpatialOps::structured::IntVec(2,1,1) );
    FieldT f2( mw2, bc, ghost, NULL );
    ftp.apply_to_field( f, f2 );
    status( f2[0] == 3.456, "Field2Point Interp (1)" );
    status( f2[1] == 1.234, "Field2Point Interp (2)" );

    f2[0] = -1.234;
    f2[1] = -3.456;
    ptf.apply_to_field( f2, f );
    status( f[2] == -3.456, "Point2Field Interp (1)" );
    status( f[4] == -1.234, "Point2Field Interp (2)" );
  }


  if( status.ok() ) return 0;
  return -1;
}
