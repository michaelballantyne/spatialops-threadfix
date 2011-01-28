#include "StencilBuilder.h"
#include "FVStaggeredOperatorTypes.h"

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/OperatorDatabase.h>

namespace SpatialOps{
namespace structured{

#define REG_BASIC_OP_TYPES( VOL )                                       \
  {                                                                     \
    typedef BasicOpTypes<VOL>  OpTypes;                                 \
    opdb.register_new_operator( new OpTypes::InterpC2FX( 0.5, 0.5 ) );  \
    opdb.register_new_operator( new OpTypes::InterpC2FY( 0.5, 0.5 ) );  \
    opdb.register_new_operator( new OpTypes::InterpC2FZ( 0.5, 0.5 ) );  \
    opdb.register_new_operator( new OpTypes::GradX( -dx, dx ) );        \
    opdb.register_new_operator( new OpTypes::GradY( -dy, dy ) );        \
    opdb.register_new_operator( new OpTypes::GradZ( -dz, dz ) );        \
    opdb.register_new_operator( new OpTypes::DivX ( -dx, dx ) );        \
    opdb.register_new_operator( new OpTypes::DivY ( -dy, dy ) );        \
    opdb.register_new_operator( new OpTypes::DivZ ( -dz, dz ) );        \
  }
  
  //------------------------------------------------------------------

  void build_stencils( const unsigned int nx,
                       const unsigned int ny,
                       const unsigned int nz,
                       const double Lx,
                       const double Ly,
                       const double Lz,
                       OperatorDatabase& opdb )
  {
    const double dx = Lx/nx;
    const double dy = Ly/ny;
    const double dz = Lz/nz;

    REG_BASIC_OP_TYPES( SVolField )  // basic operator types on a scalar volume
    REG_BASIC_OP_TYPES( XVolField )  // basic operator types on a x volume
    REG_BASIC_OP_TYPES( YVolField )  // basic operator types on a y volume
    REG_BASIC_OP_TYPES( ZVolField )  // basic operator types on a z volume

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,XVolField,YSurfXField>::type( 0.5, 0.5 ) );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,XVolField,ZSurfXField>::type( 0.5, 0.5 ) );

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,YVolField,XSurfYField>::type( 0.5, 0.5 ) );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,YVolField,ZSurfYField>::type( 0.5, 0.5 ) );

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,ZVolField,XSurfZField>::type( 0.5, 0.5 ) );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,ZVolField,YSurfZField>::type( 0.5, 0.5 ) );

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,XSurfXField>::type() );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,YSurfYField>::type() );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,ZSurfZField>::type() );

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,XVolField,SSurfXField>::type() );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,YVolField,SSurfYField>::type() );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,ZVolField,SSurfZField>::type() );


    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,XSurfYField>::type(.25,.25,.25,.25) );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,XSurfZField>::type(.25,.25,.25,.25) );

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,YSurfXField>::type(.25,.25,.25,.25) );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,YSurfZField>::type(.25,.25,.25,.25) );

    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,ZSurfXField>::type(.25,.25,.25,.25) );
    opdb.register_new_operator( new OperatorTypeBuilder<Interpolant,SVolField,ZSurfYField>::type(.25,.25,.25,.25) );

  }

  //------------------------------------------------------------------

} // namespace structured
} // namespace SpatialOps
