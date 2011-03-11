#ifndef SpatialOps_structured_FVStaggeredOpTypes_h
#define SpatialOps_structured_FVStaggeredOpTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/stencil/Stencil2.h>
#include <spatialops/structured/stencil/NullStencil.h>
#include <spatialops/structured/stencil/Stencil4.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \struct OperatorTypeBuilder
   *  \author James C. Sutherland
   *
   *  \brief Builds operator types
   *  \tparam OpT the type of operator (\c Interpolant, \c Gradient, \c Divergence)
   *  \tparam SrcT the field type that the operator acts on
   *  \tparam DestT the field type that the operator produces
   *
   *  Implementations of this struct define a public type called \c
   *  type.  This is the fully qualified operator type.
   *
   *  \par Troubleshooting
   *  Note that if the compiler fails, it is likely because the
   *  requested operator type is not supported.  There is no default
   *  implementation for this struct.  All implementations are fully
   *  specialized for the supported types.
   *
   *  \par Example Usage
   *  \code
   *  typedef OperatorTypeBuilder<Interpolant,SVolField,XVolField>::type InterpSVolXVol;
   *  typedef OperatorTypeBuilder<Divergence,XSurfYField,XVolField>::type DivX;
   *  typedef OperatorTypeBuilder<Gradient,VolT,FaceTypes<VolT>::XFace>::type GradX;
   *  \endcode
   *
   *  Note that we only provide fully specialized versions of this template
   *  so that unsupported operator types cannot be inadvertantly formed.
   */
  template<typename OpT, typename SrcT, typename DestT>
  struct OperatorTypeBuilder;


#define OP_BUILDER( STENCIL, OP, SRC, DEST )    \
  template<>                                    \
  struct OperatorTypeBuilder<OP,SRC,DEST>{      \
    typedef STENCIL<OP,SRC,DEST> type;          \
  };


  /**
   *  \struct BasicOpTypes
   *  \author James C. Sutherland
   *
   *  \brief Provides typedefs for common operator types on a given volume
   *  \tparam CellT the type of volume we are considering.
   *
   *  The following public typedefs are made:
   *   - \c InterpC2FX Interpolate to the x-surface
   *   - \c InterpC2FY Interpolate to the y-surface
   *   - \c InterpC2FZ Interpolate to the z-surface
   *   - \c GradX Calculate \f$\frac{\partial}{\partial x}$\f on the x-surface
   *   - \c GradY Calculate \f$\frac{\partial}{\partial y}$\f on the y-surface
   *   - \c GradZ Calculate \f$\frac{\partial}{\partial z}$\f on the z-surface
   *   - \c DivX Calculate the divergence from the x-surface to the volume
   *   - \c DivY Calculate the divergence from the y-surface to the volume
   *   - \c DivZ Calculate the divergence from the z-surface to the volume
   *
   *  This struct is defined for the following field types:
   *   - \c SVolField
   *   - \c XVolField
   *   - \c YVolField
   *   - \c ZVolField
   *
   *  Examples:
   *  \code
   *  typedef BasicOpTypes<SVolField>::GradY      MyGradYType;
   *  typedef BasicOpTypes<SVolField>::InterpC2FX InterpSVolSSurfX;
   *  \endcode
   */
  template< typename CellT > struct BasicOpTypes;


#define BASIC_OPTYPE_BUILDER( VOL )                                     \
  OP_BUILDER( Stencil2, Interpolant, VOL, FaceTypes<VOL>::XFace )       \
  OP_BUILDER( Stencil2, Interpolant, VOL, FaceTypes<VOL>::YFace )       \
  OP_BUILDER( Stencil2, Interpolant, VOL, FaceTypes<VOL>::ZFace )       \
  OP_BUILDER( Stencil2, Gradient,    VOL, FaceTypes<VOL>::XFace )       \
  OP_BUILDER( Stencil2, Gradient,    VOL, FaceTypes<VOL>::YFace )       \
  OP_BUILDER( Stencil2, Gradient,    VOL, FaceTypes<VOL>::ZFace )       \
  OP_BUILDER( Stencil2, Divergence,  FaceTypes<VOL>::XFace, VOL )       \
  OP_BUILDER( Stencil2, Divergence,  FaceTypes<VOL>::YFace, VOL )       \
  OP_BUILDER( Stencil2, Divergence,  FaceTypes<VOL>::ZFace, VOL )       \
  template<>                                                            \
  struct BasicOpTypes<VOL>                                              \
  {                                                                     \
    typedef OperatorTypeBuilder< Interpolant, VOL, FaceTypes<VOL>::XFace >::type InterpC2FX; \
    typedef OperatorTypeBuilder< Interpolant, VOL, FaceTypes<VOL>::YFace >::type InterpC2FY; \
    typedef OperatorTypeBuilder< Interpolant, VOL, FaceTypes<VOL>::ZFace >::type InterpC2FZ; \
    typedef OperatorTypeBuilder< Gradient,    VOL, FaceTypes<VOL>::XFace >::type GradX; \
    typedef OperatorTypeBuilder< Gradient,    VOL, FaceTypes<VOL>::YFace >::type GradY; \
    typedef OperatorTypeBuilder< Gradient,    VOL, FaceTypes<VOL>::ZFace >::type GradZ; \
    typedef OperatorTypeBuilder< Divergence,  FaceTypes<VOL>::XFace, VOL >::type DivX; \
    typedef OperatorTypeBuilder< Divergence,  FaceTypes<VOL>::YFace, VOL >::type DivY; \
    typedef OperatorTypeBuilder< Divergence,  FaceTypes<VOL>::ZFace, VOL >::type DivZ; \
  };


  BASIC_OPTYPE_BUILDER( SVolField )
  BASIC_OPTYPE_BUILDER( XVolField )
  BASIC_OPTYPE_BUILDER( YVolField )
  BASIC_OPTYPE_BUILDER( ZVolField )


  OP_BUILDER( Stencil2, Interpolant, XVolField, YSurfXField )
  OP_BUILDER( Stencil2, Gradient,    XVolField, YSurfXField )
  OP_BUILDER( Stencil2, Interpolant, XVolField, ZSurfXField )
  OP_BUILDER( Stencil2, Gradient,    XVolField, ZSurfXField )

  OP_BUILDER( Stencil2, Interpolant, YVolField, XSurfYField )
  OP_BUILDER( Stencil2, Gradient,    YVolField, XSurfYField )
  OP_BUILDER( Stencil2, Interpolant, YVolField, ZSurfYField )
  OP_BUILDER( Stencil2, Gradient,    YVolField, ZSurfYField )

  OP_BUILDER( Stencil2, Interpolant, ZVolField, XSurfZField )
  OP_BUILDER( Stencil2, Gradient,    ZVolField, XSurfZField )
  OP_BUILDER( Stencil2, Interpolant, ZVolField, YSurfZField )
  OP_BUILDER( Stencil2, Gradient,    ZVolField, YSurfZField )

  OP_BUILDER( Stencil2, Interpolant, SVolField, XVolField )
  OP_BUILDER( Stencil2, Gradient,    SVolField, XVolField )

  OP_BUILDER( Stencil2, Interpolant, SVolField, YVolField )
  OP_BUILDER( Stencil2, Gradient,    SVolField, YVolField )

  OP_BUILDER( Stencil2, Interpolant, SVolField, ZVolField )
  OP_BUILDER( Stencil2, Gradient,    SVolField, ZVolField )

  OP_BUILDER( Stencil2, Interpolant, XVolField, SVolField )
  OP_BUILDER( Stencil2, Gradient,    XVolField, SVolField )

  OP_BUILDER( Stencil2, Interpolant, YVolField, SVolField )
  OP_BUILDER( Stencil2, Gradient,    YVolField, SVolField )

  OP_BUILDER( Stencil2, Interpolant, ZVolField, SVolField )
  OP_BUILDER( Stencil2, Gradient,    ZVolField, SVolField )


  OP_BUILDER( NullStencil, Interpolant, XVolField, SSurfXField )
  OP_BUILDER( NullStencil, Interpolant, YVolField, SSurfYField )
  OP_BUILDER( NullStencil, Interpolant, ZVolField, SSurfZField )

  OP_BUILDER( NullStencil, Interpolant, SVolField, XSurfXField )
  OP_BUILDER( NullStencil, Interpolant, SVolField, YSurfYField )
  OP_BUILDER( NullStencil, Interpolant, SVolField, ZSurfZField )


  OP_BUILDER( Stencil4, Interpolant, SVolField, XSurfYField )
  OP_BUILDER( Stencil4, Interpolant, SVolField, XSurfZField )

  OP_BUILDER( Stencil4, Interpolant, SVolField, YSurfXField )
  OP_BUILDER( Stencil4, Interpolant, SVolField, YSurfZField )

  OP_BUILDER( Stencil4, Interpolant, SVolField, ZSurfXField )
  OP_BUILDER( Stencil4, Interpolant, SVolField, ZSurfYField )

} // namespace SpatialOps
} // namespace structured

#endif // SpatialOps_structured_FVStaggeredOpTypes_h
