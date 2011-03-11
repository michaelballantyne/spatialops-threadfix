#ifndef SpatialOps_1DOpTypes_h
#define SpatialOps_1DOpTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/SpatialOperator.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/FVOneDimensional.h>


namespace SpatialOps{

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
   *  typedef OperatorTypeBuilder<Divergence,SSurfXField,SVolField>::type DivX;
   *  typedef OperatorTypeBuilder<Gradient,VolT,FaceTypes<VolT>::XFace>::type GradX;
   *  \endcode
   *
   *  Note that we only provide fully specialized versions of this template
   *  so that unsupported operator types cannot be inadvertantly formed.
   */
  template< typename OpT, typename SrcT, typename DestT >
  struct OperatorTypeBuilder;

#define OP_BUILDER( OP, SRC, DEST )                             \
  template<>                                                    \
  struct OperatorTypeBuilder<OP,SRC,DEST>{                      \
    typedef SpatialOperator< LinAlg, OP, SRC, DEST > type;      \
  };

  OP_BUILDER( Interpolant, structured::SVolField,   structured::SSurfXField )
  OP_BUILDER( Interpolant, structured::SSurfXField, structured::SVolField   )
  OP_BUILDER( Gradient,    structured::SVolField,   structured::SSurfXField )
  OP_BUILDER( Divergence,  structured::SSurfXField, structured::SVolField   )

} // namespace SpatialOps

#endif // SpatialOps_1DOpTypes_h
