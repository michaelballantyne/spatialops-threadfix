#ifndef SpatialOps_FVStaggeredOps_h
#define SpatialOps_FVStaggeredOps_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/SpatialOperator.h>
#include "FVStaggeredFieldTypes.h"

namespace SpatialOps{
namespace structured{

  /**
   *  \file FVStaggeredOperatorTypes.h
   *
   *  \addtogroup structured
   *  @{
   *  \addtogroup operators
   *  @{
   */

  /**
   *  \struct OperatorTypeBuilder
   *  \brief Convenience definition for a SpatialOperator
   *
   *  Supplies a typedef defining \c type, which defines the operator.
   */
  template<typename Op, typename SrcT, typename DestT>
  struct OperatorTypeBuilder{
    typedef SpatialOps::SpatialOperator< LinAlg, Op, SrcT, DestT >  type;
  };


  /**
   *  \struct BasicOpTypes
   *  \brief provides typedefs for various operators related to the given cell type
   *
   *  \par Defined Types
   *
   *    - \b GradX - x-gradient of cell-centered quantities produced at cell x-faces
   *    - \b GradY - y-gradient of cell-centered quantities produced at cell y-faces
   *    - \b GradZ - z-gradient of cell-centered quantities produced at cell z-faces
   *
   *    - \b DivX - x-divergence of cell-centered quantities produced at cell x-faces
   *    - \b DivY - y-divergence of cell-centered quantities produced at cell y-faces
   *    - \b DivZ - z-divergence of cell-centered quantities produced at cell z-faces
   *
   *    - \b InterpC2FX - Interpolate cell-centered quantities to x-faces
   *    - \b InterpC2FY - Interpolate cell-centered quantities to y-faces
   *    - \b InterpC2FZ - Interpolate cell-centered quantities to z-faces
   *
   *    - \b InterpF2CX - Interpolate x-face quantities to cell-centered
   *    - \b InterpF2CY - Interpolate y-face quantities to cell-centered
   *    - \b InterpF2CZ - Interpolate z-face quantities to cell-centered
   *  
   *  Example:
   *  \code
   *    typedef BasicOpTypes< SVolField >  VolOps;
   *    typedef VolOps::GradX         GradX;
   *    typedef VolOps::DivX          DivX;
   *  \endcode
   */
  template< typename CellT > struct BasicOpTypes
  {
    typedef typename OperatorTypeBuilder< Gradient,
                                          CellT,
                                          typename FaceTypes<CellT>::XFace >::type	GradX;
    typedef typename OperatorTypeBuilder< Gradient,
                                          CellT,
                                          typename FaceTypes<CellT>::YFace >::type	GradY;
    typedef typename OperatorTypeBuilder< Gradient,
                                          CellT,
                                          typename FaceTypes<CellT>::ZFace >::type	GradZ;

    typedef typename OperatorTypeBuilder< Divergence,
                                          typename FaceTypes<CellT>::XFace,
                                          CellT >::type					DivX;
    typedef typename OperatorTypeBuilder< Divergence,
                                          typename FaceTypes<CellT>::YFace,
                                          CellT >::type 				DivY;
    typedef typename OperatorTypeBuilder< Divergence,
                                          typename FaceTypes<CellT>::ZFace,
                                          CellT >::type 				DivZ;

    typedef typename OperatorTypeBuilder< Interpolant,
                                          CellT,
                                          typename FaceTypes<CellT>::XFace >::type 	InterpC2FX;
    typedef typename OperatorTypeBuilder< Interpolant,
                                          CellT,
                                          typename FaceTypes<CellT>::YFace >::type 	InterpC2FY;
    typedef typename OperatorTypeBuilder< Interpolant,
                                          CellT,
                                          typename FaceTypes<CellT>::ZFace >::type 	InterpC2FZ;

    typedef typename OperatorTypeBuilder< Interpolant,
                                          typename FaceTypes<CellT>::XFace,
                                          CellT >::type 				InterpF2CX;
    typedef typename OperatorTypeBuilder< Interpolant,
                                          typename FaceTypes<CellT>::YFace,
                                          CellT >::type 				InterpF2CY;
    typedef typename OperatorTypeBuilder< Interpolant,
                                          typename FaceTypes<CellT>::ZFace,
                                          CellT >::type 				InterpF2CZ;
  };

  // #################################################################
  //
  // NOTE: the typedefs below are for convenience.  Using
  //       OperatorTypeBuilder and BasicOpTypes, one can easily construct
  //       operators that convert between field types as desired.
  //
  // #################################################################


  //-- Interpolant Operators --//

  // scalar volume to scalar surfaces (diffusion coefficients)
  typedef BasicOpTypes<SVolField>::InterpC2FX  InterpSVolSSurfX;
  typedef BasicOpTypes<SVolField>::InterpC2FY  InterpSVolSSurfY;
  typedef BasicOpTypes<SVolField>::InterpC2FZ  InterpSVolSSurfZ;

  // scalar volume to staggered surfaces (viscosity, dilatation)
  typedef OperatorTypeBuilder< Interpolant, SVolField, XSurfXField >::type  InterpSVolXSurfX;
  typedef OperatorTypeBuilder< Interpolant, SVolField, XSurfYField >::type  InterpSVolXSurfY;
  typedef OperatorTypeBuilder< Interpolant, SVolField, XSurfZField >::type  InterpSVolXSurfZ;

  typedef OperatorTypeBuilder< Interpolant, SVolField, YSurfXField >::type  InterpSVolYSurfX;
  typedef OperatorTypeBuilder< Interpolant, SVolField, YSurfYField >::type  InterpSVolYSurfY;
  typedef OperatorTypeBuilder< Interpolant, SVolField, YSurfZField >::type  InterpSVolYSurfZ;

  typedef OperatorTypeBuilder< Interpolant, SVolField, ZSurfXField >::type  InterpSVolZSurfX;
  typedef OperatorTypeBuilder< Interpolant, SVolField, ZSurfYField >::type  InterpSVolZSurfY;
  typedef OperatorTypeBuilder< Interpolant, SVolField, ZSurfZField >::type  InterpSVolZSurfZ;


  // scalar volume to staggered volumes (density)
  typedef OperatorTypeBuilder< Interpolant, SVolField, XVolField   >::type  InterpSVolXVol;
  typedef OperatorTypeBuilder< Interpolant, SVolField, YVolField   >::type  InterpSVolYVol;
  typedef OperatorTypeBuilder< Interpolant, SVolField, ZVolField   >::type  InterpSVolZVol;


  // staggered volumes to scalar surface (advecting velocities, etc.)
  typedef OperatorTypeBuilder< Interpolant, XVolField, SSurfXField >::type  InterpXVolSSurfX;
  typedef OperatorTypeBuilder< Interpolant, YVolField, SSurfYField >::type  InterpYVolSSurfY;
  typedef OperatorTypeBuilder< Interpolant, ZVolField, SSurfZField >::type  InterpZVolSSurfZ;


  // staggered volumes to staggered surfaces (momentum solution components)
  typedef BasicOpTypes<XVolField>::InterpC2FX  InterpXVolXSurfX;
  typedef BasicOpTypes<XVolField>::InterpC2FY  InterpXVolXSurfY;
  typedef BasicOpTypes<XVolField>::InterpC2FZ  InterpXVolXSurfZ;

  typedef BasicOpTypes<YVolField>::InterpC2FX  InterpYVolYSurfX;
  typedef BasicOpTypes<YVolField>::InterpC2FY  InterpYVolYSurfY;
  typedef BasicOpTypes<YVolField>::InterpC2FZ  InterpYVolYSurfZ;

  typedef BasicOpTypes<ZVolField>::InterpC2FX  InterpZVolZSurfX;
  typedef BasicOpTypes<ZVolField>::InterpC2FY  InterpZVolZSurfY;
  typedef BasicOpTypes<ZVolField>::InterpC2FZ  InterpZVolZSurfZ;


  // staggered volumes to staggered surfaces (advecting velocities)
  typedef OperatorTypeBuilder< Interpolant, XVolField, YSurfXField >::type  InterpXVolYSurfX;
  typedef OperatorTypeBuilder< Interpolant, XVolField, ZSurfXField >::type  InterpXVolZSurfX;

  typedef OperatorTypeBuilder< Interpolant, YVolField, XSurfYField >::type  InterpYVolXSurfY;
  typedef OperatorTypeBuilder< Interpolant, YVolField, ZSurfYField >::type  InterpYVolZSurfY;

  typedef OperatorTypeBuilder< Interpolant, ZVolField, XSurfZField >::type  InterpZVolXSurfZ;
  typedef OperatorTypeBuilder< Interpolant, ZVolField, YSurfZField >::type  InterpZVolYSurfZ;

  // scalar surface to staggered volumes (pressure gradients)
  typedef OperatorTypeBuilder< Interpolant, SSurfXField, XVolField >::type  InterpSSurfXXVol;
  typedef OperatorTypeBuilder< Interpolant, SSurfYField, YVolField >::type  InterpSSurfYYVol;
  typedef OperatorTypeBuilder< Interpolant, SSurfZField, ZVolField >::type  InterpSSurfZZVol;

  /*
   *  NOTE: for UNIFORM MESHES, the following DEGENERACIES exist in
   *  the interpolant operators:
   *
   *    InterpXVolSSurfX = NULL (identity operator)
   *    InterpXVolYSurfX = InterpXVolXSurfY
   *    InterpXVolZSurfX = InterpXVolXSurfZ
   *
   *    InterpYVolSSurfY = NULL (identity operator)
   *    InterpYVolXSurfY = InterpYVolYSurfX
   *    InterpYVolZSurfY = InterpYVolYSurfZ
   *
   *    InterpZVolSSurfZ = NULL (identity operator)
   *    InterpZVolXSurfZ = InterpZVolZSurfX
   *    InterpZVolYSurfZ = InterpZVolZSurfY
   *
   *  However, while since the ghosting may be different on the
   *  different fields, we require individual operators anyway.
   *
   *  For nonuniform meshes, these are all unique operators that must
   *  be independently defined.
   */




  //-- Gradient Operators --//

  typedef BasicOpTypes<SVolField>::GradX  GradSVolSSurfX;  ///< Used to form convection term in scalar equation
  typedef BasicOpTypes<SVolField>::GradY  GradSVolSSurfY;  ///< Used to form convection term in scalar equation
  typedef BasicOpTypes<SVolField>::GradZ  GradSVolSSurfZ;  ///< Used to form convection term in scalar equation

  typedef BasicOpTypes<XVolField>::GradX  GradXVolXSurfX;  ///< Used to form convection term in x-momentum equation
  typedef BasicOpTypes<XVolField>::GradY  GradXVolXSurfY;  ///< Used to form convection term in x-momentum equation
  typedef BasicOpTypes<XVolField>::GradZ  GradXVolXSurfZ;  ///< Used to form convection term in x-momentum equation

  typedef BasicOpTypes<YVolField>::GradX  GradYVolYSurfX;  ///< Used to form convection term in y-momentum equation
  typedef BasicOpTypes<YVolField>::GradY  GradYVolYSurfY;  ///< Used to form convection term in y-momentum equation
  typedef BasicOpTypes<YVolField>::GradZ  GradYVolYSurfZ;  ///< Used to form convection term in y-momentum equation

  typedef BasicOpTypes<ZVolField>::GradX  GradZVolZSurfX;  ///< Used to form convection term in z-momentum equation
  typedef BasicOpTypes<ZVolField>::GradY  GradZVolZSurfY;  ///< Used to form convection term in z-momentum equation
  typedef BasicOpTypes<ZVolField>::GradZ  GradZVolZSurfZ;  ///< Used to form convection term in z-momentum equation

  typedef OperatorTypeBuilder< Gradient, SVolField, XVolField   >::type  GradSVolXVol;    ///< Used to form pressure gradient for x-momentum equation
  typedef OperatorTypeBuilder< Gradient, SVolField, YVolField   >::type  GradSVolYVol;    ///< Used to form pressure gradient for y-momentum equation
  typedef OperatorTypeBuilder< Gradient, SVolField, ZVolField   >::type  GradSVolZVol;    ///< Used to form pressure gradient for z-momentum equation

  typedef OperatorTypeBuilder< Gradient, XVolField, YSurfXField >::type  GradXVolYSurfX;  ///< Used to form \f$\frac{\partial u}{\partial y}\f$ for \f$\tau_{yx}\f$
  typedef OperatorTypeBuilder< Gradient, XVolField, ZSurfXField >::type  GradXVolZSurfX;  ///< Used to form \f$\frac{\partial u}{\partial z}\f$ for \f$\tau_{zx}\f$

  typedef OperatorTypeBuilder< Gradient, YVolField, XSurfYField >::type  GradYVolXSurfY;  ///< Used to form \f$\frac{\partial v}{\partial x}\f$ for \f$\tau_{xy}\f$
  typedef OperatorTypeBuilder< Gradient, YVolField, ZSurfYField >::type  GradYVolZSurfY;  ///< Used to form \f$\frac{\partial v}{\partial z}\f$ for \f$\tau_{yz}\f$

  typedef OperatorTypeBuilder< Gradient, ZVolField, XSurfZField >::type  GradZVolXSurfZ;  ///< Used to form \f$\frac{\partial w}{\partial x}\f$ for \f$\tau_{zx}\f$
  typedef OperatorTypeBuilder< Gradient, ZVolField, YSurfZField >::type  GradZVolYSurfZ;  ///< Used to form \f$\frac{\partial w}{\partial y}\f$ for \f$\tau_{zy}\f$

  typedef OperatorTypeBuilder< Gradient, XVolField, SVolField   >::type  GradXVolSVol;    ///< Used to form the dilatation at scalar CV centers.
  typedef OperatorTypeBuilder< Gradient, YVolField, SVolField   >::type  GradYVolSVol;    ///< Used to form the dilatation at scalar CV centers.
  typedef OperatorTypeBuilder< Gradient, ZVolField, SVolField   >::type  GradZVolSVol;    ///< Used to form the dilatation at scalar CV centers.

  //-- Divergence Operators --//

  typedef BasicOpTypes<SVolField>::DivX  DivSSurfXSVol;
  typedef BasicOpTypes<SVolField>::DivY  DivSSurfYSVol;
  typedef BasicOpTypes<SVolField>::DivZ  DivSSurfZSVol;

  typedef BasicOpTypes<XVolField>::DivX  DivXSurfXXVol;
  typedef BasicOpTypes<XVolField>::DivY  DivXSurfYXVol;
  typedef BasicOpTypes<XVolField>::DivZ  DivXSurfZXVol;

  typedef BasicOpTypes<YVolField>::DivX  DivYSurfXYVol;
  typedef BasicOpTypes<YVolField>::DivY  DivYSurfYYVol;
  typedef BasicOpTypes<YVolField>::DivZ  DivYSurfZYVol;

  typedef BasicOpTypes<ZVolField>::DivX  DivZSurfXZVol;
  typedef BasicOpTypes<ZVolField>::DivY  DivZSurfYZVol;
  typedef BasicOpTypes<ZVolField>::DivZ  DivZSurfZZVol;


  //-- Scratch Operators --//

  typedef OperatorTypeBuilder< Scratch, SVolField, SVolField >::type  ScratchSVol; ///< Used for forming operators such as laplacian
  typedef OperatorTypeBuilder< Scratch, XVolField, XVolField >::type  ScratchXVol; ///< Used for forming operators such as laplacian
  typedef OperatorTypeBuilder< Scratch, YVolField, YVolField >::type  ScratchYVol; ///< Used for forming operators such as laplacian
  typedef OperatorTypeBuilder< Scratch, ZVolField, ZVolField >::type  ScratchZVol; ///< Used for forming operators such as laplacian

  // jcs need to think about these (and test them)
  typedef OperatorTypeBuilder< Scratch, SVolField, SSurfXField >::type ScratchSVolSSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, SVolField, SSurfYField >::type ScratchSVolSSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, SVolField, SSurfZField >::type ScratchSVolSSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef OperatorTypeBuilder< Scratch, XVolField, XSurfXField >::type ScratchXVolXSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, XVolField, XSurfYField >::type ScratchXVolXSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, XVolField, XSurfZField >::type ScratchXVolXSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef OperatorTypeBuilder< Scratch, YVolField, YSurfXField >::type ScratchYVolYSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, YVolField, YSurfYField >::type ScratchYVolYSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, YVolField, YSurfZField >::type ScratchYVolYSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef OperatorTypeBuilder< Scratch, ZVolField, ZSurfXField >::type ScratchZVolZSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, ZVolField, ZSurfYField >::type ScratchZVolZSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef OperatorTypeBuilder< Scratch, ZVolField, ZSurfZField >::type ScratchZVolZSurfZ; ///< Used for modified gradient-shaped operators on scalar CV


  // jcs need to test these.
  typedef OperatorTypeBuilder< Scratch, SSurfXField, SSurfXField >::type  ScratchSSurfX;
  typedef OperatorTypeBuilder< Scratch, SSurfYField, SSurfYField >::type  ScratchSSurfY;
  typedef OperatorTypeBuilder< Scratch, SSurfZField, SSurfZField >::type  ScratchSSurfZ;

  typedef OperatorTypeBuilder< Scratch, XSurfXField, XSurfXField >::type  ScratchXSurfX;
  typedef OperatorTypeBuilder< Scratch, XSurfYField, XSurfYField >::type  ScratchXSurfY;
  typedef OperatorTypeBuilder< Scratch, XSurfZField, XSurfZField >::type  ScratchXSurfZ;

  typedef OperatorTypeBuilder< Scratch, YSurfXField, YSurfXField >::type  ScratchYSurfX;
  typedef OperatorTypeBuilder< Scratch, YSurfYField, YSurfYField >::type  ScratchYSurfY;
  typedef OperatorTypeBuilder< Scratch, YSurfZField, YSurfZField >::type  ScratchYSurfZ;

  typedef OperatorTypeBuilder< Scratch, ZSurfXField, ZSurfXField >::type  ScratchZSurfX;
  typedef OperatorTypeBuilder< Scratch, ZSurfYField, ZSurfYField >::type  ScratchZSurfY;
  typedef OperatorTypeBuilder< Scratch, ZSurfZField, ZSurfZField >::type  ScratchZSurfZ;


  //-- Filter Operators --//
  typedef OperatorTypeBuilder< Filter, SVolField, SVolField >::type FilterSVol;

  //-- Restriction Operators --//
  typedef OperatorTypeBuilder< Restriction, SVolField, SVolField >::type RestrictSVol;


  /**
   *  @}
   *  @}
   */


} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_FVStaggeredOps_h
