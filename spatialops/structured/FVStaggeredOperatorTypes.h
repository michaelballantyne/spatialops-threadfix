#ifndef SpatialOps_FVStaggeredOps_h
#define SpatialOps_FVStaggeredOps_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOperator.h>
#include "FVStaggeredFieldTypes.h"

namespace SpatialOps{
namespace structured{

  /**
   *  \file FVStaggeredOps.h
   *
   *  \addtogroup structured
   *  @{
   *  \addtogroup operators
   *  @{
   */

  //-- Interpolant Operators --//

  // scalar volume to scalar surfaces (diffusion coefficients)
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, SSurfXField >  InterpSVolSSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, SSurfYField >  InterpSVolSSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, SSurfZField >  InterpSVolSSurfZ;

  // scalar volume to staggered surfaces (viscosity, dilatation)	      	      	      
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XSurfXField >  InterpSVolXSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XSurfYField >  InterpSVolXSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XSurfZField >  InterpSVolXSurfZ;

  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YSurfXField >  InterpSVolYSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YSurfYField >  InterpSVolYSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YSurfZField >  InterpSVolYSurfZ;

  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZSurfXField >  InterpSVolZSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZSurfYField >  InterpSVolZSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZSurfZField >  InterpSVolZSurfZ;


  // scalar volume to staggered volumes (density)
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XVolField   >  InterpSVolXVol;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YVolField   >  InterpSVolYVol;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZVolField   >  InterpSVolZVol;


  // staggered volumes to scalar surface (advecting velocities, etc.)
  typedef SpatialOperator< LinAlg, Interpolant, XVolField, SSurfXField >  InterpXVolSSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, YVolField, SSurfYField >  InterpYVolSSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, ZVolField, SSurfZField >  InterpZVolSSurfZ;


  // staggered volumes to staggered surfaces (momentum solution components)
  typedef SpatialOperator< LinAlg, Interpolant, XVolField, XSurfXField >  InterpXVolXSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, XVolField, XSurfYField >  InterpXVolXSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, XVolField, XSurfZField >  InterpXVolXSurfZ;

  typedef SpatialOperator< LinAlg, Interpolant, YVolField, YSurfXField >  InterpYVolYSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, YVolField, YSurfYField >  InterpYVolYSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, YVolField, YSurfZField >  InterpYVolYSurfZ;

  typedef SpatialOperator< LinAlg, Interpolant, ZVolField, ZSurfXField >  InterpZVolZSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, ZVolField, ZSurfYField >  InterpZVolZSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, ZVolField, ZSurfZField >  InterpZVolZSurfZ;


  // staggered volumes to staggered surfaces (advecting velocities)
  typedef SpatialOperator< LinAlg, Interpolant, XVolField, YSurfXField >  InterpXVolYSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, XVolField, ZSurfXField >  InterpXVolZSurfX;

  typedef SpatialOperator< LinAlg, Interpolant, YVolField, XSurfYField >  InterpYVolXSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, YVolField, ZSurfYField >  InterpYVolZSurfY;

  typedef SpatialOperator< LinAlg, Interpolant, ZVolField, XSurfZField >  InterpZVolXSurfZ;
  typedef SpatialOperator< LinAlg, Interpolant, ZVolField, YSurfZField >  InterpZVolYSurfZ;

  // scalar surface to staggered volumes (pressure gradients)
  typedef SpatialOperator< LinAlg, Interpolant, SSurfXField, XVolField >  InterpSSurfXXVol;
  typedef SpatialOperator< LinAlg, Interpolant, SSurfYField, YVolField >  InterpSSurfYYVol;
  typedef SpatialOperator< LinAlg, Interpolant, SSurfZField, ZVolField >  InterpSSurfZZVol;

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

  typedef SpatialOperator< LinAlg, Gradient, SVolField, SSurfXField >  GradSVolSSurfX;  ///< Used to form convection term in scalar equation
  typedef SpatialOperator< LinAlg, Gradient, SVolField, SSurfYField >  GradSVolSSurfY;  ///< Used to form convection term in scalar equation
  typedef SpatialOperator< LinAlg, Gradient, SVolField, SSurfZField >  GradSVolSSurfZ;  ///< Used to form convection term in scalar equation

  typedef SpatialOperator< LinAlg, Gradient, XVolField, XSurfXField >  GradXVolXSurfX;  ///< Used to form convection term in x-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, XVolField, XSurfYField >  GradXVolXSurfY;  ///< Used to form convection term in x-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, XVolField, XSurfZField >  GradXVolXSurfZ;  ///< Used to form convection term in x-momentum equation

  typedef SpatialOperator< LinAlg, Gradient, YVolField, YSurfXField >  GradYVolYSurfX;  ///< Used to form convection term in y-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, YVolField, YSurfYField >  GradYVolYSurfY;  ///< Used to form convection term in y-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, YVolField, YSurfZField >  GradYVolYSurfZ;  ///< Used to form convection term in y-momentum equation

  typedef SpatialOperator< LinAlg, Gradient, ZVolField, ZSurfXField >  GradZVolZSurfX;  ///< Used to form convection term in z-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, ZVolField, ZSurfYField >  GradZVolZSurfY;  ///< Used to form convection term in z-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, ZVolField, ZSurfZField >  GradZVolZSurfZ;  ///< Used to form convection term in z-momentum equation

  typedef SpatialOperator< LinAlg, Gradient, SVolField, XVolField   >  GradSVolXVol;    ///< Used to form pressure gradient for x-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, SVolField, YVolField   >  GradSVolYVol;    ///< Used to form pressure gradient for y-momentum equation
  typedef SpatialOperator< LinAlg, Gradient, SVolField, ZVolField   >  GradSVolZVol;    ///< Used to form pressure gradient for z-momentum equation

  typedef SpatialOperator< LinAlg, Gradient, XVolField, YSurfXField >  GradXVolYSurfX;  ///< Used to form \f$\frac{\partial u}{\partial y}\f$ for \f$\tau_{yx}\f$
  typedef SpatialOperator< LinAlg, Gradient, XVolField, ZSurfXField >  GradXVolZSurfX;  ///< Used to form \f$\frac{\partial u}{\partial z}\f$ for \f$\tau_{zx}\f$

  typedef SpatialOperator< LinAlg, Gradient, YVolField, XSurfYField >  GradYVolXSurfY;  ///< Used to form \f$\frac{\partial v}{\partial x}\f$ for \f$\tau_{xy}\f$
  typedef SpatialOperator< LinAlg, Gradient, YVolField, ZSurfYField >  GradYVolZSurfY;  ///< Used to form \f$\frac{\partial v}{\partial z}\f$ for \f$\tau_{yz}\f$

  typedef SpatialOperator< LinAlg, Gradient, ZVolField, XSurfZField >  GradZVolXSurfZ;  ///< Used to form \f$\frac{\partial w}{\partial x}\f$ for \f$\tau_{zx}\f$
  typedef SpatialOperator< LinAlg, Gradient, ZVolField, YSurfZField >  GradZVolYSurfZ;  ///< Used to form \f$\frac{\partial w}{\partial y}\f$ for \f$\tau_{zy}\f$

  typedef SpatialOperator< LinAlg, Gradient, XVolField, SVolField   >  GradXVolSVol;    ///< Used to form the dilatation at scalar CV centers.
  typedef SpatialOperator< LinAlg, Gradient, YVolField, SVolField   >  GradYVolSVol;    ///< Used to form the dilatation at scalar CV centers.
  typedef SpatialOperator< LinAlg, Gradient, ZVolField, SVolField   >  GradZVolSVol;    ///< Used to form the dilatation at scalar CV centers.

  //-- Divergence Operators --//

  typedef SpatialOperator< LinAlg, Divergence, SSurfXField, SVolField >  DivSSurfXSVol;
  typedef SpatialOperator< LinAlg, Divergence, SSurfYField, SVolField >  DivSSurfYSVol;
  typedef SpatialOperator< LinAlg, Divergence, SSurfZField, SVolField >  DivSSurfZSVol;

  typedef SpatialOperator< LinAlg, Divergence, XSurfXField, XVolField >  DivXSurfXXVol;
  typedef SpatialOperator< LinAlg, Divergence, XSurfYField, XVolField >  DivXSurfYXVol;
  typedef SpatialOperator< LinAlg, Divergence, XSurfZField, XVolField >  DivXSurfZXVol;

  typedef SpatialOperator< LinAlg, Divergence, YSurfXField, YVolField >  DivYSurfXYVol;
  typedef SpatialOperator< LinAlg, Divergence, YSurfYField, YVolField >  DivYSurfYYVol;
  typedef SpatialOperator< LinAlg, Divergence, YSurfZField, YVolField >  DivYSurfZYVol;

  typedef SpatialOperator< LinAlg, Divergence, ZSurfXField, ZVolField >  DivZSurfXZVol;
  typedef SpatialOperator< LinAlg, Divergence, ZSurfYField, ZVolField >  DivZSurfYZVol;
  typedef SpatialOperator< LinAlg, Divergence, ZSurfZField, ZVolField >  DivZSurfZZVol;


  //-- Scratch Operators --//

  typedef SpatialOperator< LinAlg, Scratch, SVolField, SVolField >  ScratchSVol; ///< Used for forming operators such as laplacian
  typedef SpatialOperator< LinAlg, Scratch, XVolField, XVolField >  ScratchXVol; ///< Used for forming operators such as laplacian
  typedef SpatialOperator< LinAlg, Scratch, YVolField, YVolField >  ScratchYVol; ///< Used for forming operators such as laplacian
  typedef SpatialOperator< LinAlg, Scratch, ZVolField, ZVolField >  ScratchZVol; ///< Used for forming operators such as laplacian

  // jcs need to think about these (and test them)
  typedef SpatialOperator< LinAlg, Scratch, SVolField, SSurfXField > ScratchSVolSSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, SVolField, SSurfYField > ScratchSVolSSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, SVolField, SSurfZField > ScratchSVolSSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef SpatialOperator< LinAlg, Scratch, XVolField, XSurfXField > ScratchXVolXSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, XVolField, XSurfYField > ScratchXVolXSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, XVolField, XSurfZField > ScratchXVolXSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef SpatialOperator< LinAlg, Scratch, YVolField, YSurfXField > ScratchYVolYSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, YVolField, YSurfYField > ScratchYVolYSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, YVolField, YSurfZField > ScratchYVolYSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef SpatialOperator< LinAlg, Scratch, ZVolField, ZSurfXField > ScratchZVolZSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, ZVolField, ZSurfYField > ScratchZVolZSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlg, Scratch, ZVolField, ZSurfZField > ScratchZVolZSurfZ; ///< Used for modified gradient-shaped operators on scalar CV


  // jcs need to test these.
  typedef SpatialOperator< LinAlg, Scratch, SSurfXField, SSurfXField >  ScratchSSurfX;
  typedef SpatialOperator< LinAlg, Scratch, SSurfYField, SSurfYField >  ScratchSSurfY;
  typedef SpatialOperator< LinAlg, Scratch, SSurfZField, SSurfZField >  ScratchSSurfZ;

  typedef SpatialOperator< LinAlg, Scratch, XSurfXField, XSurfXField >  ScratchXSurfX;
  typedef SpatialOperator< LinAlg, Scratch, XSurfYField, XSurfYField >  ScratchXSurfY;
  typedef SpatialOperator< LinAlg, Scratch, XSurfZField, XSurfZField >  ScratchXSurfZ;

  typedef SpatialOperator< LinAlg, Scratch, YSurfXField, YSurfXField >  ScratchYSurfX;
  typedef SpatialOperator< LinAlg, Scratch, YSurfYField, YSurfYField >  ScratchYSurfY;
  typedef SpatialOperator< LinAlg, Scratch, YSurfZField, YSurfZField >  ScratchYSurfZ;

  typedef SpatialOperator< LinAlg, Scratch, ZSurfXField, ZSurfXField >  ScratchZSurfX;
  typedef SpatialOperator< LinAlg, Scratch, ZSurfYField, ZSurfYField >  ScratchZSurfY;
  typedef SpatialOperator< LinAlg, Scratch, ZSurfZField, ZSurfZField >  ScratchZSurfZ;


  //-- Filter Operators --//
  typedef SpatialOperator< LinAlg, Filter, SVolField, SVolField > FilterSVol;

  //-- Restriction Operators --//
  typedef SpatialOperator< LinAlg, Restriction, SVolField, SVolField > RestrictSVol;


  /**
   *  @}
   *  @}
   */


} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_FVStaggeredOps_h
