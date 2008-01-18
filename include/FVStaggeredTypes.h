#ifndef FVStaggeredTypes_h
#define FVStaggeredTypes_h

#include <SpatialField.h>
#include <SpatialOperator.h>
#include <LinAlgTrilinos.h>
//#include <LinAlgUBlas.h>

typedef SpatialOps::LinAlgTrilinos LinAlg;
//typedef SpatialOps::LinAlgUBlas LinAlg;

namespace SpatialOps{
namespace FVStaggered{


  struct SVol  { enum{IsSurface=0};  typedef NODIR Dir;  typedef NODIR StagDir; };
  struct SSurfX{ enum{IsSurface=1};  typedef XDIR  Dir;  typedef NODIR StagDir; };
  struct SSurfY{ enum{IsSurface=1};  typedef YDIR  Dir;  typedef NODIR StagDir; };
  struct SSurfZ{ enum{IsSurface=1};  typedef ZDIR  Dir;  typedef NODIR StagDir; };

  struct XVol  { enum{IsSurface=0};  typedef NODIR Dir;  typedef XDIR  StagDir; };
  struct XSurfX{ enum{IsSurface=1};  typedef XDIR  Dir;  typedef XDIR  StagDir; };
  struct XSurfY{ enum{IsSurface=1};  typedef YDIR  Dir;  typedef XDIR  StagDir; };
  struct XSurfZ{ enum{IsSurface=1};  typedef ZDIR  Dir;  typedef XDIR  StagDir; };

  struct YVol  { enum{IsSurface=0};  typedef NODIR Dir;  typedef YDIR  StagDir; };
  struct YSurfX{ enum{IsSurface=1};  typedef XDIR  Dir;  typedef YDIR  StagDir; };
  struct YSurfY{ enum{IsSurface=1};  typedef YDIR  Dir;  typedef YDIR  StagDir; };
  struct YSurfZ{ enum{IsSurface=1};  typedef ZDIR  Dir;  typedef YDIR  StagDir; };

  struct ZVol  { enum{IsSurface=0};  typedef NODIR Dir;  typedef ZDIR  StagDir; };
  struct ZSurfX{ enum{IsSurface=1};  typedef XDIR  Dir;  typedef ZDIR  StagDir; };
  struct ZSurfY{ enum{IsSurface=1};  typedef YDIR  Dir;  typedef ZDIR  StagDir; };
  struct ZSurfZ{ enum{IsSurface=1};  typedef ZDIR  Dir;  typedef ZDIR  StagDir; };

  struct DefaultGhost{ enum{ NM=1, NP=1 }; };
  struct NoGhost     { enum{ NM=0, NP=0 }; };

  struct Interpolant{};
  struct Gradient{};
  struct Divergence{};
  struct Scratch{};


  //-- Field Types --//

  typedef SpatialField< LinAlg, SVol,   DefaultGhost > SVolField;
  typedef SpatialField< LinAlg, SSurfX, NoGhost      > SSurfXField;
  typedef SpatialField< LinAlg, SSurfY, NoGhost      > SSurfYField;
  typedef SpatialField< LinAlg, SSurfZ, NoGhost      > SSurfZField;
  typedef SpatialField< LinAlg, SVol,   NoGhost      > SVolRHS;

  typedef SpatialField< LinAlg, XVol,   DefaultGhost > XVolField;
  typedef SpatialField< LinAlg, XSurfX, NoGhost      > XSurfXField;
  typedef SpatialField< LinAlg, XSurfY, NoGhost      > XSurfYField;
  typedef SpatialField< LinAlg, XSurfZ, NoGhost      > XSurfZField;
  typedef SpatialField< LinAlg, XVol,   NoGhost      > XVolRHS;

  typedef SpatialField< LinAlg, YVol,   DefaultGhost > YVolField;
  typedef SpatialField< LinAlg, YSurfX, NoGhost      > YSurfXField;
  typedef SpatialField< LinAlg, YSurfY, NoGhost      > YSurfYField;
  typedef SpatialField< LinAlg, YSurfZ, NoGhost      > YSurfZField;
  typedef SpatialField< LinAlg, YVol,   NoGhost      > YVolRHS;

  typedef SpatialField< LinAlg, ZVol,   DefaultGhost > ZVolField;
  typedef SpatialField< LinAlg, ZSurfX, NoGhost      > ZSurfXField;
  typedef SpatialField< LinAlg, ZSurfY, NoGhost      > ZSurfYField;
  typedef SpatialField< LinAlg, ZSurfZ, NoGhost      > ZSurfZField;
  typedef SpatialField< LinAlg, ZVol,   NoGhost      > ZVolRHS;



  //-- Interpolant Operators --//

  // scalar volume to scalar surfaces (diffusion coefficients)
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, SSurfXField >  InterpSVolSSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, SSurfYField >  InterpSVolSSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, SSurfZField >  InterpSVolSSurfZ;

  // scalar volume to staggered surfaces (viscosity, dilatation)	      	      	      
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XSurfXField  >  InterpSVolXSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XSurfYField  >  InterpSVolXSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, XSurfZField  >  InterpSVolXSurfZ;

  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YSurfXField  >  InterpSVolYSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YSurfYField  >  InterpSVolYSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, YSurfZField  >  InterpSVolYSurfZ;

  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZSurfXField  >  InterpSVolZSurfX;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZSurfYField  >  InterpSVolZSurfY;
  typedef SpatialOperator< LinAlg, Interpolant, SVolField, ZSurfZField  >  InterpSVolZSurfZ;


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

  /*
   *  NOTE: for UNIFORM MESHES, the following DEGENERACIES exist in
   *  the interpolant operators:
   *
   *    InterpXVolSSurfX = NULL  (identity operator)
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

  typedef SpatialOperator< LinAlg, Gradient, SVolField, SSurfXField >  GradSVolSSurfX;
  typedef SpatialOperator< LinAlg, Gradient, SVolField, SSurfYField >  GradSVolSSurfY;
  typedef SpatialOperator< LinAlg, Gradient, SVolField, SSurfZField >  GradSVolSSurfZ;

  typedef SpatialOperator< LinAlg, Gradient, XVolField, XSurfXField >  GradXVolXSurfX;
  typedef SpatialOperator< LinAlg, Gradient, XVolField, XSurfYField >  GradXVolXSurfY;
  typedef SpatialOperator< LinAlg, Gradient, XVolField, XSurfZField >  GradXVolXSurfZ;

  typedef SpatialOperator< LinAlg, Gradient, YVolField, YSurfXField >  GradYVolYSurfX;
  typedef SpatialOperator< LinAlg, Gradient, YVolField, YSurfYField >  GradYVolYSurfY;
  typedef SpatialOperator< LinAlg, Gradient, YVolField, YSurfZField >  GradYVolYSurfZ;

  typedef SpatialOperator< LinAlg, Gradient, ZVolField, ZSurfXField >  GradZVolZSurfX;
  typedef SpatialOperator< LinAlg, Gradient, ZVolField, ZSurfYField >  GradZVolZSurfY;
  typedef SpatialOperator< LinAlg, Gradient, ZVolField, ZSurfZField >  GradZVolZSurfZ;


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


}// namespace FVStaggered
}// namespace SpatialOps

#endif
