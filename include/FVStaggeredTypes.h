#ifndef FVStaggeredTypes_h
#define FVStaggeredTypes_h

#include <SpatialField.h>
#include <SpatialOperator.h>
#include <LinAlgTrilinos.h>

namespace SpatialOps{
namespace FVStaggered{


  struct SVol  { enum{IsSurface=0};  enum{IsVector=0};  typedef NODIR Dir;  typedef NODIR StagDir; };
  struct SSurf { enum{IsSurface=1};  enum{IsVector=0};  typedef NODIR Dir;  typedef NODIR StagDir; };
  struct SSurfX{ enum{IsSurface=1};  enum{IsVector=1};  typedef XDIR  Dir;  typedef NODIR StagDir; };
  struct SSurfY{ enum{IsSurface=1};  enum{IsVector=1};  typedef YDIR  Dir;  typedef NODIR StagDir; };
  struct SSurfZ{ enum{IsSurface=1};  enum{IsVector=1};  typedef ZDIR  Dir;  typedef NODIR StagDir; };

  struct XVol  { enum{IsSurface=0};  enum{IsVector=0};  typedef NODIR Dir;  typedef XDIR  StagDir; };
  struct XSurf { enum{IsSurface=1};  enum{IsVector=0};  typedef NODIR Dir;  typedef XDIR  StagDir; };
  struct XSurfX{ enum{IsSurface=1};  enum{IsVector=1};  typedef XDIR  Dir;  typedef XDIR  StagDir; };
  struct XSurfY{ enum{IsSurface=1};  enum{IsVector=1};  typedef YDIR  Dir;  typedef XDIR  StagDir; };
  struct XSurfZ{ enum{IsSurface=1};  enum{IsVector=1};  typedef ZDIR  Dir;  typedef XDIR  StagDir; };

  struct YVol  { enum{IsSurface=0};  enum{IsVector=0};  typedef NODIR Dir;  typedef YDIR  StagDir; };
  struct YSurf { enum{IsSurface=1};  enum{IsVector=0};  typedef NODIR Dir;  typedef YDIR  StagDir; };
  struct YSurfX{ enum{IsSurface=1};  enum{IsVector=1};  typedef XDIR  Dir;  typedef YDIR  StagDir; };
  struct YSurfY{ enum{IsSurface=1};  enum{IsVector=1};  typedef YDIR  Dir;  typedef YDIR  StagDir; };
  struct YSurfZ{ enum{IsSurface=1};  enum{IsVector=1};  typedef ZDIR  Dir;  typedef YDIR  StagDir; };

  struct ZVol  { enum{IsSurface=0};  enum{IsVector=0};  typedef NODIR Dir;  typedef ZDIR  StagDir; };
  struct ZSurf { enum{IsSurface=1};  enum{IsVector=0};  typedef NODIR Dir;  typedef ZDIR  StagDir; };
  struct ZSurfX{ enum{IsSurface=1};  enum{IsVector=1};  typedef XDIR  Dir;  typedef ZDIR  StagDir; };
  struct ZSurfY{ enum{IsSurface=1};  enum{IsVector=1};  typedef YDIR  Dir;  typedef ZDIR  StagDir; };
  struct ZSurfZ{ enum{IsSurface=1};  enum{IsVector=1};  typedef ZDIR  Dir;  typedef ZDIR  StagDir; };

  struct DefaultGhost{ enum{ NM=1, NP=1 }; };
  struct NoGhost     { enum{ NM=0, NP=0 }; };

  struct Interpolant{};
  struct Gradient{};
  struct Divergence{};
  struct Scratch{};


  //-- Field Types --//

  typedef SpatialField< LinAlgTrilinos, SVol,   DefaultGhost > SVolField;
  typedef SpatialField< LinAlgTrilinos, SSurf,  NoGhost      > SSurfField;
  typedef SpatialField< LinAlgTrilinos, SSurfX, NoGhost      > SSurfXField;
  typedef SpatialField< LinAlgTrilinos, SSurfY, NoGhost      > SSurfYField;
  typedef SpatialField< LinAlgTrilinos, SSurfZ, NoGhost      > SSurfZField;
  typedef SpatialField< LinAlgTrilinos, SVol,   NoGhost      > SVolRHS;

  typedef SpatialField< LinAlgTrilinos, XVol,   DefaultGhost > XVolField;
  typedef SpatialField< LinAlgTrilinos, XSurf,  NoGhost      > XSurfField;
  typedef SpatialField< LinAlgTrilinos, XSurfX, NoGhost      > XSurfXField;
  typedef SpatialField< LinAlgTrilinos, XSurfY, NoGhost      > XSurfYField;
  typedef SpatialField< LinAlgTrilinos, XSurfZ, NoGhost      > XSurfZField;
  typedef SpatialField< LinAlgTrilinos, XVol,   NoGhost      > XVolRHS;

  typedef SpatialField< LinAlgTrilinos, YVol,   DefaultGhost > YVolField;
  typedef SpatialField< LinAlgTrilinos, YSurf,  NoGhost      > YSurfField;
  typedef SpatialField< LinAlgTrilinos, YSurfX, NoGhost      > YSurfXField;
  typedef SpatialField< LinAlgTrilinos, YSurfY, NoGhost      > YSurfYField;
  typedef SpatialField< LinAlgTrilinos, YSurfZ, NoGhost      > YSurfZField;
  typedef SpatialField< LinAlgTrilinos, YVol,   NoGhost      > YVolRHS;

  typedef SpatialField< LinAlgTrilinos, ZVol,   DefaultGhost > ZVolField;
  typedef SpatialField< LinAlgTrilinos, ZSurf,  NoGhost      > ZSurfField;
  typedef SpatialField< LinAlgTrilinos, ZSurfX, NoGhost      > ZSurfXField;
  typedef SpatialField< LinAlgTrilinos, ZSurfY, NoGhost      > ZSurfYField;
  typedef SpatialField< LinAlgTrilinos, ZSurfZ, NoGhost      > ZSurfZField;
  typedef SpatialField< LinAlgTrilinos, ZVol,   NoGhost      > ZVolRHS;



  //-- Interpolant Operators --//

  // scalar volume to scalar surfaces (diffusion coefficients)
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, SSurfXField >  InterpSVolSSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, SSurfYField >  InterpSVolSSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, SSurfZField >  InterpSVolSSurfZ;

  // scalar volume to staggered surfaces (viscosity, dilatation)	      	      	      
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, XSurfXField  >  InterpSVolXSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, XSurfYField  >  InterpSVolXSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, XSurfZField  >  InterpSVolXSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, YSurfXField  >  InterpSVolYSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, YSurfYField  >  InterpSVolYSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, YSurfZField  >  InterpSVolYSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, ZSurfXField  >  InterpSVolZSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, ZSurfYField  >  InterpSVolZSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, ZSurfZField  >  InterpSVolZSurfZ;


  // scalar volume to staggered volumes (density)
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, XVolField   >  InterpSVolXVol;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, YVolField   >  InterpSVolYVol;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, SVolField, ZVolField   >  InterpSVolZVol;


  // staggered volumes to scalar surface (advecting velocities, etc.)
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XVolField, SSurfXField >  InterpXVolSSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YVolField, SSurfYField >  InterpYVolSSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZVolField, SSurfZField >  InterpZVolSSurfZ;


  // staggered volumes to staggered surfaces (momentum solution components)
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XVolField, XSurfXField >  InterpXVolXSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XVolField, XSurfYField >  InterpXVolXSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XVolField, XSurfZField >  InterpXVolXSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YVolField, YSurfXField >  InterpYVolYSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YVolField, YSurfYField >  InterpYVolYSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YVolField, YSurfZField >  InterpYVolYSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZVolField, ZSurfXField >  InterpZVolZSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZVolField, ZSurfYField >  InterpZVolZSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZVolField, ZSurfZField >  InterpZVolZSurfZ;


  // staggered volumes to staggered surfaces (advecting velocities)
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XVolField, YSurfXField >  InterpXVolYSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XVolField, ZSurfXField >  InterpXVolZSurfX;

  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YVolField, XSurfYField >  InterpYVolXSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YVolField, ZSurfYField >  InterpYVolZSurfY;

  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZVolField, XSurfZField >  InterpZVolXSurfZ;
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZVolField, YSurfZField >  InterpZVolYSurfZ;

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
   *  different fields, we may require individual operators anyway.
   *
   *  For nonuniform meshes, these are all unique operators that must
   *  be independently defined.
   */




  //-- Gradient Operators --//

  typedef SpatialOperator< LinAlgTrilinos, Gradient, SVolField, SSurfXField >  GradSVolSSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, SVolField, SSurfYField >  GradSVolSSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, SVolField, SSurfZField >  GradSVolSSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Gradient, XVolField, XSurfXField >  GradXVolXSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, XVolField, XSurfYField >  GradXVolXSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, XVolField, XSurfZField >  GradXVolXSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Gradient, YVolField, YSurfXField >  GradYVolYSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, YVolField, YSurfYField >  GradYVolYSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, YVolField, YSurfZField >  GradYVolYSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZVolField, ZSurfXField >  GradZVolZSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZVolField, ZSurfYField >  GradZVolZSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZVolField, ZSurfZField >  GradZVolZSurfZ;


  //-- Divergence Operators --//

  typedef SpatialOperator< LinAlgTrilinos, Divergence, SSurfXField, SVolField >  DivSSurfXSVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, SSurfYField, SVolField >  DivSSurfYSVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, SSurfZField, SVolField >  DivSSurfZSVol;

  typedef SpatialOperator< LinAlgTrilinos, Divergence, XSurfXField, XVolField >  DivXSurfXXVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, XSurfYField, XVolField >  DivXSurfYXVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, XSurfZField, XVolField >  DivXSurfZXVol;

  typedef SpatialOperator< LinAlgTrilinos, Divergence, YSurfXField, YVolField >  DivYSurfXYVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, YSurfYField, YVolField >  DivYSurfYYVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, YSurfZField, YVolField >  DivYSurfZYVol;

  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZSurfXField, ZVolField >  DivZSurfXZVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZSurfYField, ZVolField >  DivZSurfYZVol;
  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZSurfZField, ZVolField >  DivZSurfZZVol;


  //-- Scratch Operators --//

  typedef SpatialOperator< LinAlgTrilinos, Scratch, SVolField, SVolField >  ScratchSVol; ///< Used for forming operators such as laplacian
  typedef SpatialOperator< LinAlgTrilinos, Scratch, XVolField, XVolField >  ScratchXVol; ///< Used for forming operators such as laplacian
  typedef SpatialOperator< LinAlgTrilinos, Scratch, YVolField, YVolField >  ScratchYVol; ///< Used for forming operators such as laplacian
  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZVolField, ZVolField >  ScratchZVol; ///< Used for forming operators such as laplacian

  // jcs need to think about these (and test them)
  typedef SpatialOperator< LinAlgTrilinos, Scratch, SVolField, SSurfXField > ScratchSVolSSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, SVolField, SSurfYField > ScratchSVolSSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, SVolField, SSurfZField > ScratchSVolSSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef SpatialOperator< LinAlgTrilinos, Scratch, XVolField, XSurfXField > ScratchXVolXSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, XVolField, XSurfYField > ScratchXVolXSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, XVolField, XSurfZField > ScratchXVolXSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef SpatialOperator< LinAlgTrilinos, Scratch, YVolField, YSurfXField > ScratchYVolYSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, YVolField, YSurfYField > ScratchYVolYSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, YVolField, YSurfZField > ScratchYVolYSurfZ; ///< Used for modified gradient-shaped operators on scalar CV

  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZVolField, ZSurfXField > ScratchZVolZSurfX; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZVolField, ZSurfYField > ScratchZVolZSurfY; ///< Used for modified gradient-shaped operators on scalar CV
  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZVolField, ZSurfZField > ScratchZVolZSurfZ; ///< Used for modified gradient-shaped operators on scalar CV


  // jcs need to test these.
  typedef SpatialOperator< LinAlgTrilinos, Scratch, SSurfXField, SSurfXField >  ScratchSSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, SSurfYField, SSurfYField >  ScratchSSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, SSurfZField, SSurfZField >  ScratchSSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Scratch, XSurfXField, XSurfXField >  ScratchXSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, XSurfYField, XSurfYField >  ScratchXSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, XSurfZField, XSurfZField >  ScratchXSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Scratch, YSurfXField, YSurfXField >  ScratchYSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, YSurfYField, YSurfYField >  ScratchYSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, YSurfZField, YSurfZField >  ScratchYSurfZ;

  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZSurfXField, ZSurfXField >  ScratchZSurfX;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZSurfYField, ZSurfYField >  ScratchZSurfY;
  typedef SpatialOperator< LinAlgTrilinos, Scratch, ZSurfZField, ZSurfZField >  ScratchZSurfZ;


}// namespace FVStaggered
}// namespace SpatialOps

#endif
