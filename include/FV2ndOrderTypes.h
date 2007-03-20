#ifndef UT_FV2ndOrderTypes_h
#define UT_FV2ndOrderTypes_h

/**
 *  This header is intended for inclusion in second order finite
 *  volume schemes.  It defines field and operator types required for
 *  second order finite volume spatial discretizations.
 *
 *
 */

#include <FVStaggeredSpatialOps.h>
#include <LinAlgTrilinos.h>

namespace SpatialOps{
namespace FVStaggeredUniform{



  //==================================================================
  // Spatial Field Traits

  typedef FieldTraits< Side, DefaultSideGhostingX >     XSideFieldTraits;
  typedef FieldTraits< Side, DefaultSideGhostingY >     YSideFieldTraits;
  typedef FieldTraits< Side, DefaultSideGhostingZ >     ZSideFieldTraits;

  typedef FieldTraits< Cell, DefaultCellGhosting >      CellFieldTraits;

  typedef FieldTraits< Cell, NoGhosting >               NoGhostCellFieldTraits;
  typedef FieldTraits< Side, NoGhosting >               NoGhostSideFieldTraits;

  // Spatial Field Traits
  //==================================================================




  //==================================================================
  // Spatial Field Types

  // side fields for x, y, and z.  These have ghosting of 1 on (-) side and 2 on (+) side
  typedef SpatialField< LinAlgTrilinos, XSideFieldTraits::StorageLocation, XSideFieldTraits::GhostTraits >   XSideField;
  typedef SpatialField< LinAlgTrilinos, YSideFieldTraits::StorageLocation, YSideFieldTraits::GhostTraits >   YSideField;
  typedef SpatialField< LinAlgTrilinos, ZSideFieldTraits::StorageLocation, ZSideFieldTraits::GhostTraits >   ZSideField;

  // cell fields.  These have ghosting of 1 on both sides.
  typedef SpatialField< LinAlgTrilinos, CellFieldTraits::StorageLocation,   CellFieldTraits::GhostTraits >   CellField;

  // side and cell fields without ghosting.
  typedef SpatialField< LinAlgTrilinos, NoGhostCellFieldTraits::StorageLocation, NoGhostCellFieldTraits::GhostTraits >   CellFieldNoGhost;
  typedef SpatialField< LinAlgTrilinos, NoGhostSideFieldTraits::StorageLocation, NoGhostSideFieldTraits::GhostTraits >   SideFieldNoGhost;

  // Spatial Field Types
  //==================================================================




  //==================================================================
  // Operator Assembler types

  // linear interpolant assembler - Cell to Side
  typedef LinearInterpolantAssembler< XDIR, Cell, CellFieldTraits::GhostTraits, XSideFieldTraits::GhostTraits >  InterpXC2FAssembler;
  typedef LinearInterpolantAssembler< YDIR, Cell, CellFieldTraits::GhostTraits, YSideFieldTraits::GhostTraits >  InterpYC2FAssembler;
  typedef LinearInterpolantAssembler< ZDIR, Cell, CellFieldTraits::GhostTraits, ZSideFieldTraits::GhostTraits >  InterpZC2FAssembler;

  // linear interpolant assembler - Side to Cell
  typedef LinearInterpolantAssembler< XDIR, Side, XSideFieldTraits::GhostTraits, CellFieldTraits::GhostTraits >  InterpXF2CAssembler;
  typedef LinearInterpolantAssembler< YDIR, Side, YSideFieldTraits::GhostTraits, CellFieldTraits::GhostTraits >  InterpYF2CAssembler;
  typedef LinearInterpolantAssembler< ZDIR, Side, ZSideFieldTraits::GhostTraits, CellFieldTraits::GhostTraits >  InterpZF2CAssembler;

  // divergence assembler - Side to Cell
  typedef DivergenceAssembler< XDIR, Side, XSideFieldTraits::GhostTraits, CellFieldTraits ::GhostTraits >        DivXF2CAssembler;
  typedef DivergenceAssembler< YDIR, Side, YSideFieldTraits::GhostTraits, CellFieldTraits ::GhostTraits >        DivYF2CAssembler;
  typedef DivergenceAssembler< ZDIR, Side, ZSideFieldTraits::GhostTraits, CellFieldTraits ::GhostTraits >        DivZF2CAssembler;

  // divergence assembler - Cell to Side
  typedef DivergenceAssembler< XDIR, Cell, CellFieldTraits ::GhostTraits, XSideFieldTraits::GhostTraits >        DivXC2FAssembler;
  typedef DivergenceAssembler< YDIR, Cell, CellFieldTraits ::GhostTraits, YSideFieldTraits::GhostTraits >        DivYC2FAssembler;
  typedef DivergenceAssembler< ZDIR, Cell, CellFieldTraits ::GhostTraits, ZSideFieldTraits::GhostTraits >        DivZC2FAssembler;
													         
  // gradient assembler - Cell to Side
  typedef GradientAssembler< XDIR, Cell, CellFieldTraits ::GhostTraits, XSideFieldTraits::GhostTraits >          GradXC2FAssembler;
  typedef GradientAssembler< YDIR, Cell, CellFieldTraits ::GhostTraits, YSideFieldTraits::GhostTraits >          GradYC2FAssembler;
  typedef GradientAssembler< ZDIR, Cell, CellFieldTraits ::GhostTraits, ZSideFieldTraits::GhostTraits >          GradZC2FAssembler;
													       
  // gradient assembler - Side to Cell
  typedef GradientAssembler< XDIR, Side, XSideFieldTraits::GhostTraits, CellFieldTraits::GhostTraits  >          GradXF2CAssembler;
  typedef GradientAssembler< YDIR, Side, YSideFieldTraits::GhostTraits, CellFieldTraits::GhostTraits  >          GradYF2CAssembler;
  typedef GradientAssembler< ZDIR, Side, ZSideFieldTraits::GhostTraits, CellFieldTraits::GhostTraits  >          GradZF2CAssembler;
													       
  // Operator Assembler types
  //==================================================================




  //==================================================================
  // Operator Types

  // linear interpolants - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, XSideFieldTraits, InterpXC2FAssembler >  InterpXC2F;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, YSideFieldTraits, InterpYC2FAssembler >  InterpYC2F;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, ZSideFieldTraits, InterpZC2FAssembler >  InterpZC2F;

  // linear interpolants - Side to Cell
  typedef SpatialOperator< LinAlgTrilinos, XSideFieldTraits, CellFieldTraits, InterpXF2CAssembler >  InterpXF2C;
  typedef SpatialOperator< LinAlgTrilinos, YSideFieldTraits, CellFieldTraits, InterpYF2CAssembler >  InterpYF2C;
  typedef SpatialOperator< LinAlgTrilinos, ZSideFieldTraits, CellFieldTraits, InterpZF2CAssembler >  InterpZF2C;

  // divergence operators - Side to Cell
  typedef SpatialOperator< LinAlgTrilinos, XSideFieldTraits, CellFieldTraits, DivXF2CAssembler >     DivXF2C;
  typedef SpatialOperator< LinAlgTrilinos, YSideFieldTraits, CellFieldTraits, DivYF2CAssembler >     DivYF2C;
  typedef SpatialOperator< LinAlgTrilinos, ZSideFieldTraits, CellFieldTraits, DivZF2CAssembler >     DivZF2C;
// divergence operators - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, XSideFieldTraits, DivXC2FAssembler >     DivXC2F;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, YSideFieldTraits, DivYC2FAssembler >     DivYC2F;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, ZSideFieldTraits, DivZC2FAssembler >     DivZC2F;
												        
  // gradient operators - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, XSideFieldTraits, GradXC2FAssembler >    GradXC2F;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, YSideFieldTraits, GradYC2FAssembler >    GradYC2F;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, ZSideFieldTraits, GradZC2FAssembler >    GradZC2F;

  // gradient operators - Side to Cell
  typedef SpatialOperator< LinAlgTrilinos, XSideFieldTraits, CellFieldTraits, GradXF2CAssembler >    GradXF2C;
  typedef SpatialOperator< LinAlgTrilinos, YSideFieldTraits, CellFieldTraits, GradYF2CAssembler >    GradYF2C;
  typedef SpatialOperator< LinAlgTrilinos, ZSideFieldTraits, CellFieldTraits, GradZF2CAssembler >    GradZF2C;

  // Operator Types
  //==================================================================








  //==================================================================
  // SCRATCH OPERATORS and their ASSEMBLERS

  typedef ScratchAssembler< 3, XDIR, Cell, CellFieldTraits::GhostTraits,  CellFieldTraits::GhostTraits > SxCellAssembler;
  typedef ScratchAssembler< 3, YDIR, Cell, CellFieldTraits::GhostTraits,  CellFieldTraits::GhostTraits > SyCellAssembler;
  typedef ScratchAssembler< 3, ZDIR, Cell, CellFieldTraits::GhostTraits,  CellFieldTraits::GhostTraits > SzCellAssembler;

  typedef ScratchAssembler< 2, XDIR, Cell, CellFieldTraits::GhostTraits, XSideFieldTraits::GhostTraits > SxC2FAssembler;
  typedef ScratchAssembler< 2, YDIR, Cell, CellFieldTraits::GhostTraits, YSideFieldTraits::GhostTraits > SyC2FAssembler;
  typedef ScratchAssembler< 2, ZDIR, Cell, CellFieldTraits::GhostTraits, ZSideFieldTraits::GhostTraits > SzC2FAssembler;



  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits,  CellFieldTraits, SxCellAssembler > SxCell;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits,  CellFieldTraits, SyCellAssembler > SyCell;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits,  CellFieldTraits, SzCellAssembler > SzCell;

  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, XSideFieldTraits, SxC2FAssembler  > SxCellSide;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, YSideFieldTraits, SyC2FAssembler  > SyCellSide;
  typedef SpatialOperator< LinAlgTrilinos, CellFieldTraits, ZSideFieldTraits, SzC2FAssembler  > SzCellSide;

  // SCRATCH OPERATORS and their ASSEMBLERS
  //==================================================================



} // namespace FVStaggeredUniform
} // namespace SpatialOps


#endif // UT_FV2ndOrderTypes_h
