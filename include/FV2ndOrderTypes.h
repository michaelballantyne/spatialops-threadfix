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

  /**
   * @brief Defines a Gradient type.  Doesn't need to provide any functionality.
   */
  struct Gradient{};

  /**
   * @brief Defines a Divergence type.  Doesn't need to provide any functionality.
 */
  struct Divergence{};

  /**
   * @brief Defines an Interpolant type.  Doesn't need to provide any functionality.
   */
  struct Interpolant{};

  /**
   *  @brief Defines a Scratch type with a specified number
   *  nonzeros. Doesn't need to provide any functionality.
   */
  template< int N >
  struct Scratch{ static const int NumNonZero = N; };


  //==================================================================
  // Spatial Field Traits

  typedef FieldTraits< Side, DefaultSideGhosting<XDIR> >  XSideFieldTraits;
  typedef FieldTraits< Side, DefaultSideGhosting<YDIR> >  YSideFieldTraits;
  typedef FieldTraits< Side, DefaultSideGhosting<ZDIR> >  ZSideFieldTraits;

  typedef FieldTraits< Cell, DefaultCellGhosting >        CellFieldTraits;

  typedef FieldTraits< Cell, NoGhosting >                 NoGhostCellFieldTraits;
  typedef FieldTraits< Side, NoGhosting >                 NoGhostSideFieldTraits;


  typedef FieldTraits< Edge<XDIR>, DefaultSideGhosting<YDIR> > XEdgeYDirFieldTraits;
  typedef FieldTraits< Edge<XDIR>, DefaultSideGhosting<ZDIR> > XEdgeZDirFieldTraits;

  typedef FieldTraits< Edge<YDIR>, DefaultSideGhosting<XDIR> > YEdgeXDirFieldTraits;
  typedef FieldTraits< Edge<YDIR>, DefaultSideGhosting<ZDIR> > YEdgeZDirFieldTraits;

  typedef FieldTraits< Edge<ZDIR>, DefaultSideGhosting<XDIR> > ZEdgeXDirFieldTraits;
  typedef FieldTraits< Edge<ZDIR>, DefaultSideGhosting<YDIR> > ZEdgeYDirFieldTraits;

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

  // X-edge fields
  typedef SpatialField< LinAlgTrilinos, XEdgeYDirFieldTraits::StorageLocation, XEdgeYDirFieldTraits::GhostTraits >  XEdgeYDirField;
  typedef SpatialField< LinAlgTrilinos, XEdgeZDirFieldTraits::StorageLocation, XEdgeZDirFieldTraits::GhostTraits >  XEdgeZDirField;

  // Y-edge fields
  typedef SpatialField< LinAlgTrilinos, YEdgeXDirFieldTraits::StorageLocation, YEdgeXDirFieldTraits::GhostTraits >  YEdgeXDirField;
  typedef SpatialField< LinAlgTrilinos, YEdgeZDirFieldTraits::StorageLocation, YEdgeZDirFieldTraits::GhostTraits >  YEdgeZDirField;

  // Z-edge fields
  typedef SpatialField< LinAlgTrilinos, ZEdgeXDirFieldTraits::StorageLocation, ZEdgeXDirFieldTraits::GhostTraits >  ZEdgeXDirField;
  typedef SpatialField< LinAlgTrilinos, ZEdgeYDirFieldTraits::StorageLocation, ZEdgeYDirFieldTraits::GhostTraits >  ZEdgeYDirField;

  // Spatial Field Types
  //==================================================================






  //==================================================================
  // Operator Types

  // linear interpolants - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XDIR, CellFieldTraits, XSideFieldTraits > InterpXC2F;  ///< Interpolate cell to face in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YDIR, CellFieldTraits, YSideFieldTraits > InterpYC2F;  ///< Interpolate cell to face in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZDIR, CellFieldTraits, ZSideFieldTraits > InterpZC2F;  ///< Interpolate cell to face in z-dir

  // linear interpolants - Side to Cell
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XDIR, XSideFieldTraits, CellFieldTraits > InterpXF2C;  ///< Interpolate face to cell in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YDIR, YSideFieldTraits, CellFieldTraits > InterpYF2C;  ///< Interpolate face to cell in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZDIR, ZSideFieldTraits, CellFieldTraits > InterpZF2C;  ///< Interpolate face to cell in z-dir

  // divergence operators - Side to Cell
  typedef SpatialOperator< LinAlgTrilinos, Divergence, XDIR, XSideFieldTraits, CellFieldTraits >  DivXF2C;     ///< Divergence of a face field in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Divergence, YDIR, YSideFieldTraits, CellFieldTraits >  DivYF2C;     ///< Divergence of a face field in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZDIR, ZSideFieldTraits, CellFieldTraits >  DivZF2C;     ///< Divergence of a face field in z-dir

  // divergence operators - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, Divergence, XDIR, CellFieldTraits, XSideFieldTraits >  DivXC2F;     ///< Divergence of a cell field in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Divergence, YDIR, CellFieldTraits, YSideFieldTraits >  DivYC2F;     ///< Divergence of a cell field in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZDIR, CellFieldTraits, ZSideFieldTraits >  DivZC2F;     ///< Divergence of a cell field in z-dir
									    	        
  // gradient operators - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, Gradient, XDIR, CellFieldTraits, XSideFieldTraits >    GradXC2F;     ///< Gradient of a cell field in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Gradient, YDIR, CellFieldTraits, YSideFieldTraits >    GradYC2F;     ///< Gradient of a cell field in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZDIR, CellFieldTraits, ZSideFieldTraits >    GradZC2F;     ///< Gradient of a cell field in z-dir

  // gradient operators - Side to Cell
  typedef SpatialOperator< LinAlgTrilinos, Gradient, XDIR, XSideFieldTraits, CellFieldTraits >    GradXF2C;     ///< Gradient of a cell field in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Gradient, YDIR, YSideFieldTraits, CellFieldTraits >    GradYF2C;     ///< Gradient of a cell field in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZDIR, ZSideFieldTraits, CellFieldTraits >    GradZF2C;     ///< Gradient of a cell field in z-dir


  // linear interpolant - Side to Edge
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XDIR, YSideFieldTraits, ZEdgeYDirFieldTraits >  InterpX_YF2ZE;  ///< Interpolate y-face field in x-dir to a edge z-edge
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XDIR, ZSideFieldTraits, YEdgeZDirFieldTraits >  InterpX_ZF2YE;  ///< Interpolate z-face field in x-dir to a edge y-edge
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YDIR, XSideFieldTraits, ZEdgeXDirFieldTraits >  InterpY_XF2ZE;  ///< Interpolate x-face field in y-dir to a edge z-edge
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YDIR, ZSideFieldTraits, XEdgeZDirFieldTraits >  InterpY_ZF2XE;  ///< Interpolate y-face field in y-dir to a edge x-edge
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZDIR, XSideFieldTraits, YEdgeXDirFieldTraits >  InterpZ_XF2YE;  ///< Interpolate z-face field in z-dir to a edge y-edge
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZDIR, YSideFieldTraits, XEdgeYDirFieldTraits >  InterpZ_YF2XE;  ///< Interpolate y-face field in z-dir to a edge x-edge

  // gradient operators - Side to Edge
  typedef SpatialOperator< LinAlgTrilinos, Gradient, XDIR, YSideFieldTraits, ZEdgeYDirFieldTraits >  GradX_YF2ZE;  ///< Gradient in x-dir y-face field to a edge z-edge
  typedef SpatialOperator< LinAlgTrilinos, Gradient, XDIR, ZSideFieldTraits, YEdgeZDirFieldTraits >  GradX_ZF2YE;  ///< Gradient in x-dir z-face field to a edge y-edge
  typedef SpatialOperator< LinAlgTrilinos, Gradient, YDIR, XSideFieldTraits, ZEdgeXDirFieldTraits >  GradY_XF2ZE;  ///< Gradient in y-dir x-face field to a edge z-edge
  typedef SpatialOperator< LinAlgTrilinos, Gradient, YDIR, ZSideFieldTraits, XEdgeZDirFieldTraits >  GradY_ZF2XE;  ///< Gradient in y-dir y-face field to a edge x-edge
  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZDIR, XSideFieldTraits, YEdgeXDirFieldTraits >  GradZ_XF2YE;  ///< Gradient in z-dir z-face field to a edge y-edge
  typedef SpatialOperator< LinAlgTrilinos, Gradient, ZDIR, YSideFieldTraits, XEdgeYDirFieldTraits >  GradZ_YF2XE;  ///< Gradient in z-dir y-face field to a edge x-edge

  // divergence operators - Edge to Side
  typedef SpatialOperator< LinAlgTrilinos, Divergence, YDIR, ZEdgeYDirFieldTraits, XSideFieldTraits > DivY_ZE2XF;  ///< Divergence in y-dir z-edge field to x-face field
  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZDIR, YEdgeZDirFieldTraits, XSideFieldTraits > DivZ_YE2XF;  ///< Divergence in z-dir y-edge field to x-face field
  typedef SpatialOperator< LinAlgTrilinos, Divergence, XDIR, ZEdgeXDirFieldTraits, YSideFieldTraits > DivX_ZE2YF;  ///< Divergence in x-dir z-edge field to y-face field
  typedef SpatialOperator< LinAlgTrilinos, Divergence, ZDIR, XEdgeZDirFieldTraits, YSideFieldTraits > DivZ_XE2YF;  ///< Divergence in z-dir x-edge field to y-face field
  typedef SpatialOperator< LinAlgTrilinos, Divergence, XDIR, YEdgeXDirFieldTraits, ZSideFieldTraits > DivX_YE2ZF;  ///< Divergence in x-dir y-edge field to z-face field
  typedef SpatialOperator< LinAlgTrilinos, Divergence, YDIR, XEdgeYDirFieldTraits, ZSideFieldTraits > DivY_XE2ZF;  ///< Divergence in y-dir x-edge field to z-face field

  // Operator Types
  //==================================================================








  //==================================================================
  // SCRATCH OPERATORS

  typedef SpatialOperator< LinAlgTrilinos, Scratch<3>, XDIR, CellFieldTraits, CellFieldTraits  > SxCell;
  typedef SpatialOperator< LinAlgTrilinos, Scratch<3>, YDIR, CellFieldTraits, CellFieldTraits  > SyCell;
  typedef SpatialOperator< LinAlgTrilinos, Scratch<3>, ZDIR, CellFieldTraits, CellFieldTraits  > SzCell;

  typedef SpatialOperator< LinAlgTrilinos, Scratch<1>, XDIR, XSideFieldTraits, XSideFieldTraits > SxSide;
  typedef SpatialOperator< LinAlgTrilinos, Scratch<1>, YDIR, YSideFieldTraits, YSideFieldTraits > SySide;
  typedef SpatialOperator< LinAlgTrilinos, Scratch<1>, YDIR, ZSideFieldTraits, ZSideFieldTraits > SzSide;

  typedef SpatialOperator< LinAlgTrilinos, Scratch<2>, XDIR, CellFieldTraits, XSideFieldTraits > SxCellSide;
  typedef SpatialOperator< LinAlgTrilinos, Scratch<2>, YDIR, CellFieldTraits, YSideFieldTraits > SyCellSide;
  typedef SpatialOperator< LinAlgTrilinos, Scratch<2>, ZDIR, CellFieldTraits, ZSideFieldTraits > SzCellSide;

  // SCRATCH OPERATORS
  //==================================================================



} // namespace FVStaggeredUniform
} // namespace SpatialOps



namespace SpatialOps{


  template<typename Dir, typename Location, typename SrcGhost, typename DestGhost >
  struct OpAssemblerSelector< FVStaggeredUniform::Gradient, Dir, Location, SrcGhost, DestGhost >
  {
    typedef FVStaggeredUniform::GradientAssembler<Dir,Location,SrcGhost,DestGhost>  Assembler;
  };

  template< typename Dir, typename Location, typename SrcGhost, typename DestGhost >
  struct OpAssemblerSelector< FVStaggeredUniform::Divergence, Dir, Location, SrcGhost, DestGhost >
  {
    typedef FVStaggeredUniform::DivergenceAssembler<Dir,Location,SrcGhost,DestGhost>  Assembler;
  };

  template< typename Dir, typename Location, typename SrcGhost, typename DestGhost >
  struct OpAssemblerSelector< FVStaggeredUniform::Interpolant, Dir, Location, SrcGhost, DestGhost >
  {
    typedef FVStaggeredUniform::LinearInterpolantAssembler<Dir,Location,SrcGhost,DestGhost>  Assembler;
  };

  template< int N, typename Dir, typename Location, typename SrcGhost, typename DestGhost >
  struct OpAssemblerSelector< FVStaggeredUniform::Scratch<N>, Dir, Location, SrcGhost, DestGhost >
  {
    typedef FVStaggeredUniform::ScratchAssembler< FVStaggeredUniform::Scratch<N>::NumNonZero,Dir,Location,SrcGhost,DestGhost>  Assembler;
  };


} // namespace SpatialOps

#endif // UT_FV2ndOrderTypes_h
