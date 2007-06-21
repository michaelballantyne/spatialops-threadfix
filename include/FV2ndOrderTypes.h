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
  // Operator Types

  // linear interpolants - Cell to Side
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, XDIR, CellFieldTraits, XSideFieldTraits > InterpXC2F;  ///< Interpolate cell to face in x-dir
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, YDIR, CellFieldTraits, YSideFieldTraits > InterpYC2F;  ///< Interpolate cell to face in y-dir
  typedef SpatialOperator< LinAlgTrilinos, Interpolant, ZDIR, CellFieldTraits, ZSideFieldTraits > InterpZC2F;  ///< Interpolate cell to face in z-dir
  //@}

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
