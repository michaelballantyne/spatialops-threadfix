/*
 * Copyright (c) 2011 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef FVToolsTemplates_h
#define FVToolsTemplates_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/SpatialOpsTools.h>

#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/GhostData.h>
#include <spatialops/structured/BoundaryCellInfo.h>

#include <vector>
#include <set>

namespace SpatialOps{
namespace structured{

  /**
   *  \file FVTools.h
   *
   *  \addtogroup structured
   *  @{
   *  \addtogroup fields
   *  @{
   *
   *  \brief Provides function templates useful for field and operator
   *         creation for structured meshes.  Specialize these
   *         templates to your needs.
   */

  //------------------------------------------------------------------

  /**
   *  \fn int get_nx_with_ghost( const int, const bool )
   *
   *  \brief obtain the number of points in the x direction
   *
   *  \param nxNoGhost number of points in the x-direction excluding
   *    ghost cells
   *
   *  \param ghost the GhostDataRT information
   *
   *  \return the number of points in the x-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NX for.
   */
  template<typename FieldT>
  inline int get_nx_with_ghost( const int nxNoGhost,
                                const GhostDataRT& ghost,
                                const bool hasBC )
  {
    return ( nxNoGhost>1
             ? ( nxNoGhost + ghost.get_minus(0) + ghost.get_plus(0)
                 + (hasBC ? FieldT::Location::BCExtra::x_value() : 0) )
             : 1 );
  }

  inline int get_nx_with_ghost( const int nxNoGhost,
                                const GhostDataRT& ghost,
                                const BoundaryCellInfo& bc )
  {
    return ( nxNoGhost>1
             ? ( nxNoGhost + ghost.get_minus(0) + ghost.get_plus(0)
                 + ( bc.has_bc(0) ? bc.num_extra(0) : 0) )
             : 1 );
  }

  /**
   *  \fn int get_ny_with_ghost( const int, const bool )
   *
   *  \brief obtain the number of points in the y direction
   *
   *  \param nyNoGhost number of points in the y-direction excluding
   *    ghost cells
   *
   *  \param ghost the GhostDataRT information
   *
   *  \return the number of points in the y-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NY for.
   */
  template<typename FieldT>
  inline int get_ny_with_ghost( const int nyNoGhost,
                                const GhostDataRT& ghost,
                                const bool hasBC )
  {
    return ( nyNoGhost>1
             ? ( nyNoGhost + ghost.get_minus(1) + ghost.get_plus(1)
                 + (hasBC ? FieldT::Location::BCExtra::y_value() : 0) )
             : 1 );
  }

  inline int get_ny_with_ghost( const int nyNoGhost,
                                const GhostDataRT& ghost,
                                const BoundaryCellInfo& bc )
  {
    return ( nyNoGhost>1
             ? ( nyNoGhost + ghost.get_minus(1) + ghost.get_plus(1)
                 + (bc.has_bc(1) ? bc.num_extra(1) : 0) )
             : 1 );
  }

  /**
   *  \fn int get_nz_with_ghost( const int, const bool )
   *
   *  \brief obtain the number of points in the z direction
   *
   *  \param nzNoGhost number of points in the z-direction excluding
   *    ghost cells
   *
   *  \param ghost the GhostDataRT information
   *
   *  \return the number of points in the z-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NZ for.
   */
  template<typename FieldT>
  inline int get_nz_with_ghost( const int nzNoGhost,
                                const GhostDataRT& ghost,
                                const bool hasBC )
  {
    return ( nzNoGhost>1
             ? ( nzNoGhost + ghost.get_minus(2) + ghost.get_plus(2)
                 + (hasBC ? FieldT::Location::BCExtra::z_value() : 0) )
             : 1 );
  }

  inline int get_nz_with_ghost( const int nzNoGhost,
                                const GhostDataRT& ghost,
                                const BoundaryCellInfo& bc )
  {
    return ( nzNoGhost>1
             ? ( nzNoGhost + ghost.get_minus(2) + ghost.get_plus(2)
                 + (bc.has_bc(2) ? bc.num_extra(2) : 0) )
             : 1 );
  }

  template<> inline int get_nx_with_ghost<double>( const int nxNoGhost, const GhostDataRT& ghost, const bool hasBC ){ return nxNoGhost; }
  template<> inline int get_ny_with_ghost<double>( const int nyNoGhost, const GhostDataRT& ghost, const bool hasBC ){ return nyNoGhost; }
  template<> inline int get_nz_with_ghost<double>( const int nzNoGhost, const GhostDataRT& ghost, const bool hasBC ){ return nzNoGhost; }

  /**
   *  \fn int get_dim_with_ghost( const IntVec&, const bool, const bool, const bool )
   *
   *  \brief obtain the number of points in each direction for the given field type
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param ghost the GhostDataRT information
   *
   *  \return the number of points in each direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
   */
  IntVec
  inline get_dim_with_ghost( const IntVec& dimNoGhost,
                             const GhostDataRT& ghost,
                             const BoundaryCellInfo& bc )
  {
    return IntVec( get_nx_with_ghost( dimNoGhost[0], ghost, bc ),
                   get_ny_with_ghost( dimNoGhost[1], ghost, bc ),
                   get_nz_with_ghost( dimNoGhost[2], ghost, bc ) );
  }


  /**
   *  \fn MemoryWindow get_window_with_ghost( const IntVec&, const GhostDataRT& )
   *  \brief Obtain the memory window for a field on a patch that is a single, contiguous memory block
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param ghost the GhostDataRT information
   *
   *  \param bc BoundaryCellInfo describing the behavior of a field when a (+) side
   *   boundary is present.  Note that a MemoryWindow obtained here is paired for
   *   use specifically with fields that share common BoundaryCellInfo.
   *
   *  \return the total number of points in the field, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
   */
  MemoryWindow
  inline get_window_with_ghost( const IntVec& localDim,
                                const GhostDataRT& ghost,
                                const BoundaryCellInfo& bc )
  {
    const IntVec dim( get_nx_with_ghost( localDim[0], ghost, bc ),
                      get_ny_with_ghost( localDim[1], ghost, bc ),
                      get_nz_with_ghost( localDim[2], ghost, bc ) );
    return MemoryWindow( dim );
  }

  //------------------------------------------------------------------

  /**
   *  \brief Obtain the memory window for a field on a patch that is a subset of a larger memory block
   *
   *  \param globalDim the global dimensionality of the memory block (without ghost cells)
   *  \param localDim  the dimensionality of the field in consideration (without ghost cells)
   *  \param offset    the offset (start ijk index) of the local field in the global address space.
   *
   *  \param ghost the GhostDataRT information
   */
  template<typename FieldT>
  MemoryWindow
  inline get_window_with_ghost( const IntVec& globalDim,
                                const IntVec& localDim,
                                const IntVec& offset,
                                const GhostDataRT& ghost )
  {
    const IntVec dimLoc( get_nx_with_ghost<FieldT>(localDim[0],ghost),
                         get_ny_with_ghost<FieldT>(localDim[1],ghost),
                         get_nz_with_ghost<FieldT>(localDim[2],ghost) );
    const IntVec dimGlob( get_nx_with_ghost<FieldT>(globalDim[0],ghost),
                          get_ny_with_ghost<FieldT>(globalDim[1],ghost),
                          get_nz_with_ghost<FieldT>(globalDim[2],ghost) );
    return MemoryWindow( dimGlob, offset, dimLoc );
  }

  //==================================================================

  /**
   *  @}
   *  @}
   */

}// namespace structured
}// namespace SpatialOps

#endif // FVToolsTemplates_h

