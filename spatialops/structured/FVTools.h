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
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \return the number of points in the x-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NX for.
   */
  template<typename FieldT>
  int get_nx_with_ghost( const int nxNoGhost, const bool hasPlusFaceX )
  {
    return ( nxNoGhost>1
             ? ( nxNoGhost + FieldT::Ghost::NGhostMinus::x_value() + FieldT::Ghost::NGhostPlus::x_value()
                 + (hasPlusFaceX ? FieldT::Location::BCExtra::X : 0) )
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
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \return the number of points in the y-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NY for.
   */
  template<typename FieldT>
  int get_ny_with_ghost( const int nyNoGhost, const bool hasPlusFaceY )
  {
    return ( nyNoGhost>1
             ? ( nyNoGhost + FieldT::Ghost::NGhostMinus::y_value() + FieldT::Ghost::NGhostPlus::y_value()
                 + ( hasPlusFaceY ? FieldT::Location::BCExtra::Y : 0 ) )
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
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \return the number of points in the z-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NZ for.
   */
  template<typename FieldT>
  inline int get_nz_with_ghost( const int nzNoGhost, const bool hasPlusFaceZ )
  {
    return ( nzNoGhost>1
             ? ( nzNoGhost + FieldT::Ghost::NGhostMinus::z_value() + FieldT::Ghost::NGhostPlus::z_value()
                 + ( hasPlusFaceZ ? FieldT::Location::BCExtra::Z : 0 ) )
             : 1 );
  }

  template<> inline int get_nx_with_ghost<double>( const int nxNoGhost, const bool hasPlusFaceX ){ return nxNoGhost; }
  template<> inline int get_ny_with_ghost<double>( const int nyNoGhost, const bool hasPlusFaceY ){ return nyNoGhost; }
  template<> inline int get_nz_with_ghost<double>( const int nzNoGhost, const bool hasPlusFaceZ ){ return nzNoGhost; }

  /**
   *  \fn int get_dim_with_ghost( const IntVec&, const bool, const bool, const bool )
   *
   *  \brief obtain the number of points in each direction for the given field type
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \return the number of points in each direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
   */
  template<typename FieldT>
  inline IntVec get_dim_with_ghost( const IntVec& dimNoGhost,
                             const bool hasPlusFaceX,
                             const bool hasPlusFaceY,
                             const bool hasPlusFaceZ )
  {
    return IntVec( get_nx_with_ghost<FieldT>( dimNoGhost[0], hasPlusFaceX ),
                   get_ny_with_ghost<FieldT>( dimNoGhost[1], hasPlusFaceY ),
                   get_nz_with_ghost<FieldT>( dimNoGhost[2], hasPlusFaceZ ) );
  }


  /**
   *  \fn MemoryWindow get_window_with_ghost( const IntVec&, const bool, const bool, const bool )
   *  \brief Obtain the memory window for a field on a patch that is a single, contiguous memory block
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \return the total number of points in the field, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
   */
  template<typename FieldT>
  MemoryWindow
  inline get_window_with_ghost( const IntVec& localDim, const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    const IntVec dim( get_nx_with_ghost<FieldT>(localDim[0],hasPlusFaceX),
                      get_ny_with_ghost<FieldT>(localDim[1],hasPlusFaceY),
                      get_nz_with_ghost<FieldT>(localDim[2],hasPlusFaceZ) );
    return MemoryWindow( dim, hasPlusFaceX, hasPlusFaceY, hasPlusFaceZ );
  }

  //------------------------------------------------------------------

  /**
   *  \brief Obtain the memory window for a field on a patch that is a subset of a larger memory block
   *
   *  \param globalDim the global dimensionality of the memory block (without ghost cells)
   *  \param localDim  the dimensionality of the field in consideration (without ghost cells)
   *  \param offset    the offset (start ijk index) of the local field in the global address space.
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   */
  template<typename FieldT>
  MemoryWindow
  inline get_window_with_ghost( const IntVec& globalDim, const IntVec& localDim, const IntVec& offset,
                         const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    const IntVec dimLoc( get_nx_with_ghost<FieldT>(localDim[0],hasPlusFaceX),
                         get_ny_with_ghost<FieldT>(localDim[1],hasPlusFaceY),
                         get_nz_with_ghost<FieldT>(localDim[2],hasPlusFaceZ) );
    const IntVec dimGlob( get_nx_with_ghost<FieldT>(globalDim[0],hasPlusFaceX),
                          get_ny_with_ghost<FieldT>(globalDim[1],hasPlusFaceY),
                          get_nz_with_ghost<FieldT>(globalDim[2],hasPlusFaceZ) );
    return MemoryWindow( dimGlob, offset, dimLoc, hasPlusFaceX, hasPlusFaceY, hasPlusFaceZ );
  }

  //==================================================================

  /**
   *  @}
   *  @}
   */

}// namespace structured
}// namespace SpatialOps

#endif // FVToolsTemplates_h

