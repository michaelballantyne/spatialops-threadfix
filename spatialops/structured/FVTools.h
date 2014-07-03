/*
 * Copyright (c) 2014 The University of Utah
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
   *  \addtogroup fields
   *  @{
   *
   *  \file FVTools.h
   *
   *
   *  \brief Provides function templates useful for field and operator
   *         creation for structured meshes.  Specialize these
   *         templates to your needs.
   */

  //------------------------------------------------------------------

  /**
   *  \fn int get_dim_with_ghost( const int, const int, const int, const int )
   *
   *  \brief obtain the number of points in the x direction
   *
   *  \param nNoGhost number of points in the current direction excluding
   *    ghost cells
   *
   *  \param minusGhost the number of ghost cells on the negative face
   *
   *  \param plusGhost the number of ghost cells on the positive face
   *
   *  \param bc the number of boundary cells on the positive face
   *
   *  \return the number of points in the current direction, including ghost cells
   *    and boundary cells
   *
   */
  inline int get_dim_with_ghost( const int nNoGhost,
                                 const int minusGhost,
                                 const int plusGhost,
                                 const int bc )
  {
    return ( nNoGhost > 1
             ? ( nNoGhost + minusGhost + plusGhost + bc )
             : 1 );
  }

  //------------------------------------------------------------------

  /**
   *  \fn MemoryWindow get_window_with_ghost( const IntVec&, const GhostData&, const BoundaryCellInfo& )
   *  \brief Obtain the memory window for a field on a patch that is a single, contiguous memory block
   *
   *  \param localDim number of points in each direction excluding
   *    ghost cells
   *
   *  \param ghost the GhostData information
   *
   *  \param bc BoundaryCellInfo describing the behavior of a field when a (+) side
   *   boundary is present.  Note that a MemoryWindow obtained here is paired for
   *   use specifically with fields that share common BoundaryCellInfo.
   *
   *  \return a MemoryWindow consistent with the given information
   */
  MemoryWindow
  inline get_window_with_ghost( const IntVec& localDim,
                                const GhostData& ghost,
                                const BoundaryCellInfo& bc )
  {
      return MemoryWindow( IntVec( get_dim_with_ghost( localDim[0], ghost.get_minus(0), ghost.get_plus(0), bc.has_extra(0) ),
                                   get_dim_with_ghost( localDim[1], ghost.get_minus(1), ghost.get_plus(1), bc.has_extra(1) ),
                                   get_dim_with_ghost( localDim[2], ghost.get_minus(2), ghost.get_plus(2), bc.has_extra(2) ) ) );
  }

  //------------------------------------------------------------------

  /**
   *  @}
   */

}// namespace structured
}// namespace SpatialOps

#endif // FVToolsTemplates_h

