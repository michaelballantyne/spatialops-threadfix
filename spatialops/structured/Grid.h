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

#ifndef SpatialOps_structured_Grid_h
#define SpatialOps_structured_Grid_h

#include <vector>
#include <spatialops/structured/MemoryWindow.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \class Grid
   *  \author James C. Sutherland
   *  \brief Provides a simple interface to set coordinate fields.
   */
  class Grid
  {
    const IntVec npts_;
    const std::vector<double> length_;
    std::vector<double> spacing_;
  public:
    /**
     *  \param npts the number of scalar volume cells in each direction
     *  \param length the domain length in each direction
     */
    Grid( const IntVec npts,
          const std::vector<double>& length );

    /**
     *  \brief obtain the grid spacing in the requested direction
     */
    template< typename CoordT >
    double spacing() const;

    inline const IntVec& extent() const{ return npts_; }

    inline int extent( const size_t i ) const{ return npts_[i]; }

    inline const std::vector<double>& length() const{ return length_; }
    inline double length( const size_t i ) const{ return length_[i]; }

    /**
     *  \brief set the coordinates on the given field.
     *
     *  Examples:
     *  \code
     *  SVolField xsvol, ysvol;
     *  grid.set_coord<XDIR>( xsvol );
     *  grid.set_coord<YDIR>( ysvol );
     *
     *  XSurfYField xxsy, zxsy;
     *  grid.set_coord<XDIR>( xxsy );
     *  grid.set_coord<ZDIR>( zxsy );
     *  \endcode
     */
    template< typename CoordT, typename FieldT >
    void set_coord( FieldT& f ) const;
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_Grid_h
