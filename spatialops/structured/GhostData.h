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

#ifndef SpatialOps_GhostData_h
#define SpatialOps_GhostData_h

/**
 *  \file   GhostData.h
 *
 *  \date   August, 2012
 *  \author Christopher Earl
 *
 *  \addtogroup structured
 *  @{
 */

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsTools.h>
#include <spatialops/SpatialOpsDefs.h>

#include <spatialops/structured/IntVec.h>
#include <spatialops/structured/IndexTriplet.h>

#include <iomanip>
#include <string>
#include <sstream>

#define GHOST_MAX 9001

namespace SpatialOps{
namespace structured{


  /**
   * \class GhostData
   * \date July, 2013
   * \author James C. Sutherland
   * \brief Holds information about the number of ghost cells on each side of the domain
   */
  class GhostData
  {
    IntVec minus_, plus_;
    bool isInf_;

  public:

    /**
     * @brief Construct a GhostData
     * @param nx Number of ghost cells on the -x face
     * @param px Number of ghost cells on the +x face
     * @param ny Number of ghost cells on the -y face
     * @param py Number of ghost cells on the +y face
     * @param nz Number of ghost cells on the -z face
     * @param pz Number of ghost cells on the +z face
     */
    GhostData( const int nx, const int px,
                 const int ny, const int py,
                 const int nz, const int pz );

    /**
     * @brief Construct a GhostData
     * @param minus Number of ghost cells on the (-) x, y, and z faces
     * @param plus  Number of ghost cells on the (+) x, y, and z faces
     */
    GhostData( const IntVec& minus,
                 const IntVec& plus );

    /**
     * \brief construct a GhostData with the same number of ghost cells on each face
     * @param n the number of ghost cells on each face (defaults to zero)
     */
    GhostData( const int n=0 );

    GhostData( const GhostData& );
    GhostData& operator=( const GhostData& );

    /**
     * @brief obtain the IntVec containing the number of ghost cells on the (-) faces
     */
    inline IntVec get_minus() const{ return minus_; }

    /**
     * @brief obtain the number of ghost cells on the requested (-) face (0=x, 1=y, 2=z)
     */
    inline int get_minus( const int i ) const{ return minus_[i]; }

    /**
     * @brief obtain the IntVec containing the number of ghost cells on the (+) faces
     */
    inline IntVec get_plus() const{ return plus_; }

    /**
     * @brief obtain the number of ghost cells on the requested (+) face (0=x, 1=y, 2=z)
     */
    inline int get_plus( const int i ) const{ return plus_[i]; }

    /**
     * @brief set the number of ghost cells on the requested (-) face (0=x, 1=y, 2=z)
     */
    void set_minus( const IntVec& );

    /**
     * @brief set the number of ghost cells on the requested (+) face (0=x, 1=y, 2=z)
     */
    void set_plus( const IntVec& );

    GhostData  operator+ ( const GhostData& ) const;
    GhostData& operator+=( const GhostData& );
    GhostData  operator- ( const GhostData& ) const;
    GhostData& operator-=( const GhostData& );

    bool operator==( const GhostData& ) const;
  };

  GhostData min( const GhostData&, const GhostData& );

  GhostData point_to_ghost( const IntVec& );

  std::ostream& operator<<( std::ostream&, const GhostData& );


  } // namespace structured
} // namespace SpatialOps

/**
 * @}
 */

#endif /* SpatialOps_GhostData_h */
