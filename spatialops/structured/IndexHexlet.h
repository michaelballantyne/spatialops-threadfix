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

#ifndef SpatialOps_IndexHexlet_h
#define SpatialOps_IndexHexlet_h

/**
 *  \file   IndexHexlet.h
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

#include <boost/static_assert.hpp>

namespace SpatialOps{
  namespace structured{

    /**
     *  \struct IndexHexlet
     *  \brief Used for specifying field type traits.
     *
     *  \par Key Definitions
     *  \li \c nX The negative x-component of the IndexHexlet
     *  \li \c pX The positive x-component of the IndexHexlet
     *  \li \c nY The negative y-component of the IndexHexlet
     *  \li \c pY The positive y-component of the IndexHexlet
     *  \li \c nZ The negative z-component of the IndexHexlet
     *  \li \c pZ The positive z-component of the IndexHexlet
     *  \li \c Abs The absolute value of this IndexHexlet
     *
     *  Examples:
     *  \code
     *  typedef IndexHexlet<1,2,3,4,5,6> T;
     *  assert( T::nX == 1 );
     *  assert( T::pX == 2 );
     *  assert( T::nY == 3 );
     *  assert( T::pY == 4 );
     *  assert( T::nZ == 5 );
     *  assert( T::pZ == 6 );
     *  \endcode
     */
    template<int nx, int px, int ny, int py, int nz, int pz>
    struct IndexHexlet
    {
      enum Component{
          nX=nx,  ///< The negative x component of the IndexHexlet
          pX=px,  ///< The positive x component of the IndexHexlet
          nY=ny,  ///< The negative y component of the IndexHexlet
          pY=py,  ///< The positive y component of the IndexHexlet
          nZ=nz,  ///< The negative z component of the IndexHexlet
          pZ=pz,  ///< The positive z component of the IndexHexlet
      };

      /**
       * \brief Writes the IndexHexlet to a string.
       * \return a string value representing the IndexHexlet.
       */
      static inline std::string print() {
        std::stringstream s;
        s << "("
          << "<" << std::setw(2) << nX << "," << std::setw(2) << pX << ">,"
          << "<" << std::setw(2) << nY << "," << std::setw(2) << pY << ">,"
          << "<" << std::setw(2) << nZ << "," << std::setw(2) << pZ << ">)";
        return s.str();
      }

      static inline IntVec neg_int_vec(){
        return IntVec( nX, nY, nZ );
      }

      static inline IntVec pos_int_vec(){
        return IntVec( pX, pY, pZ );
      }

    };

    struct InfiniteIndexHexlet;

    //------------------------------------------------------------------

    template<typename T1, typename T2>
    struct Invalidate;

    /**
     *  \struct Invalidate
     *  \brief Perform compile-time reduction/invalidation of ghost cells in an IndexHexlet type from an IndexHexlet type
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef IndexHexlet< 1, 0, 1, 2, 3, 6 > T1;
     *    typedef IndexHexlet< 1, 2, 3, 4, 5, 6 > T2;
     *    typedef IndexHexlet< 0, 2, 2, 2, 2, 0 > T3;
     *    typedef Invalidate< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int nx1, int px1, int nx2, int px2,
             int ny1, int py1, int ny2, int py2,
             int nz1, int pz1, int nz2, int pz2>
    struct Invalidate<IndexHexlet<nx1,px1,ny1,py1,nz1,pz1>,
                      IndexHexlet<nx2,px2,ny2,py2,nz2,pz2> > {
        IndexHexlet<nx1 - nx2, px1 - px2,
                    ny1 - ny2, py1 - py2,
                    nz1 - nz2, pz1 - pz2> typedef result;
    };

    /**
     *  \struct Invalidate
     *  \brief Perform compile-time reduction/invalidation of ghost cells in an IndexHexlet type from an IndexTriplet type
     *
     *  A negative value in the IndexTriplet invalidates the appropriate negative direction in the IndexHexlet.
     *
     *  A positive value in the IndexTriplet invalidates the appropriate positive direction in the IndexHexlet.
     *
     *  The magnitude of each invalidation is the absolute value of the values in the IndexTriplet.
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef IndexHexlet< 0, 1, 1, 1, 1, 1 > T1;
     *    typedef IndexHexlet< 1, 1, 1, 2, 1, 3 > T2;
     *    typedef IndexTriplet< -1, 1, 2 > T3;
     *    typedef Invalidate< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int x, int nx, int px,
             int y, int ny, int py,
             int z, int nz, int pz>
    struct Invalidate<IndexHexlet<nx,px,ny,py,nz,pz>,
                      IndexTriplet<x,y,z> > {
        IndexHexlet<(nx - (x < 0 ? Abs<x>::result : 0)),
                    (px - (x > 0 ? x : 0)),
                    (ny - (y < 0 ? Abs<y>::result : 0)),
                    (py - (y > 0 ? y : 0)),
                    (nz - (z < 0 ? Abs<z>::result : 0)),
                    (pz - (z > 0 ? z : 0))> typedef result;
    };

    /**
     *  \struct Invalidate
     *  \brief Perform compile-time reduction/invalidation of ghost cells in an InfiniteIndexHexlet type from any type
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult1, MyResult2, and T1 are all the same type:
     *    typedef InfiniteIndexHexlet T1;
     *    typedef IndexHexlet< 1, 1, 1, 2, 1, 3 > T2;
     *    typedef IndexTriplet< -1, 1, 2 > T3;
     *    typedef Invalidate< T1, T2 >::result  MyResult1;
     *    typedef Invalidate< T1, T3 >::result  MyResult2;
     *   \endcode
     */
    template<typename Other>
        struct Invalidate<InfiniteIndexHexlet,
                          Other> {
        InfiniteIndexHexlet typedef result;
    };

    /**
     *  \struct FromGhost
     *  \brief Perform compile-time conversion of a Ghost struct to an IndexHexlet type
     */
    template<typename Ghost>
    struct FromGhost {
        IndexHexlet<Ghost::NGhostMinus::X,
                    Ghost::NGhostPlus::X,
                    Ghost::NGhostMinus::Y,
                    Ghost::NGhostPlus::Y,
                    Ghost::NGhostMinus::Z,
                    Ghost::NGhostPlus::Z> typedef result;
    };

    template<typename T1, typename T2>
    struct Minimum;

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two IndexHexlet types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef IndexHexlet< 1, 1, 1, 1, 1, 1> T1;
     *    typedef IndexHexlet< 1, 2, 1, 2, 1, 2> T2;
     *    typedef IndexHexlet< 2, 1, 2, 1, 2, 1> T3;
     *    typedef Minimum< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int nx1, int px1, int nx2, int px2,
             int ny1, int py1, int ny2, int py2,
             int nz1, int pz1, int nz2, int pz2>
    struct Minimum<IndexHexlet<nx1,px1,ny1,py1,nz1,pz1>,
                   IndexHexlet<nx2,px2,ny2,py2,nz2,pz2> > {
        IndexHexlet<Min<nx1,nx2>::result,
                    Min<px1,px2>::result,
                    Min<ny1,ny2>::result,
                    Min<py1,py2>::result,
                    Min<nz1,nz2>::result,
                    Min<pz1,pz2>::result> typedef result;
    };

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two IndexHexlet types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef IndexHexlet< 1, 2, 1, 2, 1, 2> T1;
     *    typedef IndexHexlet< 1, 2, 1, 2, 1, 2> T2;
     *    typedef InfiniteIndexHexlet T3;
     *    typedef Minimum< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int nx, int px,
             int ny, int py,
             int nz, int pz>
    struct Minimum<IndexHexlet<nx,px,ny,py,nz,pz>,
                   InfiniteIndexHexlet> {
        IndexHexlet<nx,px,ny,py,nz,pz> typedef result;
    };

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two IndexHexlet types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef IndexHexlet< 1, 2, 1, 2, 1, 2> T1;
     *    typedef InfiniteIndexHexlet T2;
     *    typedef IndexHexlet< 1, 2, 1, 2, 1, 2> T3;
     *    typedef Minimum< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int nx, int px,
             int ny, int py,
             int nz, int pz>
    struct Minimum<InfiniteIndexHexlet,
                   IndexHexlet<nx,px,ny,py,nz,pz> > {
        IndexHexlet<nx,px,ny,py,nz,pz> typedef result;
    };

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two IndexHexlet types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef InfiniteIndexHexlet T1;
     *    typedef Minimum< T1, T1 >::result  MyResult;
     *   \endcode
     */
    template<>
    struct Minimum<InfiniteIndexHexlet,
                   InfiniteIndexHexlet> {
        InfiniteIndexHexlet typedef result;
    };
  } // namespace structured
} // namespace SpatialOps

/**
 * @}
 */

#endif /* SpatialOps_IndexHexlet_h */
