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

namespace SpatialOps{
namespace structured{


  /**
   * \class GhostDataRT
   * \date July, 2013
   * \author James C. Sutherland
   * \brief Holds information about the number of ghost cells on each side of the domain
   */
  class GhostDataRT
  {
    IntVec minus_, plus_;

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
    GhostDataRT( const int nx, const int px,
                 const int ny, const int py,
                 const int nz, const int pz );

    /**
     * @brief Construct a GhostData
     * @param minus Number of ghost cells on the (-) x, y, and z faces
     * @param plus  Number of ghost cells on the (+) x, y, and z faces
     */
    GhostDataRT( const IntVec& minus,
                 const IntVec& plus );

    /**
     * \brief construct a GhostData with the same number of ghost cells on each face
     * @param n the number of ghost cells on each face (defaults to zero)
     * @param bcx (defaults false) true if there is a physical boundary on the (+x) face
     * @param bcy (defaults false) true if there is a physical boundary on the (+y) face
     * @param bcz (defaults false) true if there is a physical boundary on the (+z) face
     */
    GhostDataRT( const int n=0 );

    GhostDataRT( const GhostDataRT& );
    GhostDataRT& operator=( const GhostDataRT& );

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
    inline void set_minus( const IntVec& minus ){ minus_ = minus; }

    /**
     * @brief set the number of ghost cells on the requested (+) face (0=x, 1=y, 2=z)
     */
    inline void set_plus( const IntVec& plus ){ plus_ = plus; }

    GhostDataRT  operator+ ( const GhostDataRT& ) const;
    GhostDataRT& operator+=( const GhostDataRT& );
    GhostDataRT  operator- ( const GhostDataRT& ) const;
    GhostDataRT& operator-=( const GhostDataRT& );
  };

  std::ostream& operator<<( std::ostream&, const GhostDataRT& );










    /**
     *  \struct GhostData
     *  \brief Used for specifying field type traits.
     *
     *  \par Key Definitions
     *  \li \c nX The negative x-component of the GhostData
     *  \li \c pX The positive x-component of the GhostData
     *  \li \c bX The positive (boundary condition) x-component of the GhostData
     *  \li \c nY The negative y-component of the GhostData
     *  \li \c pY The positive y-component of the GhostData
     *  \li \c bY The positive (boundary condition) y-component of the GhostData
     *  \li \c nZ The negative z-component of the GhostData
     *  \li \c pZ The positive z-component of the GhostData
     *  \li \c bZ The positive (boundary condition) z-component of the GhostData
     *
     *  Examples:
     *  \code
     *  typedef GhostData<1,2,3,4,5,6,7,8,9> T;
     *  assert( T::nX == 1 );
     *  assert( T::pX == 2 );
     *  assert( T::bX == 3 );
     *  assert( T::nY == 4 );
     *  assert( T::pY == 5 );
     *  assert( T::bY == 6 );
     *  assert( T::nZ == 7 );
     *  assert( T::pZ == 8 );
     *  assert( T::bZ == 9 );
     *  \endcode
     */
      template<int nx, int px, int bx,
               int ny, int py, int by,
               int nz, int pz, int bz>
    struct GhostData
    {
      enum Component{
          nX=nx,  ///< The negative x component of the GhostData
          pX=px,  ///< The positive x component of the GhostData
          bX=bx,  ///< The positive (boundary condition) x component of the GhostData
          nY=ny,  ///< The negative y component of the GhostData
          pY=py,  ///< The positive y component of the GhostData
          bY=by,  ///< The positive (boundary condition) y component of the GhostData
          nZ=nz,  ///< The negative z component of the GhostData
          pZ=pz,  ///< The positive z component of the GhostData
          bZ=bz   ///< The positive (boundary condition) z component of the GhostData
      };

      /**
       * \brief Writes the GhostData to a string.
       * \return a string value representing the GhostData.
       */
      static inline std::string print() {
        std::stringstream s;
        s << "("
          << "<" << std::setw(2) << nX << "," << std::setw(2) << pX << "," << std::setw(2) << bX << ">,"
          << "<" << std::setw(2) << nY << "," << std::setw(2) << pY << "," << std::setw(2) << bY << ">,"
          << "<" << std::setw(2) << nZ << "," << std::setw(2) << pZ << "," << std::setw(2) << bZ << ">)";
        return s.str();
      }

      static inline IntVec neg_int_vec(){
        return IntVec( nX, nY, nZ );
      }

      static inline IntVec pos_int_vec(){
        return IntVec( pX, pY, pZ );
      }

      static inline IntVec bc_pos_int_vec(){
        return IntVec( bX, bY, bZ );
      }

      static inline IntVec real_pos_int_vec(IntVec const bc){
        return IntVec( (bc[0] ? bX : pX),
                       (bc[1] ? bY : pY),
                       (bc[2] ? bZ : pZ) );
      }

      //returns the change to offset caused by ghost cells
      static inline IntVec offset() {
        return neg_int_vec();
      }

      //returns the change to extent caused by ghost cells (and boundary conditions)
      static inline IntVec extent(IntVec const bc) {
        return neg_int_vec() + real_pos_int_vec(bc);
      }


    };

    struct InfiniteGhostData;

    //------------------------------------------------------------------

    template<typename T1, typename T2>
    struct Invalidate;

    /**
     *  \struct Invalidate
     *  \brief Perform compile-time reduction/invalidation of ghost cells in an GhostData type from an GhostData type
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef GhostData< 1, 0, 1, 2, 3, 6, 4, 5, 6 > T1;
     *    typedef GhostData< 1, 2, 3, 4, 5, 6, 7, 8, 9 > T2;
     *    typedef GhostData< 0, 2, 2, 2, 2, 0, 3, 3, 3 > T3;
     *    typedef Invalidate< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int nx1, int px1, int bx1, int nx2, int px2, int bx2,
             int ny1, int py1, int by1, int ny2, int py2, int by2,
             int nz1, int pz1, int bz1, int nz2, int pz2, int bz2>
    struct Invalidate<GhostData<nx1,px1,bx1,ny1,py1,by1,nz1,pz1,bz1>,
                      GhostData<nx2,px2,bx2,ny2,py2,by2,nz2,pz2,bz2> > {
        GhostData<nx1 - nx2, px1 - px2, bx1 - bx2,
                  ny1 - ny2, py1 - py2, by1 - by2,
                  nz1 - nz2, pz1 - pz2, bz1 - bz2> typedef result;
    };

    /**
     *  \struct Invalidate
     *  \brief Perform compile-time reduction/invalidation of ghost cells in an GhostData type from an IndexTriplet type
     *
     *  A negative value in the IndexTriplet invalidates the appropriate negative direction in the GhostData.
     *
     *  A positive value in the IndexTriplet invalidates the appropriate positive direction in the GhostData.
     *
     *  A positive value in the IndexTriplet also invalidates the appropriate positive (boundary condition) direction in the GhostData.
     *
     *  The magnitude of each invalidation is the absolute value of the values in the IndexTriplet.
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef GhostData< 0, 1, 1, 1, 1, 1, 0, -1, 0 > T1;
     *    typedef GhostData< 1, 1, 1, 2, 1, 3, 0, 0, 2 > T2;
     *    typedef IndexTriplet< -1, 1, 2 > T3;
     *    typedef Invalidate< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int x, int nx, int px, int bx,
             int y, int ny, int py, int by,
             int z, int nz, int pz, int bz>
    struct Invalidate<GhostData<nx,px,bx,ny,py,by,nz,pz,bz>,
                      IndexTriplet<x,y,z> > {
        GhostData<(nx + (x < 0 ? x : 0)),
                  (px - (x > 0 ? x : 0)),
                  (bx - (x > 0 ? x : 0)),
                  (ny + (y < 0 ? y : 0)),
                  (py - (y > 0 ? y : 0)),
                  (by - (y > 0 ? y : 0)),
                  (nz + (z < 0 ? z : 0)),
                  (pz - (z > 0 ? z : 0)),
                  (bz - (z > 0 ? z : 0))> typedef result;
    };

    /**
     *  \struct Invalidate
     *  \brief Perform compile-time reduction/invalidation of ghost cells in an InfiniteGhostData type from any type
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult1, MyResult2, and T1 are all the same type:
     *    typedef InfiniteGhostData T1;
     *    typedef GhostData< 1, 1, 1, 2, 1, 3, 0, 2, 1 > T2;
     *    typedef IndexTriplet< -1, 1, 2 > T3;
     *    typedef Invalidate< T1, T2 >::result  MyResult1;
     *    typedef Invalidate< T1, T3 >::result  MyResult2;
     *   \endcode
     */
    template<typename Other>
    struct Invalidate<InfiniteGhostData,
                      Other> {
        InfiniteGhostData typedef result;
    };

    /**
     *  \struct GhostFromField
     *  \brief Perform compile-time computation to pull GhostData from a FieldType
     */
    template<typename FieldType>
    struct GhostFromField {
        typename FieldType::Ghost::NGhostMinus typedef GM;
        typename FieldType::Ghost::NGhostPlus  typedef GP;
        typename FieldType::Location::BCExtra  typedef BC;

        GhostData<GM::X,
                  GP::X,
                  (GP::X + BC::X),
                  GM::Y,
                  GP::Y,
                  (GP::Y + BC::Y),
                  GM::Z,
                  GP::Z,
                  (GP::Z + BC::Z)> typedef result;
    };

    /**
     *  \struct MinimumGhostFromField
     *  \brief Perform compile-time computation to pull Minimum GhostData from a FieldType
     *
     *  MinimumGhostFromField is the inverse of GhostFromField.
     *  MinimumGhostFromField returns the minimum number of ghost cells needed to maintain the interior cells of the field.
     *  Practically speaking, the minimum is just the boundary condition cells.
     */
    template<typename FieldType>
    struct MinimumGhostFromField {
        typename FieldType::Location::BCExtra typedef BC;

        GhostData<0,
                  0,
                  BC::X,
                  0,
                  0,
                  BC::Y,
                  0,
                  0,
                  BC::Z> typedef result;
    };

    /**
     *  \struct TestField
     *  \brief Provides a way to test GhostFromField without a full SpatialField
     */
    template<typename G, typename BC>
    struct TestField {
        G typedef Ghost;
        struct Location { BC typedef BCExtra; };
    };

    template<typename T1, typename T2>
    struct Minimum;

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two GhostData types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef GhostData< 1, 1, 1, 1, 1, 1, 0, 1, 1> T1;
     *    typedef GhostData< 1, 2, 1, 2, 1, 2, 0, 1, 2> T2;
     *    typedef GhostData< 2, 1, 2, 1, 2, 1, 1, 1, 1> T3;
     *    typedef Minimum< T2, T3 >::result  MyResult;
     *   \endcode
     */
    template<int nx1, int px1, int bx1, int nx2, int px2, int bx2,
             int ny1, int py1, int by1, int ny2, int py2, int by2,
             int nz1, int pz1, int bz1, int nz2, int pz2, int bz2>
    struct Minimum<GhostData<nx1,px1,bx1,ny1,py1,by1,nz1,pz1,bz1>,
                   GhostData<nx2,px2,bx2,ny2,py2,by2,nz2,pz2,bz2> > {
        GhostData<Min<nx1,nx2>::result,
                  Min<px1,px2>::result,
                  Min<bx1,bx2>::result,
                  Min<ny1,ny2>::result,
                  Min<py1,py2>::result,
                  Min<by1,by2>::result,
                  Min<nz1,nz2>::result,
                  Min<pz1,pz2>::result,
                  Min<bz1,bz2>::result> typedef result;
    };

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two GhostData types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef GhostData< 1, 2, 1, 2, 1, 2, 2, 0, 1> T1;
     *    typedef InfiniteGhostData T2;
     *    typedef Minimum< T1, T2 >::result  MyResult;
     *   \endcode
     */
    template<int nx, int px, int bx,
             int ny, int py, int by,
             int nz, int pz, int bz>
    struct Minimum<GhostData<nx,px,bx,ny,py,by,nz,pz,bz>,
                   InfiniteGhostData> {
        GhostData<nx,px,bx,ny,py,by,nz,pz,bz> typedef result;
    };

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two GhostData types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef GhostData< 1, 2, 1, 2, 1, 2, 2, 0, 1> T1;
     *    typedef InfiniteGhostData T2;
     *    typedef Minimum< T2, T1 >::result  MyResult;
     *   \endcode
     */
    template<int nx, int px, int bx,
             int ny, int py, int by,
             int nz, int pz, int bz>
    struct Minimum<InfiniteGhostData,
                   GhostData<nx,px,bx,ny,py,by,nz,pz,bz> > {
        GhostData<nx,px,bx,ny,py,by,nz,pz,bz> typedef result;
    };

    /**
     *  \struct Minimum
     *  \brief Perform compile-time minimum of two GhostData types
     *
     *  Example usage:
     *   \code
     *    // In the following, MyResult and T1 are the same type:
     *    typedef InfiniteGhostData T1;
     *    typedef Minimum< T1, T1 >::result  MyResult;
     *   \endcode
     */
    template<>
    struct Minimum<InfiniteGhostData,
                   InfiniteGhostData> {
        InfiniteGhostData typedef result;
    };
  } // namespace structured
} // namespace SpatialOps

/**
 * @}
 */

#endif /* SpatialOps_GhostData_h */
