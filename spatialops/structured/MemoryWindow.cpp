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

#include <spatialops/structured/MemoryWindow.h>

#include <ostream>
using namespace std;

namespace SpatialOps{
  namespace structured{

    inline bool check_positive( const IntVec& v ){ return (v[0]>=0) & (v[1]>=0) & (v[2]>=0); }

#ifndef NDEBUG
    bool sanity_check( const IntVec& nglob, const IntVec& offset, const IntVec& extent ){
      return check_positive( nglob  ) &&
             check_positive( offset ) &&
             check_positive( extent ) &&
             extent[0] <= nglob[0]    &&
             extent[1] <= nglob[1]    &&
             extent[2] <= nglob[2];
    }
#endif

    MemoryWindow::MemoryWindow( const int npts[3],
                                const int offset[3],
                                const int extent[3],
                                const bool bcx,
                                const bool bcy,
                                const bool bcz )
    : nptsGlob_( npts ),
      offset_( offset ),
      extent_( extent ),
      bc_( bcx, bcy, bcz )
    {
#   ifndef NDEBUG
      assert( sanity_check( nptsGlob_, offset_, extent_ ) );
#   endif
    }

    MemoryWindow::MemoryWindow( const IntVec& npts,
                                const IntVec& offset,
                                const IntVec& extent,
                                const bool bcx,
                                const bool bcy,
                                const bool bcz )
    : nptsGlob_( npts ),
      offset_( offset ),
      extent_( extent ),
      bc_( bcx, bcy, bcz )
    {
#   ifndef NDEBUG
      assert( sanity_check( nptsGlob_, offset_, extent_ ) );
#   endif
    }

    MemoryWindow::MemoryWindow( const int npts[3],
                                const bool bcx,
                                const bool bcy,
                                const bool bcz )
    : nptsGlob_( npts ), offset_(0,0,0), extent_( npts ),
      bc_( bcx, bcy, bcz )
    {
#   ifndef NDEBUG
      assert( sanity_check( nptsGlob_, offset_, extent_ ) );
#   endif
    }

    MemoryWindow::MemoryWindow( const IntVec& npts,
                                const bool bcx,
                                const bool bcy,
                                const bool bcz )
    : nptsGlob_( npts ), offset_(0,0,0), extent_( npts ),
      bc_( bcx, bcy, bcz )
    {
#   ifndef NDEBUG
      assert( sanity_check( nptsGlob_, offset_, extent_ ) );
#   endif
    }

    MemoryWindow::MemoryWindow( const MemoryWindow& other )
    : nptsGlob_( other.nptsGlob_ ),
      offset_  ( other.offset_   ),
      extent_  ( other.extent_   ),
      bc_      ( other.bc_       )
    {}

    MemoryWindow::~MemoryWindow()
    {}

    std::vector<MemoryWindow>
    MemoryWindow::split( const IntVec splitPattern,
                         const IntVec nGhostMinus,
                         const IntVec nGhostPlus,
                         const IntVec bcExtents ) const
    {
      const IntVec extent = IntVec((extent_[0] == 1) ? 1 : extent_[0] - nGhostMinus[0] - nGhostPlus[0] - (bc_[0] ? bcExtents[0] : 0),
                                   (extent_[1] == 1) ? 1 : extent_[1] - nGhostMinus[1] - nGhostPlus[1] - (bc_[1] ? bcExtents[1] : 0),
                                   (extent_[2] == 1) ? 1 : extent_[2] - nGhostMinus[2] - nGhostPlus[2] - (bc_[2] ? bcExtents[2] : 0));
      const IntVec offset = IntVec((extent_[0] == 1) ? 0 : offset_[0] + nGhostMinus[0],
                                   (extent_[1] == 1) ? 0 : offset_[1] + nGhostMinus[1],
                                   (extent_[2] == 1) ? 0 : offset_[2] + nGhostMinus[2]);

#     ifndef NDEBUG
      for( size_t i=0; i<3; ++i ){
        assert( extent[i] >= splitPattern[i] );
        assert( extent[i] + offset[i] <= nptsGlob_[i] );
        assert( splitPattern[i] > 0 );
        assert( nGhostMinus[i] >= 0);
        assert( nGhostPlus[i] >= 0);
      }
#     endif

      // try to make windows close to same size
      const int nextra[3] = { extent[0] % splitPattern[0],
                              extent[1] % splitPattern[1],
                              extent[2] % splitPattern[2] };
      typedef std::vector<int>  Ints;
      Ints nxyz[3] = { Ints( splitPattern[0], extent[0] / splitPattern[0] ),
                       Ints( splitPattern[1], extent[1] / splitPattern[1] ),
                       Ints( splitPattern[2], extent[2] / splitPattern[2] ) };

      for( int idir=0; idir<3; ++idir ){
        for( int isp=0; isp<splitPattern[idir]; ++isp ){
          nxyz[idir][isp] += ( isp<nextra[idir] ? 1 : 0 );
        }
      }

      IntVec cumOffset = offset;  // keep track of how far each chunk is offset
      std::vector<MemoryWindow> children;
      for( int k=0; k<splitPattern[2]; ++k ){
        const bool bcz = (k==splitPattern[2]-1) ? bc_[2] : false;
        cumOffset[1] = offset[1];
        for( int j=0; j<splitPattern[1]; ++j ){
          const bool bcy = (j==splitPattern[1]-1) ? bc_[1] : false;
          cumOffset[0] = offset[0];
          for( int i=0; i<splitPattern[0]; ++i ){
            const bool bcx = (i==splitPattern[0]-1) ? bc_[0] : false;
            children.push_back( MemoryWindow(nptsGlob_,
                                             IntVec((extent_[0] == 1) ? 0 : cumOffset[0] - nGhostMinus[0],
                                                    (extent_[1] == 1) ? 0 : cumOffset[1] - nGhostMinus[1],
                                                    (extent_[2] == 1) ? 0 : cumOffset[2] - nGhostMinus[2]),
                                             IntVec((extent_[0] == 1) ? 1 : nxyz[0][i] + nGhostMinus[0] + nGhostPlus[0] + (bcx ? bcExtents[0] : 0),
                                                    (extent_[1] == 1) ? 1 : nxyz[1][j] + nGhostMinus[1] + nGhostPlus[1] + (bcy ? bcExtents[1] : 0),
                                                    (extent_[2] == 1) ? 1 : nxyz[2][k] + nGhostMinus[2] + nGhostPlus[2] + (bcz ? bcExtents[2] : 0)),
                                             bcx,
                                             bcy,
                                             bcz) );
            cumOffset[0] += nxyz[0][i];
          }
          cumOffset[1] += nxyz[1][j];
        }
        cumOffset[2] += nxyz[2][k];
      }
      return children;
    }

    ostream& operator<<(ostream& os, const MemoryWindow& w ){
      os << w.nptsGlob_ << w.offset_ << w.extent_ << w.bc_;
      return os;
    }

  } // namespace structured
} // namespace SpatialOps
