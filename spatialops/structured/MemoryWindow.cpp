#include <spatialops/structured/MemoryWindow.h>

#include <iostream>
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
                         const IntVec npad,
                         const IntVec bcExtents ) const
    {
      const IntVec extent = extent_ - npad * 2 - bcExtents;
      const IntVec offset = offset_ + npad;

#     ifndef NDEBUG
      for( size_t i=0; i<3; ++i ){
    	std::cout << "(" << i << ") Extent: " << extent[i] << " splitPattern: " << splitPattern[i] << std::endl;
        assert( extent[i] >= splitPattern[i] );
        assert( extent[i] + offset[i] <= nptsGlob_[i] );
        assert( splitPattern[i] > 0 );
        assert( npad[i] >= 0);
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
                                             cumOffset - npad,
                                             IntVec(nxyz[0][i] + npad[0] * 2 + (bcx ? bcExtents[0] : 0),
                                                    nxyz[1][j] + npad[1] * 2 + (bcy ? bcExtents[1] : 0),
                                                    nxyz[2][k] + npad[2] * 2 + (bcz ? bcExtents[2] : 0)),
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
