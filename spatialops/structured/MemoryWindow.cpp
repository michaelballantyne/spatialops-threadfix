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

#include <spatialops/structured/MemoryWindow.h>

#include <ostream>
using namespace std;

namespace SpatialOps{
  namespace structured{

    inline bool check_positive( const IntVec& v ){ return (v[0]> 0) & (v[1]> 0) & (v[2]> 0); }
    inline bool check_ge_zero ( const IntVec& v ){ return (v[0]>=0) & (v[1]>=0) & (v[2]>=0); }

    //---------------------------------------------------------------

    MemoryWindow::MemoryWindow( const int npts[3],
                                const int offset[3],
                                const int extent[3] )
    : nptsGlob_( npts ),
      offset_( offset ),
      extent_( extent )
    {
#   ifndef NDEBUG
      assert( sanity_check() );
#   endif
    }

    MemoryWindow::MemoryWindow( const IntVec& npts,
                                const IntVec& offset,
                                const IntVec& extent )
    : nptsGlob_( npts ),
      offset_( offset ),
      extent_( extent )
    {
#   ifndef NDEBUG
      assert( sanity_check() );
#   endif
   }

    MemoryWindow::MemoryWindow( const int npts[3] )
    : nptsGlob_( npts ), offset_(0,0,0), extent_( npts )
    {
#   ifndef NDEBUG
      assert( sanity_check() );
#   endif
    }

    MemoryWindow::MemoryWindow( const IntVec& npts )
    : nptsGlob_( npts ), offset_(0,0,0), extent_( npts )
    {
#   ifndef NDEBUG
      assert( sanity_check() );
#   endif
    }

    //---------------------------------------------------------------

    MemoryWindow&
    MemoryWindow::operator=( const MemoryWindow& other )
    {
      nptsGlob_ = other.nptsGlob_;
      offset_   = other.offset_;
      extent_   = other.extent_;
      return *this;
    }

    //---------------------------------------------------------------

    MemoryWindow::~MemoryWindow()
    {}

    //---------------------------------------------------------------

    MemoryWindow
    MemoryWindow::refine( const IntVec& splitPattern,
                          const IntVec& location ) const
    {
#     ifndef NDEBUG
      for( size_t i=0; i<3; ++i ){
        assert( extent_[i] >= splitPattern[i] );
        assert( splitPattern[i] > 0 );
        assert( extent_[i] + offset_[i] <= nptsGlob_[i] );
        assert( location[i] < splitPattern[i] );
        assert( location[i] >= 0 );
      }
#     endif
      // try to make windows close to same size
      const IntVec nextra = IntVec( extent_[0] % splitPattern[0],
                                    extent_[1] % splitPattern[1],
                                    extent_[2] % splitPattern[2] );
      const IntVec stdExtent = IntVec( extent_[0] / splitPattern[0],
                                       extent_[1] / splitPattern[1],
                                       extent_[2] / splitPattern[2] );
      const IntVec baseOffset = IntVec( offset_[0] + stdExtent[0] * location[0] + (location[0] < nextra[0] ? location[0] : nextra[0]),
                                        offset_[1] + stdExtent[1] * location[1] + (location[1] < nextra[1] ? location[1] : nextra[1]),
                                        offset_[2] + stdExtent[2] * location[2] + (location[2] < nextra[2] ? location[2] : nextra[2]) );
      const IntVec baseExtent = IntVec( stdExtent[0] + (location[0] < nextra[0] ? 1 : 0),
                                        stdExtent[1] + (location[1] < nextra[1] ? 1 : 0),
                                        stdExtent[2] + (location[2] < nextra[2] ? 1 : 0) );

      return MemoryWindow( nptsGlob_, baseOffset, baseExtent );
    };

    //---------------------------------------------------------------

    ostream& operator<<(ostream& os, const MemoryWindow& w ){
      os << w.nptsGlob_ << w.offset_ << w.extent_;
      return os;
    }

    bool
    MemoryWindow::sanity_check() const
    {
      return check_positive( nptsGlob_ ) &&
             check_ge_zero ( offset_   ) &&
             check_ge_zero ( extent_   ) &&
             check_ge_zero ( nptsGlob_ - extent_ );
    }

  } // namespace structured
} // namespace SpatialOps
