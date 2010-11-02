#include <spatialops/structured/MemoryWindow.h>

#include <ostream>
using namespace std;

namespace SpatialOps{
namespace structured{

  bool check_positive( const IntVec& v )
  {
    return (v[0]>=0) & (v[1]>=0) & (v[2]>=0);
  }
	
  MemoryWindow::MemoryWindow( const int npts[3],
                              const int offset[3],
                              const int extent[3] )
    : nptsGlob_( npts ),
      offset_( offset ),
      extent_( extent )
  {
#   ifndef NDEBUG
    assert( check_positive( nptsGlob_ ) );
    assert( check_positive( offset_   ) );
    assert( check_positive( extent_   ) );
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
    assert( check_positive( nptsGlob_ ) );
    assert( check_positive( offset_   ) );
    assert( check_positive( extent_   ) );
#   endif
  }
	
  MemoryWindow::MemoryWindow( const int npts[3] )
    : nptsGlob_( npts ), offset_(0,0,0), extent_( npts )
  {
#   ifndef NDEBUG
    assert( check_positive( nptsGlob_ ) );
    assert( check_positive( offset_   ) );
    assert( check_positive( extent_   ) );
#   endif
  }
	
  MemoryWindow::MemoryWindow( const IntVec& npts )
    : nptsGlob_( npts ), offset_(0,0,0), extent_( npts )
  {
#   ifndef NDEBUG
    assert( check_positive( nptsGlob_ ) );
    assert( check_positive( offset_   ) );
    assert( check_positive( extent_   ) );
#   endif
  }
	
  MemoryWindow::MemoryWindow( const MemoryWindow& other )
    : nptsGlob_( other.nptsGlob_ ),
      offset_  ( other.offset_   ),
      extent_  ( other.extent_   )
  {}
	
  MemoryWindow::~MemoryWindow()
  {}
	
  ostream& operator<<(ostream& os, const MemoryWindow& w )
  {
    os << w.nptsGlob_ << w.offset_ << w.extent_;
    return os;
  }

  ostream& operator<<(ostream& os, const IntVec& v )
  {
    os << "[ " << v[0]
       << ","  << v[1]
       << ","  << v[2]
       << " ]";
    return os;
  }

} // namespace structured
} // namespace SpatialOps
