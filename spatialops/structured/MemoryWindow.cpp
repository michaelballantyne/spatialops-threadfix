#include <spatialops/structured/MemoryWindow.h>

#include <ostream>
using namespace std;

namespace SpatialOps{
namespace structured{

  ostream& operator<<(ostream& os, const MemoryWindow& w )
  {
    os << w.nptsGlob_ << w.offset_ << w.extent_;
    return os;
  }

  void write( ostream& os, const MemoryWindow& w )
  {
    write( os, w.nptsGlob_ );
    write( os, w.offset_   );
    write( os, w.extent_   );
  }

  ostream& operator<<(ostream& os, const IntVec& v )
  {
    os << "[ " << v[0]
       << ","  << v[1]
       << ","  << v[2]
       << " ]" << endl;
    return os;
  }

  void write(ostream& os, const IntVec& v )
  {
    os.write( reinterpret_cast<const char*>(&(v.ijk[0])), 3*sizeof(int) );
  }

} // namespace structured
} // namespace SpatialOps
