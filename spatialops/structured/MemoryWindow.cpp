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

  ostream& operator<<(ostream& os, const IntVec& v )
  {
    os << "[ " << v[0]
       << ","  << v[1]
       << ","  << v[2]
       << " ]" << endl;
    return os;
  }

} // namespace structured
} // namespace SpatialOps
