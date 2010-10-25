#ifndef SpatialOps_FieldWriter_h
#define SpatialOps_FieldWriter_h

#include <ostream>

namespace SpatialOps{
namespace structured{

  template<typename FieldT>
  void write( std::ostream& os,
              const std::string& name,
              const FieldT& field,
              const bool includeGhost=false )
  {
    const MemoryWindow& w = includeGhost ? field.window_with_ghost() : field.window_without_ghost();
    const int n = name.size();
    os.write( reinterpret_cast<const char*>(&n), sizeof(int) );
    os.write( name.c_str(), n*sizeof(char) );
    write( os, w );
    for( size_t k=0; k<w.extent(2); ++k ){
      for( size_t j=0; j<w.extent(1); ++j ){
        // should be contiguous in x dir
        const typename FieldT::AtomicT& val = field(0,j,k);
        os.write( reinterpret_cast<const char*>(&val),
                  sizeof(typename FieldT::AtomicT) * w.extent(0) );
      }
    }
  }

  /*
    % to read in matlab:
    fid = fopen('fields.bin','r+');

    n = fread(fid,1,'int')
    nam = char(fread(fid,n,'char')')

    nglob  = fread(fid,3,'int')
    offset = fread(fid,3,'int')
    extent = fread(fid,3,'int')

    fclose(fid);
   */

//   template<typename FieldT>
//   std::ostream& operator<<( std::ostream& os, const FieldT& f )
//   {
//     const MemoryWindow& w = f.window_without_ghost();
//     os << w;
//     for( size_t k=0; k<w.extent(2); ++k ){
//       for( size_t j=0; j<w.extent(1); ++j ){
//         for( size_t i=0; i<w.extent(0); ++i ){
//           os << f(0,j,k) << ", ";
//         }
//         os << std::endl;
//       }
//       os << std::endl;
//     }
//     return os;
//   }

} // namespace structured
} // namespace SpatialOps


#endif // SpatialOps_FieldWriter_h
