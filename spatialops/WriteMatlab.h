#ifndef SpatialOps_WriteMatlab_h
#define SpatialOps_WriteMatlab_h

#include <fstream>
#include <iomanip>

namespace SpatialOps{

  /**
   *  \file WriteMatlab
   *  \function write_matlab
   *  \brief writes a field to a matlab file
   */
  template<typename FieldT>
  void write_matlab( const FieldT& field,
                     const std::string prefix,
                     const bool includeGhost=false )
  {
    const std::string fname = "load_"+prefix+".m";
    std::ofstream fout( fname.c_str() );
    fout << "function x = load_" << prefix << "()" << std::endl;
    fout << std::scientific;
    fout.precision( 14 );
    if( includeGhost ){
      typename FieldT::const_iterator i=field.begin(), iend=field.end();
      fout << "x = [ " << *i;
      ++i;
      for( ; i!=iend; ++i )  fout << ", " << *i << std::endl;
      fout << " ];" << std::endl;
    }
    else{
      typename FieldT::const_interior_iterator i=field.interior_begin(), iend=field.interior_end();
      fout << "x = [ " << *i;
      ++i;
      for( ; i!=iend; ++i )  fout << ", " << *i << std::endl;
      fout << " ];" << std::endl;
    }
    fout.close();
  }

} // namespace SpatialOps

#endif // SpatialOps_WriteMatlab_h
