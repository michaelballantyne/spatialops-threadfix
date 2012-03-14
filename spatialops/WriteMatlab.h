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

#ifndef SpatialOps_WriteMatlab_h
#define SpatialOps_WriteMatlab_h

#include <string>
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
      typename FieldT::const_iterator
        i    = field.begin(),
        iend = field.end();
      fout << "x = [ " << *i;
      ++i;
      for( ; i!=iend; ++i )  fout << std::endl << ", " << *i;
      fout << " ];" << std::endl;
    }
    else{
      typename FieldT::const_interior_iterator
        i    = field.interior_begin(),
        iend = field.interior_end();
      fout << "x = [ " << *i;
      ++i;
      for( ; i!=iend; ++i )  fout << std::endl << ", " << *i;
      fout << " ];" << std::endl;
    }
    fout.close();
  }

} // namespace SpatialOps

#endif // SpatialOps_WriteMatlab_h
