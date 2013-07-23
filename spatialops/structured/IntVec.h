/**
 *  \file   IntVec.h
 *
 *  \date   Sep 28, 2011
 *  \author James C. Sutherland
 *
 *
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

#ifndef SpatialOps_IntVec_h
#define SpatialOps_IntVec_h

#include <ostream>
#include <vector>
#include <cassert>


namespace SpatialOps{
  namespace structured{

    /**
     *  \class IntVec
     *  \author James C. Sutherland
     *  \ingroup structured
     *  \brief provides a lightweight class to deal with a 3D vector of integers.
     */
    class IntVec
    {
      friend std::ostream& operator<<(std::ostream&, const IntVec&);

      int ijk[3];

    public:

      IntVec(){ ijk[0]=0; ijk[1]=0; ijk[2]=0; }

      inline IntVec( const int i, const int j, const int k ){
        ijk[0]=i; ijk[1]=j; ijk[2]=k;
      }

      inline IntVec( const int vec[3] ){
        ijk[0]=vec[0]; ijk[1]=vec[1]; ijk[2]=vec[2];
      }

      IntVec( const std::vector<int>& vec ){
        ijk[0]=vec[0]; ijk[1]=vec[1]; ijk[2]=vec[2];
      }

      inline IntVec( const IntVec& x ){
        ijk[0]=x.ijk[0];  ijk[1]=x.ijk[1];  ijk[2]=x.ijk[2];
      }

      inline int  operator[](const size_t i) const{ assert(i<3); return ijk[i]; }
      inline int& operator[](const size_t i) { assert(i<3); return ijk[i]; }

      IntVec& operator=(const IntVec& x){
        for( size_t i=0; i<3; ++i ) ijk[i] = x.ijk[i];
        return *this;
      }

      inline bool operator==(const IntVec& v) const{
        return (ijk[0]==v.ijk[0]) & (ijk[1]==v.ijk[1]) & (ijk[2]==v.ijk[2]);
      }
      inline bool operator!=(const IntVec& v) const{
        return (ijk[0]!=v.ijk[0]) | (ijk[1]!=v.ijk[1]) | (ijk[2]!=v.ijk[2]);
      }

      inline IntVec operator+( const IntVec& v ) const{
        return IntVec( ijk[0] + v.ijk[0],
                       ijk[1] + v.ijk[1],
                       ijk[2] + v.ijk[2] );
      }
      inline IntVec operator-( const IntVec& v ) const{
        return IntVec( ijk[0] - v.ijk[0],
                       ijk[1] - v.ijk[1],
                       ijk[2] - v.ijk[2] );
      }
      inline IntVec operator*( const IntVec& v ) const{
        return IntVec( ijk[0] * v.ijk[0],
                       ijk[1] * v.ijk[1],
                       ijk[2] * v.ijk[2] );
      }
      inline IntVec operator/( const IntVec& v ) const{
        return IntVec( ijk[0] / v.ijk[0],
                       ijk[1] / v.ijk[1],
                       ijk[2] / v.ijk[2] );
      }

      inline IntVec operator*( const int v ) const{
        return IntVec(ijk[0] * v,
                      ijk[1] * v,
                      ijk[2] * v);
      }

      inline IntVec& operator+=( const IntVec& v ){
        ijk[0] += v.ijk[0];
        ijk[1] += v.ijk[1];
        ijk[2] += v.ijk[2];
        return *this;
      }
      inline IntVec& operator-=( const IntVec& v ){
        ijk[0] -= v.ijk[0];
        ijk[1] -= v.ijk[1];
        ijk[2] -= v.ijk[2];
        return *this;
      }
    };

    inline std::ostream& operator<<( std::ostream& os, const IntVec& v ){
      os << "[ " << v[0] << ","  << v[1] << ","  << v[2] << " ]";
      return os;
    }

    /**
     * @param dim the 3D layout of points
     * @param loc the 3D index
     * @return the 1D (flat) index
     * Assumes a fortran-style layout (e.g. first index varies fastest)
     */
    inline int ijk_to_flat( const IntVec& dim, const IntVec& loc ){
      assert( loc[0] < dim[0] && loc[1] < dim[1] && loc[2] < dim[2] );
      assert( loc[0] >= 0 && loc[1] >= 0 && loc[2] >= 0 );
      return loc[0] + loc[1]*dim[0] + loc[2]*dim[0]*dim[1];
    }

    /**
     * @param dim the 3D layout of points
     * @param pt the 1d (flat) index
     * @return the 3D index
     * Assumes a fortran-style layout (e.g. first index varies fastest)
     */
    inline IntVec flat_to_ijk( const IntVec& dim, const int pt ){
      assert( pt >= 0 );
      assert( pt < dim[0]*dim[1]*dim[2] );
      return IntVec( pt % dim[0], (pt / dim[0]) % dim[1], pt / (dim[0]*dim[1]) );
    }

  } // namespace structured
} // namespace SpatialOps

#endif /* SpatialOps_IntVec_h */
