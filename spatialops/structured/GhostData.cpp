/**
 *  \file   GhostDataRT.cpp
 *  \date   Jul 8, 2013
 *  \author "James C. Sutherland"
 *
 *
 * The MIT License
 *
 * Copyright (c) 2013 The University of Utah
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
 *
 */


#include "GhostData.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <sstream>

namespace SpatialOps{ namespace structured{

  inline void check_valid( const IntVec& minus, const IntVec& plus )
  {
#   ifndef NDEBUG
    for( int i=0; i<3; ++i ){
      assert( minus[i] >= 0 );
      assert(  plus[i] >= 0 );
    }
#   endif
  }

  //-----------------------------------------------------------------

  // If any value is infinite (>= GHOST_MAX), then return true
  //  Infinite ghost data is used within Nebo to account for scalars
  inline bool is_IntVec_infinite(const IntVec & values) {
    return (values[0] >= GHOST_MAX &&
            values[1] >= GHOST_MAX &&
            values[2] >= GHOST_MAX);
  }

  // If any value is infinite (>= GHOST_MAX), then return true
  //  Infinite ghost data is used within Nebo to account for scalars
  inline bool is_ghost_infinite(const IntVec & minus, const IntVec & plus) {
    return (is_IntVec_infinite(minus) &&
            is_IntVec_infinite(plus));
  }

  //-----------------------------------------------------------------

  GhostDataRT::GhostDataRT( const int nx, const int px,
                            const int ny, const int py,
                            const int nz, const int pz )
  : minus_( nx, ny, nz ),
    plus_ ( px, py, pz ),
    isInf_(is_ghost_infinite(minus_,plus_))
  {
    check_valid(minus_,plus_);
  }

  GhostDataRT::GhostDataRT( const IntVec& minus,
                            const IntVec& plus )
  : minus_( minus ),
    plus_ ( plus  ),
    isInf_(is_ghost_infinite(minus_,plus_))
  {
    check_valid(minus_,plus_);
  }

  GhostDataRT::GhostDataRT( const int n )
  : minus_( n, n, n ),
    plus_ ( n, n, n ),
    isInf_(is_ghost_infinite(minus_,plus_))
  {
    check_valid(minus_,plus_);
  }

  //-----------------------------------------------------------------

  GhostDataRT&
  GhostDataRT::operator=( const GhostDataRT& rhs )
  {
    minus_ = rhs.minus_;
    plus_  = rhs.plus_;
    check_valid(minus_,plus_);
    isInf_ = rhs.isInf_;
    return *this;
  }

  GhostDataRT::GhostDataRT( const GhostDataRT& rhs )
  {
    minus_ = rhs.minus_;
    plus_  = rhs.plus_;
    check_valid(minus_,plus_);
    isInf_ = rhs.isInf_;
  }

  //-----------------------------------------------------------------

  void
  GhostDataRT::set_minus( const IntVec& minus )
  {
    minus_ = minus;
    isInf_ = is_ghost_infinite(minus_,plus_);
  }

  void
  GhostDataRT::set_plus( const IntVec& plus )
  {
    plus_ = plus;
    isInf_ = is_ghost_infinite(minus_,plus_);
  }

  //-----------------------------------------------------------------

  GhostDataRT
  GhostDataRT::operator-( const GhostDataRT& rhs ) const
  {
    GhostDataRT g(*this);
    g -= rhs;
    return g;
  }

  GhostDataRT
  GhostDataRT::operator+( const GhostDataRT& rhs ) const
  {
    GhostDataRT g(*this);
    g += rhs;
    return g;
  }

  //-----------------------------------------------------------------

  GhostDataRT&
  GhostDataRT::operator-=( const GhostDataRT& rhs )
  {
    if(rhs.isInf_) {
      throw(std::runtime_error("Cannot use infinite ghost data on the right-hand side of subtraction."));
    } else if(!isInf_) {
      minus_ -= rhs.minus_;
      plus_  -= rhs.plus_;
      check_valid(minus_,plus_);
    }
    return *this;
  }

  GhostDataRT&
  GhostDataRT::operator+=( const GhostDataRT& rhs )
  {
    if(!isInf_) {
      if(rhs.isInf_) {
        *this = rhs;
      } else {
        minus_ += rhs.minus_;
        plus_  += rhs.plus_;
        check_valid(minus_,plus_);
      }
    }
    return *this;
  }

  //-----------------------------------------------------------------

  bool
  GhostDataRT::operator==( const GhostDataRT& rhs ) const
  {
    return (minus_ == rhs.minus_) && (plus_ == rhs.plus_);
  }

  //-----------------------------------------------------------------

  GhostDataRT
  min( const GhostDataRT& first, const GhostDataRT& second )
  {
      return GhostDataRT(min(first.get_minus(), second.get_minus()),
                         min(first.get_plus(), second.get_plus()));
  }

  //-----------------------------------------------------------------

  GhostDataRT
  point_to_ghost( const IntVec& given )
  {
      return GhostDataRT((given[0] < 0 ? - given[0] : 0),
                         (given[0] > 0 ?   given[0] : 0),
                         (given[1] < 0 ? - given[1] : 0),
                         (given[1] > 0 ?   given[1] : 0),
                         (given[2] < 0 ? - given[2] : 0),
                         (given[2] > 0 ?   given[2] : 0));
  }

  //-----------------------------------------------------------------

  std::ostream& operator<<( std::ostream& out, const GhostDataRT& gd )
  {
    out << "{ " << gd.get_minus() << " " << gd.get_plus() << " }";
    return out;
  }

  //-----------------------------------------------------------------
}
}
