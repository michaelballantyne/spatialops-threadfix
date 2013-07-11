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

  //-----------------------------------------------------------------

  GhostDataRT::GhostDataRT( const int nx, const int px,
                            const int ny, const int py,
                            const int nz, const int pz )
  : minus_( nx, ny, nz ),
    plus_ ( px, py, pz )
  {}

  GhostDataRT::GhostDataRT( const IntVec& minus,
                            const IntVec& plus )
  : minus_( minus ),
    plus_ ( plus  )
  {}

  GhostDataRT::GhostDataRT( const int n )
  : minus_( n, n, n ),
    plus_ ( n, n, n )
  {}

  //-----------------------------------------------------------------

  GhostDataRT&
  GhostDataRT::operator=( const GhostDataRT& rhs )
  {
    minus_ = rhs.minus_;
    plus_  = rhs.plus_;
    return *this;
  }

  GhostDataRT::GhostDataRT( const GhostDataRT& rhs )
  {
    minus_ = rhs.minus_;
    plus_  = rhs.plus_;
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
    minus_ -= rhs.minus_;
    plus_  -= rhs.plus_;
    return *this;
  }

  GhostDataRT&
  GhostDataRT::operator+=( const GhostDataRT& rhs )
  {
    minus_ += rhs.minus_;
    plus_  += rhs.plus_;
    return *this;
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
