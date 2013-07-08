/**
 *  \file   Stencil4_def.h
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

#include <spatialops/SpatialOpsConfigure.h>

#ifdef ENABLE_CUDA
#include "CudaStencil4.h"
#else
#include "Stencil4.h"
#endif

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/Nebo.h>


#ifndef STENCIL4_DEF_H_
#define STENCIL4_DEF_H_

namespace SpatialOps{
  namespace structured{

    //---------------------------------------------------------------

    template< typename OpT, typename SrcT, typename DestT >
    Stencil4<OpT,SrcT,DestT>::
    Stencil4( const double coef1,
        const double coef2,
        const double coef3,
        const double coef4 )
        : coef1_( coef1 ),
          coef2_( coef2 ),
          coef3_( coef3 ),
          coef4_( coef4 ),
          coefCollection_( build_four_point_coef_collection(coef1, coef2, coef3, coef4) )
          {}

    //---------------------------------------------------------------

    template< typename OpT, typename SrcT, typename DestT >
    void Stencil4<OpT,SrcT,DestT>::apply_to_field( const SrcT& src, DestT& dest ) const
    {   switch( dest.memory_device_type() ){
    case LOCAL_RAM:
      dest <<= operator()(src);
      break;
#ifdef ENABLE_CUDA
    case EXTERNAL_CUDA_GPU:
      cuda_stencil_4_apply_to_field_helper<OpT,SrcT,DestT>( src, dest, coef1_, coef2_, coef3_, coef4_ );
      break;
#endif
    default:{
      std::ostringstream msg;
      msg << "Destination field has unsupported device type ( "
          << DeviceTypeTools::get_memory_type_description( dest.memory_device_type() )
      << " )\n";
      msg << "\t - " << __FILE__ << " : " << __LINE__;
      throw(std::runtime_error(msg.str()));
    }
    }
    }

  //---------------------------------------------------------------

  }
}

#endif /* STENCIL4_DEF_H_ */
