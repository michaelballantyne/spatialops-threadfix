/*
 * Stencil2_def.h
 * Defines the implementation of stencil2.  This allows for further
 * instantiation by applications that use this library.  Basic Stencil2
 * instantiations are contained in Stencil2.cpp
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
 *
 */

#ifndef STENCIL2_DEF_H_
#define STENCIL2_DEF_H_


#include <spatialops/SpatialOpsConfigure.h>
#ifdef ENABLE_CUDA
#include "CudaStencil2.h"
#else
#include "Stencil2.h"
#endif
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/FieldExpressions.h>

namespace SpatialOps{ namespace structured{


  template< typename OperatorT, typename SrcT, typename DestT >
  Stencil2<OperatorT,SrcT,DestT>::
  Stencil2( const double coefLo, const double coefHi )
    : coefLo_( coefLo ),
      coefHi_( coefHi ),
      coefList_( build_two_point_coef_list(coefLo, coefHi) )
  {}

  //------------------------------------------------------------------

  template< typename OperatorT, typename SrcT, typename DestT >
  Stencil2<OperatorT,SrcT,DestT>::
  ~Stencil2()
  {}

  //------------------------------------------------------------------

  template< typename OperatorT, typename SrcT, typename DestT >
  void
  Stencil2<OperatorT,SrcT,DestT>::
  apply_to_field( const SrcT& src, DestT& dest ) const
  {
    switch( dest.memory_device_type() ){
      case LOCAL_RAM:
          {
              dest <<= operator()(src);
          }
        break;

#ifdef ENABLE_CUDA
      case EXTERNAL_CUDA_GPU:
        cuda_stencil_2_apply_to_field_helper<OperatorT, SrcT, DestT>( src, dest, coefLo_, coefHi_ );
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

  //==================================================================
  // Explicit template instantiation
  //
#define DECLARE_STENCIL( OP, SRC, DEST )        \
  template class Stencil2< OP, SRC, DEST >;

#define DECLARE_BASIC_VARIANTS( VOL )                          \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::XFace )   \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::YFace )   \
  DECLARE_STENCIL( Interpolant, VOL, FaceTypes<VOL>::ZFace )   \
  DECLARE_STENCIL( Interpolant, FaceTypes<VOL>::XFace, VOL )   \
  DECLARE_STENCIL( Interpolant, FaceTypes<VOL>::YFace, VOL )   \
  DECLARE_STENCIL( Interpolant, FaceTypes<VOL>::ZFace, VOL )   \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::XFace )   \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::YFace )   \
  DECLARE_STENCIL( Gradient,    VOL, FaceTypes<VOL>::ZFace )   \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::XFace, VOL )   \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::YFace, VOL )   \
  DECLARE_STENCIL( Divergence,  FaceTypes<VOL>::ZFace, VOL )

} // namespace structured
} // namespace SpatialOps

#endif /* STENCIL2_DEF_H_ */
