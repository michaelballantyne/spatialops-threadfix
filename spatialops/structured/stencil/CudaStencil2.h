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

#ifndef CUDASTENCIL2_CUH_
#define CUDASTENCIL2_CUH_

#include <spatialops/structured/stencil/CudaStencil2Bridge.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/stencil/Stencil2.h>

// Used to avoid obnoxious error parsing in eclipse. __CDT_PARSER__ is an
// eclipse defined variable used by the syntax parser.

namespace SpatialOps {
  namespace structured {
    template< class OperatorType, class SrcType, class DestType >
      inline void
      cuda_stencil_2_apply_to_field_helper( const SrcType& src,
                                            DestType& dest,
                                            double const low,
                                            double const high )
      {
        //Shortcut!
        s2detail::ExtentsAndOffsets<SrcType, DestType> typedef Extents;

        // Gather destination window information to generate blocking info
        const MemoryWindow& wdest = dest.window_with_ghost();
        IntVec dOFF  = Extents::DestOffset::int_vec();
        IntVec s1OFF = Extents::Src1Offset::int_vec();
        IntVec s2OFF = Extents::Src2Offset::int_vec();

        IntVec wEX = wdest.extent();
        IntVec dEX = wEX + Extents::DestExtent::int_vec()
                   + wdest.has_bc() * Extents::DestExtentBC::int_vec();

        //Call interface function -- hack to avoid nvcc meta-template failures
        cuda_stencil_2_apply_to_field< typename DestType::AtomicT, typename Extents::Dir >(
            dest.ext_field_values(),
            src.field_values_consumer(EXTERNAL_CUDA_GPU, dest.device_index()),
            low, high,                  	//Stencil Coeffcients
            wEX[0], wEX[1], wEX[2],			// Global field dimensions
            dEX[0], dEX[1], dEX[3],         //Destination extents dEX <= wdest.extent
            dOFF[0], dOFF[1], dOFF[2],      //Destination point offsets
            s1OFF[0], s1OFF[1], s1OFF[2],   //Source 1 point offsets
            s2OFF[0], s2OFF[1], s2OFF[2]    //Source 2 point offsets
        );
      }
  }  // structured
}  // SpatialOps

#undef Index3D

#endif /* CUDASTENCIL2_H_ */
