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

#ifndef CUDASTENCIL4_CUH_
#define CUDASTENCIL4_CUH_

#include <spatialops/structured/stencil/CudaStencil4Bridge.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/stencil/Stencil4.h>


namespace SpatialOps {
  namespace structured {
    template <typename OperatorType, typename SrcType, typename DestType>
      inline void cuda_stencil_4_apply_to_field_helper( const SrcType& src,
 					                DestType& dest,
						        double const Coef1,
						        double const Coef2,
                                                        double const Coef3,
    							double const Coef4 )

{
 //Shortcut !
s4detail::ExtentsAndOffsets<SrcType, DestType> typedef Extents;

// Gather destination window information to generate blocking info
const MemoryWindow& wdest = dest.window_with_ghost();
const MemoryWindow& ws    =  src.window_with_ghost();

IntVec dOFF  = Extents::DestOffset::int_vec();
IntVec s1OFF = Extents::Src1Offset::int_vec();
IntVec s2OFF = Extents::Src2Offset::int_vec();
IntVec s3OFF = Extents::Src3Offset::int_vec();
IntVec s4OFF = Extents::Src4Offset::int_vec();

IntVec wEX = wdest.extent();
IntVec dEX = wEX + Extents::DestOffset::int_vec() + dest.get_ghost_data().has_bc() * Extents::DestOffset::int_vec();
IntVec sEX = ws.glob_dim();

//Call interface function -- hack to avoid nvcc meta-template failures
cuda_stencil_4_apply_to_field< typename DestType::AtomicT, typename Extents::Dir1 >(
    dest.field_values( EXTERNAL_CUDA_GPU, dest.device_index() ),
    src.field_values( EXTERNAL_CUDA_GPU, dest.device_index() ),
    Coef1, Coef2, Coef3, Coef4,              // Stencil 4 Coefficients
    wEX[0],   wEX[1],   wEX[2],              // Global field dimensions
    sEX[0],   sEX[1],   sEX[2],              // Global Source Extents
    dEX[0],   dEX[1],   dEX[2],              // Destination extents dEX <= wdest.extent
    dOFF[0],  dOFF[1],  dOFF[2],             // Destination point offsets
    s1OFF[0], s1OFF[1], s1OFF[2],            // Source 1 point offsets
    s2OFF[0], s2OFF[1], s2OFF[2],            // Source 2 point offsets
    s3OFF[0], s3OFF[1], s3OFF[2],            // Source 3 point offsets
    s4OFF[0], s4OFF[1], s4OFF[2]             // Source 4 point offsets
    );
  }
 } // structured
} // SpatialOps

#undef Index3D

#endif /* CUDASTENCIL4_CUH */
