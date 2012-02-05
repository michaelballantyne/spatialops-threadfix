/*
 * CudaStencil2.h
 *
 *  Created on: Jan 24, 2012
 *      Author: Devin Robison
 */

#ifndef CUDASTENCIL2_CUH_
#define CUDASTENCIL2_CUH_

#include <spatialops/structured/stencil/CudaStencil2Bridge.h>
#include <spatialops/structured/MemoryTypes.h>
#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/stencil/stencil2.h>

// Used to avoid obnoxious error parsing in eclipse. __CDT_PARSER__ is an
// eclipse defined variable used by the syntax parser.
#ifdef __CDT_PARSER__
#define __host__
#define __global__
#define __device__
#define __shared__
#endif

namespace SpatialOps {
  namespace structured {
    template< class OperatorType, class SrcType, class DestType >
      inline void
      cuda_stencil_2_apply_to_field_helper( SrcType src,
                                            DestType dest,
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
            src.ext_field_values_consumer(EXTERNAL_CUDA_GPU, dest.device_index()),
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
