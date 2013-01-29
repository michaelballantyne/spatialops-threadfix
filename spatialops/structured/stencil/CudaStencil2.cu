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

#include <spatialops/structured/stencil/CudaStencil2Bridge.h>
#include <spatialops/SpatialOpsDefs.h>
//#include <iostream>
//using std::cout;

#define BLOCK_DIM 16
#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))

namespace SpatialOps {
	namespace structured {

		using namespace SpatialOps;
		using namespace SpatialOps::structured;

		template< class DataType, class Dir>
		__global__ void __cuda_stencil_2_apply_to_field( DataType* dest, const DataType* src,
									 DataType low,   DataType high,
									 const int nx,      const int ny,      const int nz,
                                                                         const int sEX_x,   const int sEX_y,   const int sEX_z, 
									 const int dEX_x,   const int dEX_y,   const int dEX_z,
									 const int dOFF_x,  const int dOFF_y,  const int dOFF_z,
									 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z,
									 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z )
		{
			const uint i  = ( blockIdx.x * blockDim.x + threadIdx.x );
			const uint j  = ( blockIdx.y * blockDim.y + threadIdx.y );
 
			const int di  = i + dOFF_x;
			const int dj  = j + dOFF_y;
  
			const int s1i = i + s1OFF_x;
			const int s1j = j + s1OFF_y;
			const int s2i = i + s2OFF_x;
			const int s2j = j + s2OFF_y;

			//////////////////////////////////////////////////////////////////////////////////////
			//          _________________________                    _________________________
			//         / _____________________  /|                  / _____________________  /|
			//        / / ___________________/ / |                 / / ___________________/ / |
			//       / / /| |               / /  |                / / /| |               / /  |
			//      / / / | |    *(s2i,s2j,s2k). |               / / / | |              / / | |
			//     / / /| | |             / / /| |              / / /| | | *(di,dj,dk) / / /| |
			//    / / / | |*(s1i,s1j,s1k)/ / / | |             / / / | | |            / / / | |
			//   / / /  | | |           / / /| | |       \    / / /  | | |           / / /| | |
			//  / /_/__________________/ / / | | | ====== \  / /_/__________________/ / / | | |
			// /________________________/ /  | | | ====== / /________________________/ /  | | |
			// | ______________________ | |  | | |       /  | ______________________ | |  | | |
			// | | |    | | |_________| | |__| | |          | | |    | | |_________| | |__| | |
			// | | |    | |___________| | |____| |          | | |    | |___________| | |____| |
			// | | |   / / ___________| | |_  / /           | | |   / / ___________| | |_  / /
			// | | |  / / /           | Z |/ / /            | | |  / / /           | Z |/ / /
			// | | | / / /   COMPUTE  | | | / /             | | | / / /   COMPUTE  | | | / /
			// | | |/ / /     PLANE   | | |/ Y              | | |/ / /     PLANE   | | |/ Y
			// | | | / /              | | ' /               | | | / /              | | ' /
			// | | |/_/_______________| |  /                | | |/_/_______________| |  /
			// | |____________________| | /                 | |____________________| | /
			// |____________X___________|/                  |____________X___________|/
			///////////////////////////////////////////////////////////////////////////////////////

			if( di < dEX_x && dj < dEX_y ) {
			  int dIdx  = Index3D( nx, ny, di, dj, dOFF_z );
			  int s1Idx = Index3D( sEX_x, sEX_y, s1i, s1j, s1OFF_z );
			  int s2Idx = Index3D( sEX_x, sEX_y, s2i, s2j, s2OFF_z );
			  int nxny  = nx*ny;                               // Size of the overall base plane.
                          int sxsy = sEX_x*sEX_y; 

			  // For each destination valid index, ( di, dj ) in the plane, we march up
			  // the z axis and compute the resulting values per thread.
			  for( int dk = dOFF_z; dk < dEX_z; ++dk ) {
				//Destination value for a 2 point stencil
				// d[i] = low*s1[i] + high*s2[i]
				dest[ dIdx ] = ( low * src[ s1Idx ] + high * src[ s2Idx ] );

				//Bump flat index pointer by the size of a full XY plane.
				dIdx  += nxny;
				s1Idx += sxsy;
				s2Idx += sxsy;
			  }
			}
		  }

		template< class DataType, class Dir>
		void cuda_stencil_2_apply_to_field( DataType* dest, const DataType* src,
										 DataType low,   DataType high,
										 const int nx,      const int ny,      const int nz,
                                                                                 const int sEX_x,   const int sEX_y,   const int sEX_z,
										 const int dEX_x,   const int dEX_y,   const int dEX_z,
										 const int dOFF_x,  const int dOFF_y,  const int dOFF_z,
										 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z,
										 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z ){
			//Compute Grid and Block dimensions
			unsigned int gDimx = (dEX_x / BLOCK_DIM ) + ( dEX_x % BLOCK_DIM > 0 ? 1 : 0 );
			unsigned int gDimy = (dEX_y / BLOCK_DIM ) + ( dEX_y % BLOCK_DIM > 0 ? 1 : 0 );
			dim3 dimBlock( BLOCK_DIM, BLOCK_DIM );
			dim3 dimGrid( gDimx, gDimy );

			//Launch kernel
              		__cuda_stencil_2_apply_to_field<DataType, Dir><<<dimGrid, dimBlock>>>(
										 dest, src,
										 low, high,
										 nx, ny, nz,
                                                                                 sEX_x, sEX_y, sEX_z,
										 dEX_x, dEX_y, dEX_z,
										 dOFF_x, dOFF_y, dOFF_z,
										 s1OFF_x, s1OFF_y, s1OFF_z,
										 s2OFF_x, s2OFF_y, s2OFF_z );

			cudaDeviceSynchronize();
		}

#define DECLARE_STENCIL( TYPE, DIR) template void \
		cuda_stencil_2_apply_to_field< TYPE , DIR >( TYPE* dest, const TYPE* src, \
										 TYPE low,   TYPE high, \
										 const int nx,      const int ny,      const int nz, \
                                                                                 const int sEX_x, const int sEX_y, const int sEX_z, \
                                                                                 const int dEX_x,   const int dEX_y,   const int dEX_z, \
										 const int dOFF_x,  const int dOFF_y,  const int dOFF_z, \
										 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z, \
										 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z );

#define DECLARE_STENCIL_CUDA( TYPE, DIR) template void \
		__global__ __cuda_stencil_2_apply_to_field< TYPE , DIR >( TYPE* dest, const TYPE* src, \
										 TYPE low,   TYPE high, \
										 const int nx,      const int ny,      const int nz, \
                                                                                 const int sEX_x,   const int sEX_y,   const int sEX_z,\
                                                                                 const int dEX_x,   const int dEX_y,   const int dEX_z, \
										 const int dOFF_x,  const int dOFF_y,  const int dOFF_z, \
										 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z, \
										 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z );

		DECLARE_STENCIL(float, NODIR);
		DECLARE_STENCIL(float, XDIR);
		DECLARE_STENCIL(float, YDIR);
		DECLARE_STENCIL(float, ZDIR);

		DECLARE_STENCIL(double, NODIR);
		DECLARE_STENCIL(double, XDIR);
		DECLARE_STENCIL(double, YDIR);
		DECLARE_STENCIL(double, ZDIR);

		DECLARE_STENCIL_CUDA(float, NODIR);
		DECLARE_STENCIL_CUDA(float, XDIR);
		DECLARE_STENCIL_CUDA(float, YDIR);
		DECLARE_STENCIL_CUDA(float, ZDIR);

		DECLARE_STENCIL_CUDA(double, NODIR);
		DECLARE_STENCIL_CUDA(double, XDIR);
		DECLARE_STENCIL_CUDA(double, YDIR);
		DECLARE_STENCIL_CUDA(double, ZDIR);
	}
}
