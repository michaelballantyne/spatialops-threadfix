#include <spatialops/structured/stencil/CudaStencil2Bridge.h>
#include <spatialops/SpatialOpsDefs.h>

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
									 const int dEX_x,   const int dEX_y,   const int dEX_z,
									 const int dOFF_x,  const int dOFF_y,  const int dOFF_z,
									 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z,
									 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z )
		{
			const int i   = ( blockIdx.x * blockDim.x + threadIdx.x );
			const int j   = ( blockIdx.y * blockDim.y + threadIdx.y );
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
			  int s1Idx = Index3D( nx, ny, s1i, s1j, s1OFF_z );
			  int s2Idx = Index3D( nx, ny, s2i, s2j, s2OFF_z );
			  int nxny  = nx*ny; // Size of the overall base plane.

			  // For each destination valid index, ( di, dj ) in the plane, we march up
			  // the z axis and compute the resulting values per thread.
			  for( int dk = dOFF_z; dk < dEX_z; ++dk ) {
				//Destination value for a 2 point stencil
				// d[i] = low*s1[i] + high*s2[i]
				dest[ dIdx ] = ( low * src[ s1Idx ] + high * src[ s2Idx ] );

				//Bump flat index pointer by the size of a full XY plane.
				dIdx += nxny;
				s1Idx += nxny;
				s2Idx += nxny;
			  }
			}
		  }

		template< class DataType, class Dir>
		void cuda_stencil_2_apply_to_field( DataType* dest, const DataType* src,
										 DataType low,   DataType high,
										 const int nx,      const int ny,      const int nz,
										 const int dEX_x,   const int dEX_y,   const int dEX_z,
										 const int dOFF_x,  const int dOFF_y,  const int dOFF_z,
										 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z,
										 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z ){
		//Compute blocking dimensions
		dim3 dimBlock( BLOCK_DIM, BLOCK_DIM );
		dim3 dimGrid( 1, 1 );

		//Launch kernel
		__cuda_stencil_2_apply_to_field<DataType, Dir><<<dimGrid, dimBlock>>>(
									 dest, src,
									 low, high,
									 nx, ny, nz,
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
										 const int dEX_x,   const int dEX_y,   const int dEX_z, \
										 const int dOFF_x,  const int dOFF_y,  const int dOFF_z, \
										 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z, \
										 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z );

#define DECLARE_STENCIL_CUDA( TYPE, DIR) template void \
		__cuda_stencil_2_apply_to_field< TYPE , DIR >( TYPE* dest, const TYPE* src, \
										 TYPE low,   TYPE high, \
										 const int nx,      const int ny,      const int nz, \
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
