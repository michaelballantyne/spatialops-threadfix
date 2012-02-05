/*
 * CudaStencil2Bridge.h
 *
 *  Created on: Jan 29, 2012
 *      Author: Devin Robison
 */

#ifndef CUDASTENCIL2BRIDGE_H_
#define CUDASTENCIL2BRIDGE_H_

namespace SpatialOps {
	namespace structured {
		template< class DataType, class Dir>
		void cuda_stencil_2_apply_to_field( DataType* dest, const DataType* src,
									 DataType low,   DataType high,
									 const int nx,      const int ny,      const int nz,
									 const int dEX_x,   const int dEX_y,   const int dEX_z,
									 const int dOFF_x,  const int dOFF_y,  const int dOFF_z,
									 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z,
									 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z );

	}
}
#endif /* CUDASTENCIL2BRIDGE_H_ */
