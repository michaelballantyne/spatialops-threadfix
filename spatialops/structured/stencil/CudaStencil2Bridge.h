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

#ifndef CUDASTENCIL2BRIDGE_H_
#define CUDASTENCIL2BRIDGE_H_

namespace SpatialOps {
	namespace structured {
		template< class DataType, class Dir>
		void cuda_stencil_2_apply_to_field( DataType* dest, const DataType* src,
									 DataType low,   DataType high,
									 const int nx,      const int ny,      const int nz,
                                                                         const int sEX_x,   const int sEX_y,   const int sEX_z, 
									 const int dEX_x,   const int dEX_y,   const int dEX_z,
									 const int dOFF_x,  const int dOFF_y,  const int dOFF_z,
									 const int s1OFF_x, const int s1OFF_y, const int s1OFF_z,
									 const int s2OFF_x, const int s2OFF_y, const int s2OFF_z );

	}
}
#endif /* CUDASTENCIL2BRIDGE_H_ */
