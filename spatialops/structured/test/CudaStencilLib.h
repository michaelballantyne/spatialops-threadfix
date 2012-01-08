/*
 * CudaStencilLib.h
 *
 *  Created on: Jan 5, 2012
 *      Author: devin
 */

//TODO
// dx,dy,dz specified because I'm not sure where else to get them right now.
#ifndef CUDASTENCILLIB_H_
#define CUDASTENCILLIB_H_

#include <cuda.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <spatialops/structured/SpatialField.h>

// Used to avoid obnoxious error parsing in eclipse. __CDT_PARSER__ is an
// eclipse defined variable used by the syntax parser.
#ifdef __CDT_PARSER__
#define __host__
#define __global__
#define __device__
#define __shared__
#endif




// ( 256 threads per block )
#define BLOCK_DIM 16
#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))
#define Index3DG(_nx, _ny, _g, _i, _j, _k) (((_i)+(_g)) + _nx*((_j)+(_g) + _ny*((_k)+(_g)) ))

void __global__ _div_float_slow( float* global_data_in,
                                 float* global_data_out,
                                 float dx, float dy, float dz,
                                 int nx, int ny, int nz){
  int i = blockIdx.x*blockDim.x+threadIdx.x;
  int j = blockIdx.y*blockDim.y+threadIdx.y;

  //Temporarily -> don't deal with computing boundary points
  if( i>0 && j>0 && (i<nx-1) && (j<ny-1) )
  {
          for(int k=1;k<nz-1;k++)/* For each frame along the z-axis */
          {
                  float tijk = 2*global_data_in[ Index3D( nx, ny, i,j,k) ];
                  global_data_out[Index3D (nx, ny, i, j, k)] =
                  // X direction
                  (
                      global_data_in[Index3D (nx, ny, i - 1, j, k)] +
                      global_data_in[Index3D (nx, ny, i + 1, j, k)] -
                      tijk
                  ) / ( dx * dx )
                  +
                  // Y direction
                  (
                      global_data_in[Index3D (nx, ny, i, j - 1, k)] +
                      global_data_in[Index3D (nx, ny, i, j + 1, k)] -
                      tijk
                  ) / ( dy * dy )
                  +
                  // Z direction
                  (
                      global_data_in[Index3D (nx, ny, i, j, k - 1)] +
                      global_data_in[Index3D (nx, ny, i, j, k + 1)] -
                      tijk
                  ) / ( dz * dz );
          }
  }
}

void __global__ _div_float_opt1( float* global_data_in,
                                 float* global_data_out,
                                 float dx, float dy, float dz,
                                 int nx, int ny, int nz){
  int i = blockIdx.x*blockDim.x+threadIdx.x;
  int j = blockIdx.y*blockDim.y+threadIdx.y;

  //Temporarily -> don't deal with computing boundary points
  float zm1, zc, zp1;
  zc = global_data_in[Index3D( nx, ny, i, j, 0)];
  zp1 = global_data_in[Index3D( nx, ny, i, j, 1 )];

  if( i>0 && j>0 && (i<nx-1) && (j<ny-1) )
  {
          for(int k=1;k<nz-1;k++)/* For each frame along the z-axis */
          {
                  zm1 = zc;
                  zc = zp1;
                  zp1 = global_data_in[Index3D( nx, ny, i, j, k + 1 )];
                  float tijk = 2*zc;

                  global_data_out[Index3D (nx, ny, i, j, k)] =
                  // X direction
                  (
                      global_data_in[Index3D (nx, ny, i - 1, j, k)] +
                      global_data_in[Index3D (nx, ny, i + 1, j, k)] -
                      tijk
                  ) / ( dx * dx )
                  +
                  // Y direction
                  (
                      global_data_in[Index3D (nx, ny, i, j - 1, k)] +
                      global_data_in[Index3D (nx, ny, i, j + 1, k)] -
                      tijk
                  ) / ( dy * dy )
                  +
                  // Z direction
                  (
                      zm1 +
                      zp1 -
                      tijk
                  ) / ( dz * dz );
          }
  }
}

void __global__ _div_float_opt2( float* global_data_in,
                                 float* global_data_out,
                                 float dx, float dy, float dz,
                                 int nx, int ny, int nz){

  int ghost = 1; // static, because shared memory allocation has to be fixed.

  /** Compute inner tile variables **/
  int inner_dim    = BLOCK_DIM;                         // Dimensions of the inner tile
  int inner_x      = ( BlockIdx.x * inner_dim ) + ghost;// Inner tile absolute 'x' val
  int inner_y      = ( BlockIdx.y * inner_dim ) + ghost;// Inner tile absolute 'y' val

  /** Compute outer tile variables **/
  int outer_dim    = BLOCK_DIM + 2 * ghost;             // Dimensions of the outer tile
  int outer_x      = inner_x - ghost;                   // Outer tile absolute 'x' val
  int outer_y      = inner_y - ghost;                   // Outer tile absolute 'y' val

  // Cache the vertical stencil values zm1, z, zp1
  float zlow;
  float zcenter;
  float zhigh;

  __shared__ float local_data_in [BLOCK_DIM+2][BLOCK_DIM+2];

  for( int k = 1; k < nz - 1; ++k ) { // Iterate frames along the Z axis
    int inner_offset = Index3D(nx, ny, inner_x, inner_y, k);
    int outer_offset = Index3D(nx, ny, outer_x, outer_y, k);
    float tijk = 0.01;

    local_data_in[ThreadIdx.x][ThreadIdx.y] =
        global_data_in[outer_offset + ThreadIdx.x + nx * ThreadIdx.y];

    if( (ThreadIdx.x < BLOCK_DIM) && (ThreadIdx.y < BLOCK_DIM) ){

      global_data_out[ inner_offset + ThreadIdx.x + ThreadIdx.y ] =
          // X direction
          (
              data_in[Index3D (nx, ny, i - 1, j, k)] +
              data_in[Index3D (nx, ny, i + 1, j, k)] -
              tijk
          ) / ( dx * dx )
          +
          // Y direction
          (
              data_in[Index3D (nx, ny, i, j - 1, k)] +
              data_in[Index3D (nx, ny, i, j + 1, k)] -
              tijk
          ) / ( dy * dy )
          +
          // Z direction
          (
              data_in[Index3D (nx, ny, i, j, k - 1)] +
              data_in[Index3D (nx, ny, i, j, k + 1)] -
              tijk
          ) / ( dz * dz );
    }
  }
}

template< class FieldT >
void __host__ divergence_float_gpu(FieldT* f, float dx, float dy, float dz){
    /** Determine dimensions
     *          - grab x, y extents ( we should get the entire memory window ) ?
     *          - decide how many 16x16 blocks are needed to tile the base plane.
     *            This might be a value that is not a multiple of 16.
     */
    using namespace SpatialOps::structured;
    cudaError err;

    MemoryWindow window = f->window_with_ghost();
    const IntVec& extent = window.extent(); // Get the dimensions of this window.
    unsigned int blockx = ( extent[0] / BLOCK_DIM + ( extent[0] % BLOCK_DIM == 0 ? 0 : 1) );
    unsigned int blocky = ( extent[1] / BLOCK_DIM + ( extent[1] % BLOCK_DIM == 0 ? 0 : 1) );

    float* h_workspace = (float*)malloc(f->get_data_size());
    float* d_workspace;
    cudaMalloc((void**)&d_workspace, f->get_data_size() );
    cudaMemcpy((void*)d_workspace, (void*)f->get_ext_pointer(), f->get_data_size(), cudaMemcpyDeviceToDevice );

    dim3 dimBlock( BLOCK_DIM + 2, BLOCK_DIM + 2 );
    dim3 dimGrid( blockx, blocky );

    std::cout << "Executing with " << blockx << " X " << blocky << " block grid" << std::endl;
    std::cout << "               " << BLOCK_DIM << " X " << BLOCK_DIM << " block threads\n";

    std::cout << "Extent: " << extent[0] << " " << extent[1] << " " << extent[2] << std::endl;
    if( cudaSuccess != ( err = cudaSetDevice( f->device_index() ) ) ){
      throw( std::runtime_error( cudaGetErrorString(err) ) );
    }

    clock_t tick = clock();
    for( int i = 0; i < 1e3; ++i ){
      _div_float_opt2<<<dimGrid, dimBlock, 0, 0>>>( d_workspace,
                                                  f->get_ext_pointer(),
                                                  dx, dy, dz,
                                                  extent[0], extent[1], extent[2]);

      cudaMemcpy(h_workspace, f->get_ext_pointer(), f->get_data_size(), cudaMemcpyDeviceToHost);
    }
    clock_t tock = clock() - tick;

    //cudaThreadSynchronize();
    cudaFree(d_workspace);


    std::cout << "Elapsed: " << (double)tock / (double)CLOCKS_PER_SEC << std::endl;

    //Dump first computed plane in the z axis
    std::ofstream file;
    file.open("cudaout.txt");
    for(int i = 0; i < extent[0]; ++i){
      for( int j = 0; j < extent[1]; ++j){
        float q = h_workspace[Index3D(extent[0], extent[1], i, j, 1)];
        file << std::setw(8) << q << " ";
      }
      file << std::endl;
    }
    file.close();
}

template< class FieldT >
void __host__ divergence_double(FieldT* f, double dx, double dy, double dz ){
 throw;
}

#endif /* CUDASTENCILLIB_H_ */
