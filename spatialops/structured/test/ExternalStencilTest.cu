/*
 * ExternalStencilTest.cpp
 *
 *  Created on: Jan 5, 2012
 *      Author: Devin Robison
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include "CudaStencilLib.h"
#include "ExternalStencilTest.h"

using namespace SpatialOps;
using namespace structured;
using namespace Point;
using namespace ema::cuda;

#define SMALL 16
//#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))
//#define Index3DG(_nx, _ny, _g, _i, _j, _k) (((_i)+(_g)) + _nx*(((_j)+(_g) + _ny*((_k)+(_g)) ))

void stencil7_cpu( float* data_in, float* data_out, int nx, int ny, int nz ){
  int dx = 1, dy = 1, dz = 1;
  for( int i = 1; i < SMALL-1; ++i)
    for( int j = 1; j < SMALL-1; ++j)
      for( int k = 1; k < SMALL-1; ++k){
        float tijk = 2*data_in[ Index3D( nx, ny, i,j,k) ];
        data_out[Index3D (nx, ny, i, j, k)] =
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

int main(int argc, char** argv){
  const IntVec sset(SMALL,SMALL,SMALL);
  const IntVec mset(64,64,64);
  const IntVec lset(128,128,128);

  float* svals = new float[SMALL*SMALL*SMALL];
  float* data_out = new float[SMALL*SMALL*SMALL];

  srand(0);

  for( int i = 1; i < SMALL-1; ++i)
    for( int j = 1; j < SMALL-1; ++j)
      for( int k = 1; k < SMALL-1; ++k)
         svals[Index3D(SMALL,SMALL, i,j,k)] = (float)rand() / (float)RAND_MAX;

  CUDADeviceInterface& CDI = CUDADeviceInterface::self();

  MemoryWindow swin(sset);
  MemoryWindow mwin(mset);
  MemoryWindow lwin(lset);


  PointFloatField pff(swin, NULL, InternalStorage, LOCAL_RAM, 0);
  PointFloatField spff(swin, svals, InternalStorage, EXTERNAL_CUDA_GPU, 0);

  CDI.print_device_info();

  clock_t tick = clock();
  for( int i = 0; i < 1e3; ++i )
  stencil7_cpu(svals, data_out, SMALL, SMALL, SMALL);
  clock_t tock = clock() - tick;

  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++)
      std::cout << Index3DG(4,4,1,j,i,0) << " ";
    std::cout << std::endl;
  }

  std::cout << "Elapsed: " << (double)tock / (double)CLOCKS_PER_SEC << std::endl;

  //Dump first computed plane in the z axis
  std::ofstream file;
  file.open("cpuout.txt");
  for(int i = 0; i < SMALL; ++i){
    for( int j = 0; j < SMALL; ++j){
      float q = data_out[Index3D(SMALL, SMALL, i, j, 1)];
      file << std::setw(8) << q << " ";
    }
    file << std::endl;
  }
  file.close();

  divergence_float_gpu(&spff, 1, 1, 1);
  pff = spff;

  delete[] svals;
}
