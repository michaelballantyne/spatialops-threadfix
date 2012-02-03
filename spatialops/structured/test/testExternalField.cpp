/*
 * testExternalField.cpp
 *
 *  Created on: Dec 21, 2011
 *      Author: Devin Robison
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "testExternalField.h"

using namespace SpatialOps;
using namespace structured;
using namespace Point;
using namespace ema::cuda;

int main(int argc, char** argv) {
  unsigned int bytes = 128 * 128 * 128;
  double* T1 = (double*) malloc(sizeof(double) * (bytes));
  double* T2 = (double*) malloc(sizeof(double) * (bytes));
  const IntVec npts(128, 128, 128);
  const IntVec badpts(1024, 1024, 1024);

  CUDADeviceInterface& CDI = CUDADeviceInterface::self();

  const MemoryWindow window(npts);

  /** Excercise constructors **/
  memset(T1, 1, bytes * sizeof(double));
  try {
    std::cout << "Checking simple allocation, " << bytes * sizeof(double)
        << " bytes: ";
    PointField p(window, T1, InternalStorage, EXTERNAL_CUDA_GPU, 0);
    std::cout << "PASS\n";
  } catch (...) {
    std::cout << "FAIL\n";
    return -1;
  }

  try {
    std::cout << "Checking proper value initialization: ";
    PointField p(window, T1, InternalStorage, EXTERNAL_CUDA_GPU, 0);
    CDI.memcpy_from((void*) T2, p.ext_field_values(), bytes * sizeof(double),
        p.device_index());
    for (unsigned int k = 0; k < bytes; ++k) {
      if (T2[k] != T1[k]) {
        throw;
      }
    }
    std::cout << "PASS\n";
  } catch (...) {
    return -1;
  }

  try {
    std::cout
        << "Attempting to construct field with invalid memory count ( This shouldn't work ): ";
    // try to allocate more points then we have available memory for.
    CDI.update_memory_statistics();
    CUDAMemStats CMS;
    CDI.get_memory_statistics(CMS, 0);
    const IntVec badpts(1, 1, CMS.f / sizeof(double) + 1);
    const MemoryWindow badwindow(badpts);

    PointField fail(badwindow, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
    std::cout << "FAIL\n";
    return -1;
  } catch (std::runtime_error e) {
    std::cout << "PASS\n";
  }

  try {
    std::cout
        << "Attempting to construct field with invalid field memory type ( this shouldn't work ): ";
    PointField fail(npts, T1, InternalStorage, DEBUG_TEST_OPT, 0);
    std::cout << "FAIL\n";
    return -1;
  } catch (std::runtime_error e) {
    std::cout << "PASS\n";
  }

  /** Check assignment and comparison operations **/

  try {
    PointField p(window, T1, InternalStorage, EXTERNAL_CUDA_GPU, 0);
    PointField q(window, T1, InternalStorage, LOCAL_RAM, 0);

    memset(T2, rand(), bytes * sizeof(double));
    PointField r(window, T2, InternalStorage, EXTERNAL_CUDA_GPU, 0);

    std::cout << "P, Q, R assignment and comparison testing: ";
    if (!(p == q)) {
      throw;
    }

    if ((p == q) != !(p != q)) {
      throw;
    }

    if (p == r) {
      throw;
    }

    if (q == r) {
      throw;
    }

    q = r;

    if (q != r) {
      throw;
    }

    if ((p == q) || !(p != q)) {
      throw;
    }

    std::cout << "PASS\n";
  } catch (std::runtime_error e) {
    std::cout << "FAIL\n";
    return -1;
  }

  try {
    std::cout << "Checking SpatialFieldPointer creation: ";
    PointField q(window, T1, InternalStorage, LOCAL_RAM, 0);
    //SpatFldPtr<PointField> p = SpatialFieldStore<PointField>::self().get(window);
    //SpatFldPtr<PointField> r = SpatialFieldStore<PointField>::self().get(q);
    SpatFldPtr<PointField> t = SpatialFieldStore<PointField>::self().get(window, EXTERNAL_CUDA_GPU, 0);

    //p.detach();
    //r.detach();
    t.detach();
    std::cout << "PASS\n";
  } catch ( std::runtime_error e) {
    std::cout << "FAIL\n";
    std::cout << e.what() << std::endl;
    return -1;
  }

  try {
    std::cout << "Checking proper failure conditions ( bad device type ): ";
    SpatFldPtr<PointField> t = SpatialFieldStore<PointField>::self().get(window, DEBUG_TEST_OPT, 0);
    std::cout << "FAIL\n";
    return -1;
  } catch ( std::runtime_error e){
    std::cout << "PASS\n";
  }

  try {
    std::cout << "Checking proper failure conditions ( bad device index ): ";
    SpatFldPtr<PointField> t = SpatialFieldStore<PointField>::self().get(window, EXTERNAL_CUDA_GPU, 75);
    std::cout << "FAIL\n";
    return -1;
  } catch ( std::runtime_error e){
    std::cout << "PASS\n";
  }
  std::cout << std::endl << "All tests passed: Success\n";
}
