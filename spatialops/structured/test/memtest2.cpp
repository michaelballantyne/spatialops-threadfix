#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

#define DEBUG_CUDA_VERBOSE
#include <spatialops/structured/ExternalAllocators.h>

using namespace ema::cuda;

int main(int argc, char** argv) {
  CUDADeviceInterface& CDI = CUDADeviceInterface::self();
  CUDASharedPointer p;
  int temp;
  bool should_fail = false;

  try {
    for (int device = 0; device < CDI.get_device_count(); ++device) {
      CDI.update_memory_statistics();
      const CUDAMemStats& CMS = CDI.get_memory_statistics(device);

      std::cout << "Testing pointer creation, allocating " << CMS.f
          << " bytes, on device " << device << "...";
      p = CDI.get((CMS.f / 2), device);

      if (p.get_refcount() != 1) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found invalid reference count on pointer P for device: "
            << device;
        throw(std::runtime_error(msg.str()));
      }

      if (p.get_deviceID() != device) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found incorrect deviceID on pointer P for device: " << device;
        throw(std::runtime_error(msg.str()));
      }
      std::cout << " PASS\n";

      std::cout << "Testing pointer assignment operations...";
      CUDASharedPointer q;
      q = p;

      if (q.get_refcount() != 2) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found invalid reference count, " << q.get_refcount()
            << ", on pointer Q for device: " << device;
        throw(std::runtime_error(msg.str()));
      }

      if (q.get_deviceID() != device) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found incorrect deviceID on pointer Q for device: " << device;
        throw(std::runtime_error(msg.str()));
      }
      std::cout << " PASS\n";

      std::cout << "Testing pointer detach and reassignment...";
      q.detatch();

      if (p.get_refcount() != 1) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found invalid reference count, " << p.get_refcount()
            << ", on pointer P for device: " << device;
        throw(std::runtime_error(msg.str()));
      }

      if (p.get_deviceID() != device) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found incorrect deviceID on pointer P for device: " << device;
        throw(std::runtime_error(msg.str()));
      }

      for (int i = 0; i < 10; i++) {
        CUDASharedPointer x;
        x = p;
        CUDASharedPointer y = x;
        CUDASharedPointer z(y);
      }

      if (p.get_refcount() != 1) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found invalid reference count," << p.get_refcount()
            << ", on pointer P for device: " << device;
        throw(std::runtime_error(msg.str()));
      }

      if (p.get_deviceID() != device) {
        std::cout << " FAIL\n";
        std::ostringstream msg;
        msg << "Found incorrect deviceID on pointer P for device: " << device;
        throw(std::runtime_error(msg.str()));
      }
      std::cout << " PASS\n";
    }
  }
  catch( std::runtime_error e) {
    std::cout << e.what() << std::endl;
    exit(1);
  }

  printf("Success\n");

  return 0;
}
