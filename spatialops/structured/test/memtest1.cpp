#include <iostream>
#include <string>
#include <stdexcept>
#include <string.h>
#include <stdlib.h>

#define DEBUG_CUDA_VERBOSE
#include <spatialops/structured/ExternalAllocators.h>

using namespace ema::cuda;

int main(int argc, char** argv) {
  CUDADeviceInterface& CDI = CUDADeviceInterface::self();
  CUDASharedPointer p;
  int temp;

  for (int device = 0; device < CDI.get_device_count(); ++device) {
    try {
      for (int i = 1;; i *= 2) {
        CDI.update_memory_statistics();
        const CUDAMemStats& CMS = CDI.get_memory_statistics(device);

        std::cout << "Attempting allocation of size " << i << " bytes, on device " << device;
        std::cout << "\n\t Free memory: " << CMS.f << " / " << CMS.t << std::endl;

        p = CDI.get(i, 0);
        char* bytesin = (char*)malloc(sizeof(char)*i);
        char* bytesout = (char*)malloc(sizeof(char)*i);

        bzero(bytesin, i);
        bzero(bytesout, i);
        std::cout << "Checking zeros read/write... ";
        CDI.copy_to(p, (void*)bytesin, i);
        CDI.copy_from((void*)bytesout, p, i);

        if( memcmp(bytesin, bytesout, i) ){
           std::cout << "failed -> Zero byte pattern does do not match\n";
           exit(1);
        }

        std::cout << "OK\n";

        memset(bytesin, 1, i);
        bzero(bytesout, i);
        std::cout << "Checking ones read/write... ";
        CDI.copy_to(p, (void*)bytesin, i);
        CDI.copy_from((void*)bytesout, p, i);

        if( memcmp(bytesin, bytesout, i) ){
           std::cout << "failed -> Ones byte pattern does not match\n";
           exit(1);
        }

        std::cout << "OK\n";

        srand(0);
        for(int k = 0; i < i; ++k){
          bytesin[k] = rand();
        }
        bzero(bytesout, i);
        std::cout << "Checking random read/write... ";
        CDI.copy_to(p, (void*)bytesin, i);
        CDI.copy_from((void*)bytesout, p, i);

        if( memcmp(bytesin, bytesout, i) ){
           std::cout << "failed -> Random byte pattern does not match\n";
           exit(1);
        }

        std::cout << "OK\n";
      }
    }
    catch ( std::runtime_error e ) {
      //Note: malloc will fail at some point, this is expected
      std::cout << e.what() << std::endl;
    }
  }

  printf("Success\n");

  return 0;
}
