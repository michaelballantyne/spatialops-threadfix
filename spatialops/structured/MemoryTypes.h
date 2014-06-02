/*
 * DeviceTypes.h
 *
 *  Created on: Dec 15, 2011
 *      Author: Devin Robison
 */

#ifndef MEMORYTYPES_H_
#define MEMORYTYPES_H_
#include <string>
namespace SpatialOps{

#define CPU_INDEX -1
#define GPU_INDEX  0

#define IS_CPU_INDEX(INDEX)   (INDEX == CPU_INDEX)
#define IS_GPU_INDEX(INDEX)   (INDEX >= GPU_INDEX)
#define IS_VALID_INDEX(INDEX) (INDEX >= CPU_INDEX)

  /**
   *  \class DeviceTypeTools
   *  \brief Provide descriptions
   */
  namespace DeviceTypeTools {
  inline std::string get_memory_type_description( short int deviceIndex )
  {
    if( deviceIndex == CPU_INDEX ){
       return std::string("(Locally allocated, generic system RAM)");
     }
     else if( IS_GPU_INDEX(deviceIndex) ){
       return std::string("(Externally allocated, CUDA GPU device)");
     }
     else{
       return std::string("(Unknown or Invalid)");
     }
  }

  } // namespace DeviceTypeTools
} // namespace SpatialOps



#endif /* DEVICETYPES_H_ */
