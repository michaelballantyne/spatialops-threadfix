get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( src
  ${dir}/Grid.cpp
  ${dir}/MemoryWindow.cpp
  )

if( ENABLE_TESTS )
  add_library(
    spatialops-structured
    STATIC
    ${src}
    )
  add_subdirectory( test )
endif( ENABLE_TESTS )

add_subdirectory( onedim )
add_subdirectory( stencil )

if( ENABLE_CUDA )
    include_directories( ${CUDA_INCLUDE_DIRS} )
    set( TPL_LIBRARIES ${TPL_LIBRARIES} ${CUDA_LIBRARIES} )
    set( src ${src}
      ${dir}/CudaMemoryAllocator.cu
    )
endif( ENABLE_CUDA )

set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )

install( FILES
  FVStaggeredBCTools.h
  FVStaggeredFieldTypes.h
  FVTools.h
  Grid.h
  IndexTriplet.h
  IntVec.h
  ExternalAllocators.h
  CudaMemoryAllocator.h
  MemoryTypes.h
  MemoryWindow.h
  SpatialField.h
  SpatialFieldStore.h
  DESTINATION include/spatialops/structured
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )
