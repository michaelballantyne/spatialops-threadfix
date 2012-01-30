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

set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )

if( ENABLE_CUDA )
  if( CUDA_CUDA_LIBRARY )
		include_directories( ${CUDA_INCLUDE_DIRS} )
		set( cuda_src
		  ${dir}/CudaMemoryAllocator.cu
		)
		cuda_add_library( spatialops_cuda STATIC ${cuda_src} )
		target_link_libraries( spatialops_cuda ${CUDA_LIBRARIES} ${Boost_LIBRARIES} )
		    
		install( TARGETS spatialops_cuda EXPORT so_cuda DESTINATION lib/spatialops )
		install( EXPORT so_cuda DESTINATION lib/spatialops )
		      
		#Configure CMM 
		configure_file( 
		  ${PROJECT_SOURCE_DIR}/spatialops/structured/so_cudaConfig.cmake.in
		  ${PROJECT_BINARY_DIR}/config/so_cudaConfig.cmake
		)   
		    
		install( FILES
		  ${PROJECT_BINARY_DIR}/config/so_cudaConfig.cmake
		    DESTINATION ${SpatialOps_CONFIG_INSTALL}
		)
	endif( CUDA_CUDA_LIBRARY )
endif( ENABLE_CUDA )

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
