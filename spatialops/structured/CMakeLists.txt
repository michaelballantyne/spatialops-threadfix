set( src
     Grid.cpp
     MemoryWindow.cpp
    )

if( ENABLE_CUDA )
    set( src ${src}
         CudaMemoryAllocator.cu
         GPUtest.cu
        )
    cuda_add_library( spatialops-structured STATIC ${src} )
else( ENABLE_CUDA )
    add_library( spatialops-structured STATIC ${src} )
endif( ENABLE_CUDA )

if( ENABLE_TESTS )
    add_subdirectory( test )
endif( ENABLE_TESTS )

add_subdirectory( onedim )
target_link_libraries( spatialops-structured spatialops-onedim )

add_subdirectory( stencil )
target_link_libraries( spatialops-structured spatialops-stencil )

install( FILES
  FVStaggeredBCTools.h
  FVStaggeredFieldTypes.h
  FVTools.h
  Grid.h
  GhostData.h
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

install( TARGETS spatialops-structured
         EXPORT SpatialOps
         DESTINATION lib/spatialops/
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
        )
