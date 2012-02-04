get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( src
  ${dir}/BoxFilter.cpp
  ${dir}/NullStencil.cpp
  ${dir}/Stencil2.cpp
  ${dir}/Stencil4.cpp
  ${dir}/StencilBuilder.cpp
  )

##Stencil2.(cpp/cu) switch
if( ENABLE_CUDA )
    set(src ${src} ${dir}/CudaStencil2.cu)
    cuda_add_library(
        spatialops-stencil
        STATIC
        ${src}
    )
else( ENABLE_CUDA )
    add_library(
        spatialops-stencil
        STATIC
        ${src}
    )
endif( ENABLE_CUDA )

if( ENABLE_TESTS )
  add_subdirectory( test )
endif( ENABLE_TESTS )

set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )

install( FILES
  BoxFilter.h
  FVStaggeredBCOp.h
  NullStencil.h
  CudaStencil2.h
  CudaStencil2Bridge.h
  Stencil2.h
  Stencil4.h
  StencilBuilder.h
  DESTINATION include/spatialops/structured/stencil/
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

install( FILES
  FVStaggered.h
  FVStaggeredOperatorTypes.h
  DESTINATION include/spatialops/structured/
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

install( TARGETS
  spatialops-stencil
  EXPORT SpatialOps
  DESTINATION lib/spatialops/
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )
