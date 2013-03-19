set( src
     BoxFilter.cpp
     FDStencil2.cpp
     NullStencil.cpp
     Stencil2.cpp
     Stencil4.cpp
     StencilBuilder.cpp
    )

if( ENABLE_CUDA )
    set( src ${src}
         CudaStencil2.cu
         CudaStencil4.cu
        )
    cuda_add_library( spatialops-stencil STATIC ${src} )
else( ENABLE_CUDA )
    add_library( spatialops-stencil STATIC ${src} )
endif( ENABLE_CUDA )

if( ENABLE_TESTS )
  add_subdirectory( test )
endif( ENABLE_TESTS )

install( FILES
  BoxFilter.h
  FDStencil2.h
  FVStaggeredBCOp.h
  NullStencil.h
  CudaStencil2.h
  CudaStencil2Bridge.h
  CudaStencil4.h
  CudaStencil4Bridge.h
  Stencil2.h
  Stencil2_def.h
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

install( TARGETS spatialops-stencil
         EXPORT SpatialOps
         DESTINATION lib/spatialops/
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
        )
