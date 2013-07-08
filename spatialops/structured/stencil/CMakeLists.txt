nebo_cuda_prep_dir()

set( src
     BoxFilter.cpp
     StencilBuilder.cpp
    )

nebo_set_add_file( src FDStencil2.cpp )
nebo_set_add_file( src NullStencil.cpp )
nebo_set_add_file( src Stencil2.cpp )
nebo_set_add_file( src Stencil4.cpp )

if( ENABLE_CUDA )
    set( src ${src}
         CudaStencil2.cu
         CudaStencil4.cu
        )
endif( ENABLE_CUDA )

nebo_add_library( spatialops-stencil src )

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
  Stencil4_def.h
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
