nebo_cuda_prep_dir()

set( src
     BoxFilter.cpp
     StencilBuilder.cpp
    )

nebo_add_library( spatialops-stencil src )

if( ENABLE_TESTS )
  add_subdirectory( test )
endif( ENABLE_TESTS )

install( FILES
  BoxFilter.h
  FVStaggeredBCOp.h
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
