get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( src
  ${dir}/NullStencil.cpp
  ${dir}/Stencil2.cpp
  ${dir}/Stencil4.cpp
  ${dir}/StencilBuilder.cpp
  )

if( ENABLE_TESTS )
  add_library(
    spatialops-stencil
    STATIC
    ${src}
    )
  add_subdirectory( test )
endif( ENABLE_TESTS )


# only install if enabled
if( ENABLE_STENCIL )
  set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )

  install( FILES
    FVStaggeredBCOp.h
    FVStaggeredOperatorTypes.h
    NullStencil.h
    Stencil2.h
    Stencil4.h
    StencilBuilder.h
    DESTINATION include/spatialops/structured/stencil/
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )

  install( FILES
    FVStaggered.h
    DESTINATION include/spatialops/structured/
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )

endif( ENABLE_STENCIL )