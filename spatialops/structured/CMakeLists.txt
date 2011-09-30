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

if( ENABLE_MATRIX )
  add_subdirectory( matrix )
endif( ENABLE_MATRIX )

if( ENABLE_STENCIL )
  add_subdirectory( stencil )
endif( ENABLE_STENCIL )

set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )

install( FILES
  FVStaggeredBCTools.h
  FVStaggeredFieldTypes.h
  FVTools.h
  Grid.h
  IntVec.h
  MemoryWindow.h
  SpatialField.h
  SpatialFieldStore.h
  DESTINATION include/spatialops/structured
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )