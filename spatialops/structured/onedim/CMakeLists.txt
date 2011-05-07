get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( src
  ${dir}/LagrangePoly.cpp
  )

if( ENABLE_MATRIX )
  set( src ${src} ${dir}/FVOneDimensional.cpp )
endif( ENABLE_MATRIX )

if( ENABLE_TESTS )
  add_library( spatialops-onedim ${src} )
  add_subdirectory( test )
endif( ENABLE_TESTS )

if( ENABLE_MATRIX )
  install( FILES
    FVStaggered1DOperatorTypes.h
    DESTINATION include/spatialops/structured/
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ 
    )

  install( FILES
    FVOneDimensional.h
    DESTINATION include/spatialops
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )

endif( ENABLE_MATRIX )

install( FILES
  LagrangePoly.h
  DESTINATION include/spatialops
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )
