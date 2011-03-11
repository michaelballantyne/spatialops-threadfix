get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( src
  ${dir}/LagrangePoly.cpp
  ${dir}/FVOneDimensional.cpp
  )

if( ENABLE_TESTS )
  add_library( spatialops-onedim ${src} )
  add_subdirectory( test )
endif( ENABLE_TESTS )

install( FILES
  FVStaggered1DOperatorTypes.h
  DESTINATION include/spatialops/structured/
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ 
  )

if( ENABLE_MATRIX )
  set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )

  install( FILES
    LagrangePoly.h
    FVOneDimensional.h
    DESTINATION include/spatialops
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )

endif( ENABLE_MATRIX )
