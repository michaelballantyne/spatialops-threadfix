get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( src
  ${dir}/LagrangePoly.cpp
  )

if( ENABLE_TESTS )
  add_library( spatialops-onedim ${src} )
  add_subdirectory( test )
endif( ENABLE_TESTS )

install( FILES
  LagrangePoly.h
  DESTINATION include/spatialops
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

set( SO_SRC ${SO_SRC} ${src} PARENT_SCOPE )
