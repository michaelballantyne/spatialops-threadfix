add_subdirectory( structured )
add_subdirectory( particles  )
add_subdirectory( pointfield )


# setup the build for the spatialops library.  Note that the linker
# order is important, and future revisions to libraries such as
# trilinos may require a modification of link order.

if( LINALG_TRILINOS )
  set( SO_Sources ${SO_Sources}
    LinAlgTrilinos.cpp
    LinearSystem.cpp
    )
endif()
if( LINALG_UBLAS )
  set( SO_Sources ${SO_Sources} LinAlgUBlas.cpp )
endif()

# add in collected dependencies from subdirs
set( SO_Sources ${SO_Sources} ${SO_SRC} )

add_library( spatialops
  STATIC
  ${SO_Sources}
  )
target_link_libraries( spatialops ${TPL_LIBRARIES} )

set( SO_HEADERS
  FieldExpressions.h
  FieldExpressionsExtended.h
  FieldReductions.h
  FieldFunctions.h
  OperatorDatabase.h
  SpatialOpsDefs.h
  SpatialOpsTools.h
  WriteMatlab.h
  )

if( ENABLE_MATRIX )
  set( SO_HEADERS ${SO_HEADERS}
    LinearSystem.h
    SpatialOperator.h
    )
  if( LINALG_TRILINOS )
    set( SO_HEADERS ${SO_HEADERS} LinAlgTrilinos.h )
  endif()
  if( LINALG_UBLAS )
    set( SO_HEADERS ${SO_HEADERS} LinAlgUBlas.h )
  endif()
endif( ENABLE_MATRIX )

install( FILES
  ${SO_HEADERS}
  DESTINATION include/spatialops
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

install( TARGETS
  spatialops
  EXPORT SpatialOps
  DESTINATION lib/spatialops/
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )
install( EXPORT SpatialOps DESTINATION lib/spatialops )