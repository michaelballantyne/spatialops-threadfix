add_subdirectory( structured )
add_subdirectory( particles  )
add_subdirectory( pointfield )

# add in collected dependencies from subdirs
set( SO_Sources ${SO_Sources} ${SO_SRC} )

set( SO_HEADERS
  FieldExpressions.h
  FieldReductions.h
  FieldFunctions.h
  OperatorDatabase.h
  SpatialOpsDefs.h
  SpatialOpsTools.h
  WriteMatlab.h
  FieldExpressionsExtended.h
  )
  
if ( ENABLE_THREADS )
   set( SO_Sources ${SO_Sources} ThreadPool.cpp )
   set( SO_HEADERS ${SO_HEADERS} ThreadPool.h )
   add_subdirectory( threadpool )
endif( ENABLE_THREADS )

if( ENABLE_CUDA )
  cuda_add_library( spatialops
  STATIC
  ${SO_Sources}
  )
else( ENABLE_CUDA )
  add_library( spatialops
  STATIC
  ${SO_Sources}
  )
endif( ENABLE_CUDA )


target_link_libraries( spatialops ${TPL_LIBRARIES} )

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
