# import the build settings from SpatialOps
include( CMakeImportBuildSettings )
cmake_import_build_settings( ${SpatialOps_BUILD_SETTINGS_FILE} )

# advertise location of header files
include_directories(
  ${SpatialOps_INCLUDE_DIR}
  ${SpatialOps_TPL_INCLUDE_DIRS}
  )

set( SpatialOps_LIBRARIES ${SpatialOps_LIBRARIES} ${spatialops_LIB_DEPENDS} )

message( STATUS
  "SpatialOps include paths: "
  ${SpatialOps_INCLUDE_DIR}
  ${SpatialOps_TPL_INCLUDE_DIRS}
  )