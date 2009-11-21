cmake_minimum_required(VERSION 2.6)
project( SpatialOps CXX C Fortran )

set( MAJOR_VERSION "0" )
set( MINOR_VERSION "1" )
set( PATCH_VERSION "0" )

#=====================================================================
#---------------------------------------------------------------------
# BEGIN USER MODIFIABLE BLOCK


option( LINALG_TRILINOS "use Trilinos for linear algebra - also requires definition of {Trilinos_path}"  OFF )
option( LINALG_UBLAS    "use UBLAS for linear algebra - cannot solve linear systems"  ON )


set( CMAKE_BUILD_TYPE  # choose one type below.
  Release
# Debug
# RelWithDebInfo
CACHE STRING "set build type"
  )


# Set the path for the Daixtrose installation on your machine
# daixtrose is used for expression templates - we support only a few
# thus far.  If empty then SpatialOps will be built without Daixtrose
# support
set( Daixtrose
  ""
  CACHE PATH "path to daixtrose installation (optional)"
  )

# Set the location for the installation.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX
    ${CMAKE_CURRENT_BINARY_DIR}
    CACHE PATH "Installation directory" FORCE
    )
endif()


set( SpatialOps_CONFIG_INSTALL
  ${CMAKE_INSTALL_PREFIX}/share
  CACHE PATH
  "Location to install configuration files"
  )

# END USER MODIFIABLE BLOCK
#---------------------------------------------------------------------
#=====================================================================


if( LINALG_TRILINOS )
# Set the path for the trilinos installation on your machine.  If
# empty, SpatialOps will be built without trilinos (using UBLAS
# instead) comment this out to build without trilinos
set( Trilinos_path
  ""
  CACHE PATH "path to trilinos intallation"
  )
endif()


#--------------------------------------
# resolve required libraries

# setup the boost library
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREAD ON)
set(Boost_ADDITIONAL_VERSIONS "1.39" "1.39.0" )
find_package( Boost REQUIRED 
  COMPONENTS  # don't require any compiled boost libs - only headers.
  NO_MODULE
  ) 
set( TPL_INCLUDE_DIRS ${TPL_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} )
set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Boost_LIBRARIES} )

if( LINALG_TRILINOS )
  message( STATUS "Using Trilinos for linear algebra" )
  set( Trilinos_INCLUDE_DIR ${Trilinos_path}/include )
  set( Trilinos_LIBRARY_DIR ${Trilinos_path}/lib )
  set( Trilinos_libs
    ${Trilinos_LIBRARY_DIR}/libaztecoo.a
    ${Trilinos_LIBRARY_DIR}/libepetraext.a
    ${Trilinos_LIBRARY_DIR}/libepetra.a
    ${Trilinos_LIBRARY_DIR}/libteuchos.a
    )
  set( TPL_INCLUDE_DIRS ${TPL_INCLUDE_DIRS} ${Trilinos_INCLUDE_DIR} )
  set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Trilinos_libs} )

  # find BLAS and LAPACK (required for trilinos)
  find_package( BLAS REQUIRED )
  find_package( LAPACK REQUIRED )
  set( TPL_LIBRARIES ${TPL_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    )

  set( LINALG_TRILINOS ON )
  set( LINALG_UBLAS OFF )
else( Trilinos_path )
  message( STATUS "Using UBLAS for linear algebra" )
  set( LINALG_UBLAS ON )
  set( LINALG_TRILINOS OFF )
endif()

# until we get a package finder for daixtrose, we will need this:
if( Daixtrose )
  set( Daixtrose_INCLUDE_DIR ${Daixtrose}/include )
  set( TPL_INCLUDE_DIRS ${TPL_INCLUDE_DIRS} ${Daixtrose_INCLUDE_DIR} )
  set( HAVE_DAIXTROSE ON )
endif()
#--------------------------------------


include_directories(
  include
  ${TPL_INCLUDE_DIRS}
  )


# build the configure file using the information defined here.
configure_file(
  ${PROJECT_SOURCE_DIR}/SpatialOpsConfigure.h.in
  ${PROJECT_BINARY_DIR}/config/spatialops/SpatialOpsConfigure.h
  )
include_directories( ${PROJECT_BINARY_DIR}/config )


# setup the build for the spatialops library.  Note that the linker
# order is important, and future revisions to libraries such as
# trilinos may require a modification of link order.
set( SO_Sources
  src/LagrangePoly.cpp
  src/FVOneDimensional.cpp
  src/FVTools.cpp
  include/spatialops/particles/ParticleFieldTypes.cpp
  )
if( LINALG_TRILINOS )
  set( SO_Sources ${SO_Sources}
    src/LinAlgTrilinos.cpp
    src/LinearSystem.cpp
    )
endif()
if( LINALG_UBLAS )
  set( SO_Sources ${SO_Sources} src/LinAlgUBlas.cpp )
endif()


if( CMAKE_BUILD_TYPE STREQUAL Release )
  add_definitions( -O3 -fexpensive-optimizations )
  message( STATUS "turning on aggressive optimizations" )
endif( CMAKE_BUILD_TYPE STREQUAL Release )

set( SpatialOps_LIBRARY spatialops )
add_library( ${SpatialOps_LIBRARY}
  STATIC
  ${SO_Sources}
  )
target_link_libraries( ${SpatialOps_LIBRARY} ${TPL_LIBRARIES} )


# this will compile all of the test cases.  Use "make test" to execute
# them and see a summary.
enable_testing()
add_subdirectory( test )



#--- INSTALLATION SECTION ---
# "make install" performs these commands, along with any install commands inside child CMakeLists.txt files.
message( STATUS "Files will be installed to: " ${CMAKE_INSTALL_PREFIX} )

# install spatialops headers
file( GLOB so_headers "include/spatialops/*.h" )
file( GLOB so_structured_headers "include/spatialops/structured/*.h" )
file( GLOB so_particle_headers "include/spatialops/particles/*.h" )
set( so_headers
  ${so_headers}
  ${PROJECT_BINARY_DIR}/config/spatialops/SpatialOpsConfigure.h
  )
install( FILES
  ${so_headers}
  DESTINATION include/spatialops
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )
install( FILES
  ${so_structured_headers}
  DESTINATION include/spatialops/structured
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )
install( FILES
  ${so_particle_headers}
  DESTINATION include/spatialops/particles
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

# install the spatialops library
install( TARGETS
  ${SpatialOps_LIBRARY}
  ARCHIVE
  DESTINATION lib/spatialops/
  )



# ------------------------------------------------------------------
# configuration for projects using SpatialOps this section defines
# variables that are used for building the configuration file for the
# SpatialOps installation.  This can be used by other CMake build
# systems that use the SpatialOps library.

# get the location for the installed library
get_target_property( SOLocation ${SpatialOps_LIBRARY} LOCATION )

# get the full name for the installed library
get_filename_component( SOInstallLib ${SOLocation} NAME )

set( SpatialOps_INCLUDE ${CMAKE_INSTALL_PREFIX}/include )

# at this point, all of the required variables are defined to generate
# the configure file
configure_file(
  ${PROJECT_SOURCE_DIR}/SpatialOpsConfig.cmake.in
  ${PROJECT_BINARY_DIR}/config/SpatialOpsConfig.cmake
  @ONLY
  )

install( TARGETS ${SpatialOps_LIBRARY} EXPORT SpatialOps DESTINATION lib/spatialops )
install( EXPORT SpatialOps DESTINATION lib/spatialops )

# install all of these files into the spatial ops lib dir
install( FILES
  ${PROJECT_BINARY_DIR}/config/SpatialOpsConfig.cmake
   DESTINATION ${SpatialOps_CONFIG_INSTALL}
   )
#-----------------------------------------------------------------


# Information for creating packages
#set( CPACK_PACKAGE_VENDOR "James C. Sutherland (james.sutherland@utah.edu" )
#set( CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION} )
#set( CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION} )
#set( CPACK_PACKAGE_VERSION_PATCH ${PATCH_VERSION} )
#include( CPack )
