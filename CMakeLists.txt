cmake_minimum_required(VERSION 2.6)
project( SpatialOps CXX C )

set( MAJOR_VERSION "1" )
set( MINOR_VERSION "3" )
set( PATCH_VERSION "0" )


option( ENABLE_BOOST_SERIALIZATION "allow serialization of information via the boost library" OFF )
option( ENABLE_TESTS "Enable building of regression tests (recommended)" ON )
option( ENABLE_THREADS "Enable multi-threaded support for stencil operations" OFF )
option( ENABLE_CUDA   "Build Spatial Ops with CUDA support" OFF )

set( NTHREADS 1 CACHE STRING "Number of threads to use if ENABLE_THREADS is ON" )

if( ENABLE_THREADS )
  set( TILL_THREADS ON )
  set( FIELD_EXPRESSION_THREADS ON )
  set( STENCIL_THREADS ON )
else()
  set( NTHREADS 1 )  
endif()

# default to release builds
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )

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



#--------------------------------------
# resolve required libraries

# setup the boost library
set( boost_comps "program_options" )

if( ENABLE_THREADS )
   set( boost_comps ${boost_comps} "thread" )
endif( ENABLE_THREADS )

set(Boost_USE_MULTITHREAD ON)
if( ENABLE_BOOST_SERIALIZATION )
  set( SOPS_BOOST_SERIALIZATION ON )
  set( boost_comps ${boost_comps} "serialization" )
else()
  set( SOPS_BOOST_SERIALIZATION OFF )
endif()

set(BOOST_MIN_VERSION "1.44.0" )
set(Boost_USE_STATIC_LIBS ON)
find_package( Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${boost_comps} )
if( NOT Boost_FOUND )
  message( FATAL_ERROR "Boost was not found.  Try setting BOOST_ROOT to the path where boost is installed and try again" )
else()
  message( STATUS "Boost header installation found at: " ${Boost_INCLUDE_DIR} )
endif()
set( TPL_INCLUDE_DIRS ${TPL_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  )
if( boost_comps )
  set( TPL_LIBRARIES ${TPL_LIBRARIES} ${Boost_LIBRARIES} )
  message( STATUS ${Boost_LIBRARIES} )
endif( boost_comps )

#-- CUDA
if( ENABLE_CUDA )
  message( STATUS "SpatialOps will be compiled with CUDA Support" )
  FIND_PACKAGE( CUDA REQUIRED )
  SET( CUDA_ARCHITECTURE "1.0" CACHE STRING "1.0" )
  STRING( REGEX REPLACE "([1-9])\\.([0-9])" "\\1\\2" CUDA_ARCH_STR ${CUDA_ARCHITECTURE} )   
  SET( CUDA_NVCC_FLAGS "-arch=sm_${CUDA_ARCH_STR}" CACHE STRING "Compiler flags passed to NVCC" FORCE )
  SET( CUDA_NVCC_FLAGS_RELEASE "--optimize=03" CACHE STRING "Release build compiler flags" FORCE )
else( ENABLE_CUDA )
#avoids confusion if someone tries to build with cuda support enabled and then reverts without clearing
#the build directory
	SET( CUDA_CUDA_LIBRARY OFF )
endif( ENABLE_CUDA )

#--------------------------------------


include_directories(
  ${PROJECT_SOURCE_DIR}
  ${TPL_INCLUDE_DIRS}
  )

# look for git.  This is used to configure version information into the
# executable and also to build upstream dependencies if necessary
set( SOPS_REPO_DATE "\"DATE NOT FOUND\"" )
set( SOPS_REPO_HASH "\"HASH NOT FOUND\"" )
find_package( Git )
if( GIT_FOUND )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%H\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE RESULT
    OUTPUT_VARIABLE SOPS_REPO_HASH 
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if( NOT ${RESULT} EQUAL 0 )
    set( SOPS_REPO_HASH "\"HASH NOT FOUND\"" )
  endif()
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%cd\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE SOPS_REPO_DATE 
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
   if( NOT ${RESULT} EQUAL 0 )
     set( SOPS_REPO_DATE "\"DATE NOT FOUND\"" )
   endif()
endif( GIT_FOUND )


# build the configure file using the information defined here.
configure_file(
  ${PROJECT_SOURCE_DIR}/SpatialOpsConfigure.h.in
  ${PROJECT_BINARY_DIR}/config/spatialops/SpatialOpsConfigure.h
  )
include_directories( ${PROJECT_BINARY_DIR}/config )

if( ENABLE_TESTS )
  enable_testing()
endif( ENABLE_TESTS )

add_subdirectory( spatialops )


# this will compile all of the test cases.  Use "make test" to execute
# them and see a summary.
if( ENABLE_TESTS )
  add_subdirectory( test )
endif()



#--- INSTALLATION SECTION ---
# "make install" performs these commands, along with any install commands inside child CMakeLists.txt files.
message( STATUS "Files will be installed to: " ${CMAKE_INSTALL_PREFIX} )

# install spatialops headers
install( FILES
  ${PROJECT_BINARY_DIR}/config/spatialops/SpatialOpsConfigure.h
  DESTINATION include/spatialops
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )


# ------------------------------------------------------------------
# configuration for projects using SpatialOps this section defines
# variables that are used for building the configuration file for the
# SpatialOps installation.  This can be used by other CMake build
# systems that use the SpatialOps library.

# get the location for the installed library
get_target_property( SOLocation spatialops LOCATION )

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
