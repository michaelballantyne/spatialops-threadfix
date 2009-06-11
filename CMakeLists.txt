cmake_minimum_required(VERSION 2.6)
project(spatialops)

set( MAJOR_VERSION "0" )
set( MINOR_VERSION "1" )
set( PATCH_VERSION "0" )

#=====================================================================
#---------------------------------------------------------------------
# BEGIN USER MODIFIABLE BLOCK

set( CMAKE_BUILD_TYPE  # choose one type below.
  Release
# Debug
# RelWithDebInfo
  )

#-- Set the path for the trilinos installation on your machine
set( Trilinos
  /jcs/software/trilinos902_opt
  )

#-- Set the path for the Daixtrose installation on your machine
set( Daixtrose
  /jcs/software/daixtrose-0.0.3/jcs_install
  )

#-- Set the location for the installation.  Comment this out to install to /usr/local
set( CMAKE_INSTALL_PREFIX
  ${CMAKE_CURRENT_BINARY_DIR}
  )

# END USER MODIFIABLE BLOCK
#---------------------------------------------------------------------
#=====================================================================




#--------------------------------------
# resolve required libraries

 # setup the boost library
set(Boost_USE_MULTITHREAD ON)
set(Boost_ADDITIONAL_VERSIONS "1.37.0" "1.37" )
find_package( Boost REQUIRED 
  COMPONENTS  # don't require any compiled boost libs - only headers.
  ) 

# until we get a package finder for trilinos, we will need these:
set( Trilinos_INCLUDE_DIR ${Trilinos}/include )
set( Trilinos_LIBRARY_DIR ${Trilinos}/lib )

# until we get a package finder for daixtrose, we will need this:
set( Daixtrose_INCLUDE_DIR ${Daixtrose}/include )
#--------------------------------------

include_directories( ${Trilinos_INCLUDE_DIR} ${Daixtrose_INCLUDE_DIR} ${Boost_INCLUDE_DIR} include )
link_directories( ${Trilinos_LIBRARY_DIR} )


# setup the build for the spatialops library.  Note that the linker
# order is important, and future revisions to libraries such as
# trilinos may require a modification of link order.
file( GLOB SO_Sources "src/*.cpp" )
add_definitions( -DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR )
if( CMAKE_BUILD_TYPE STREQUAL Release )
  add_definitions( -O3 -fexpensive-optimizations )
  message( STATUS "turning on aggressive optimizations" )
endif( CMAKE_BUILD_TYPE STREQUAL Release )
add_library( spatialops STATIC ${SO_Sources} )
target_link_libraries( spatialops aztecoo epetraext epetra blas lapack teuchos ${EXTRA_LIBS} )

# this will compile all of the test cases.  Use "make test" to execute
# them and see a summary.
enable_testing()
add_subdirectory( test )


#--- INSTALLATION SECTION ---
# "make install" performs these commands, along with any install commands inside child CMakeLists.txt files.
message( STATUS "Files will be installed to: " ${CMAKE_INSTALL_PREFIX} )

# package spatialops headers
file( GLOB so_headers "include/spatialops/*.h" )
install( FILES ${so_headers}
  DESTINATION include/spatialops
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

# package the spatialops library
install( TARGETS spatialops ARCHIVE DESTINATION lib )

# package daixtrose locally so that we can keep it locally
install( DIRECTORY ${Daixtrose_INCLUDE_DIR}/daixtrose
  DESTINATION include
  FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

# package trilinos components required by spatialops
file( GLOB trilinos_headers
  ${Trilinos_INCLUDE_DIR}/Epetra*.h
  ${Trilinos_INCLUDE_DIR}/Teuchos*.h*
  ${Trilinos_INCLUDE_DIR}/Az*.h
  )
set( trilinos_libs
  ${Trilinos_LIBRARY_DIR}/libaztecoo.a
  ${Trilinos_LIBRARY_DIR}/libepetra.a
  ${Trilinos_LIBRARY_DIR}/libepetraext.a
  ${Trilinos_LIBRARY_DIR}/libteuchos.a
  )
install( FILES ${trilinos_headers}
  DESTINATION include
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )
install( FILES ${trilinos_libs} DESTINATION lib )



# Information for creating packages
set( CPACK_PACKAGE_VENDOR "James C. Sutherland (james.sutherland@utah.edu" )
set( CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION} )
set( CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION} )
set( CPACK_PACKAGE_VERSION_PATCH ${PATCH_VERSION} )
include( CPack )
