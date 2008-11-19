cmake_minimum_required(VERSION 2.6)
project(spatialops)

#=====================================================================
#---------------------------------------------------------------------
# BEGIN USER MODIFIABLE BLOCK

#-- Set the path for the trilinos installation on your machine
set( Trilinos
  /jcs/software/trilinos
  )

#-- Set the path for the Daixtrose installation on your machine
set( Daixtrose
  /jcs/software/daixtrose-0.0.3/jcs_install
  )

# hack for mac. If you are building for a mac, you will need to
# specify the path to your f95 library here.
if(APPLE)
  link_directories(
    /sw/lib/gcc-lib/i386-apple-darwin8/4.0.3
    )
  enable_language( Fortran C CXX )
  message( STATUS "adding -lf95 to linker line" )
  set( EXTRA_LIBS f95 )
endif(APPLE)

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
add_library( spatialops STATIC ${SO_Sources} )
target_link_libraries( spatialops aztecoo epetraext epetra blas lapack teuchos ${EXTRA_LIBS} )

# this will compile all of the test cases.  Use "make test" to execute
# them and see a summary.
enable_testing()
add_subdirectory( test )


#--- INSTALLATION SECTION ---
# "make install" performs these commands, along with any install commands inside child CMakeLists.txt files.
# Installation - headers
file( GLOB so_headers "include/spatialops/*.h" )
install( FILES ${so_headers}
  DESTINATION include/spatialops )

# INSTALLATION - libraries
install( TARGETS spatialops
  ARCHIVE DESTINATION lib )


message( STATUS "Files will be installed to: " ${CMAKE_INSTALL_PREFIX} )
