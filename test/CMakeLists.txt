include_directories( ${spatialops_SOURCE_DIR}/include . )
#link_directories( ${spatialops_BINARY_DIR} )

# commonly required file:
set( common_SRC buildOps.cpp )

set( lagrange_SRC testLagrange.cpp )
add_executable( lagrange ${lagrange_SRC} )
target_link_libraries( lagrange ${SpatialOps_LIBRARY} )
add_test( lagrange_coefficients lagrange )

set( test_SRC testNew.cpp )
add_executable( testlong ${common_SRC} ${test_SRC} )
target_link_libraries( testlong ${SpatialOps_LIBRARY} )
#add_test( misc_tests testlong )

set( bctest_SRC bcTest.cpp )
add_executable( bctest ${bctest_SRC} ${common_SRC} )
target_link_libraries( bctest ${SpatialOps_LIBRARY} )
add_test( boundary_condition_tests bctest )

set( test1d_SRC testOneDim.cpp )
add_executable( test1d ${test1d_SRC} )
target_link_libraries( test1d ${SpatialOps_LIBRARY} )
add_test( 1D_operator_tests test1d )

if( LINALG_TRILINOS )
  set( poisson_SRC testPoisson.cpp )
  add_executable( poisson ${poisson_SRC} ${common_SRC} )
  target_link_libraries( poisson ${SpatialOps_LIBRARY} )
  add_test( poisson_system poisson )
endif( LINALG_TRILINOS )

set( filter_SRC tophat_test.cpp )
add_executable( testfilter ${filter_SRC} ${common_SRC} )
target_link_libraries( testfilter ${SpatialOps_LIBRARY} )
add_test( test_tophat_filter testfilter )

add_executable( restriction restriction_test.cpp )
target_link_libraries( restriction ${SpatialOps_LIBRARY} )
add_test( test_restriction restriction )