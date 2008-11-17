include_directories( ${spatialops_SOURCE_DIR}/include . )
link_directories( ${spatialops_BINARY_DIR} )

# commonly required file:
set( common_SRC buildOps.cpp )

set( lagrange_SRC testLagrange.cpp )
add_executable( lagrange ${lagrange_SRC} )
target_link_libraries( lagrange spatialops )
add_test( lagrange_coefficients lagrange )

set( test_SRC testNew.cpp )
add_executable( test ${common_SRC} ${test_SRC} )
target_link_libraries( test spatialops )
#add_test( operator_suite_test test )

set( bctest_SRC bcTest.cpp )
add_executable( bctest ${bctest_SRC} ${common_SRC} )
target_link_libraries( bctest spatialops )
add_test( boundary_condition_tests bctest )

set( test1d_SRC testOneDim.cpp )
add_executable( test1d ${test1d_SRC} )
target_link_libraries( test1d spatialops )
add_test( 1D_operator_tests test1d )

set( poisson_SRC testPoisson.cpp )
add_executable( poisson ${poisson_SRC} ${common_SRC} )
target_link_libraries( poisson spatialops )
add_test( poisson_system poisson )
