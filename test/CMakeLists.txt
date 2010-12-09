include_directories( ${spatialops_SOURCE_DIR}/include )

set( libs ${SpatialOps_LIBRARY} ${TPL_LIBRARIES} )

# commonly required file:
set( common_SRC buildOps.cpp Grid.cpp )
add_library( testlib STATIC ${common_SRC} )

set( lagrange_SRC testLagrange.cpp )
add_executable( lagrange ${lagrange_SRC} )
target_link_libraries( lagrange ${libs} )
add_test( lagrange_coefficients lagrange )

add_executable( bcOpTest testOpBC.cpp )
target_link_libraries( bcOpTest testlib ${libs} )
add_test( test_bc_op bcOpTest )

set( test_SRC testNew.cpp )
add_executable( longtest ${test_SRC} )
target_link_libraries( longtest testlib ${libs} )
#add_test( misc_tests longtest )

set( bctest_SRC bcTest.cpp )
add_executable( bctest ${bctest_SRC} )
target_link_libraries( bctest testlib ${libs} )
add_test( boundary_condition_tests bctest )

configure_file( results_test1d.gold ${CMAKE_BINARY_DIR}/test COPYONLY )
set( test1d_SRC testOneDim.cpp )
add_executable( test1d ${test1d_SRC} )
target_link_libraries( test1d ${libs} )
add_test( 1D_operator_tests test1d )

if( LINALG_TRILINOS )
  set( poisson_SRC testPoisson.cpp )
  add_executable( poisson ${poisson_SRC}  )
  target_link_libraries( poisson testlib ${libs} )
  add_test( poisson_system poisson )
endif( LINALG_TRILINOS )

set( filter_SRC tophat_test.cpp )
add_executable( testfilter ${filter_SRC}  )
target_link_libraries( testfilter testlib ${libs} )
add_test( test_tophat_filter testfilter )

add_executable( restriction restriction_test.cpp )
target_link_libraries( restriction ${libs} )
add_test( test_restriction restriction )

add_executable( fielddim field_dim_test.cpp )
target_link_libraries( fielddim ${libs} )
add_test( field_dims fielddim )


add_executable( particleOps TestParticleOperators.cpp )
target_link_libraries( particleOps ${libs} )
add_test( particle_ops particleOps )


add_executable( simple_field SimpleFieldTest.cpp )
target_link_libraries( simple_field ${libs} )
add_test( simple_field simple_field )

add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )