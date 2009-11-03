include_directories( ${spatialops_SOURCE_DIR}/include )

# commonly required file:
set( common_SRC buildOps.cpp )
add_library( testlib STATIC ${common_SRC} )

set( lagrange_SRC testLagrange.cpp )
add_executable( lagrange ${lagrange_SRC} )
target_link_libraries( lagrange ${SpatialOps_LIBRARY} )
add_test( lagrange_coefficients lagrange )

set( test_SRC testNew.cpp )
add_executable( testlong ${test_SRC} )
target_link_libraries( testlong ${SpatialOps_LIBRARY} testlib )
#add_test( misc_tests testlong )

set( bctest_SRC bcTest.cpp )
add_executable( bctest ${bctest_SRC} )
target_link_libraries( bctest ${SpatialOps_LIBRARY} testlib )
add_test( boundary_condition_tests bctest )

set( test1d_SRC testOneDim.cpp )
add_executable( test1d ${test1d_SRC} )
target_link_libraries( test1d ${SpatialOps_LIBRARY} )
add_test( 1D_operator_tests test1d )

if( LINALG_TRILINOS )
  set( poisson_SRC testPoisson.cpp )
  add_executable( poisson ${poisson_SRC}  )
  target_link_libraries( poisson ${SpatialOps_LIBRARY} testlib )
  add_test( poisson_system poisson )
endif( LINALG_TRILINOS )

set( filter_SRC tophat_test.cpp )
add_executable( testfilter ${filter_SRC}  )
target_link_libraries( testfilter ${SpatialOps_LIBRARY} testlib )
add_test( test_tophat_filter testfilter )

add_executable( restriction restriction_test.cpp )
target_link_libraries( restriction ${SpatialOps_LIBRARY} )
add_test( test_restriction restriction )

add_executable( fielddim field_dim_test.cpp )
target_link_libraries( fielddim ${SpatialOps_LIBRARY} )
add_test( field_dims fielddim )


add_executable( particleOps TestParticleOperators.cpp )
target_link_libraries( particleOps ${SpatialOps_LIBRARY} )
add_test( particle_ops particleOps )