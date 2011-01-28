set( libs spatialops-structured spatialops ${TPL_LIBRARIES} )


# commonly required file:
set( common_SRC buildOps.cpp Grid.cpp )
add_library( testlib STATIC ${common_SRC} )


add_executable( test_memw_ops testOperator.cpp )
target_link_libraries( test_memw_ops testlib ${libs} )
add_test( test_memw_ops test_memw_ops )

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
