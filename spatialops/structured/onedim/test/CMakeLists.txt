set( libs spatialops-onedim spatialops ${TPL_LIBRARIES} )

if( ENABLE_MATRIX )
  configure_file( results_test1d.gold ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  add_executable( test1d testOneDim.cpp )
  target_link_libraries( test1d ${libs} )
  add_test( 1D_operator_tests test1d )
endif( ENABLE_MATRIX )

add_executable( lagrange testLagrange.cpp )
target_link_libraries( lagrange ${libs} )
add_test( lagrange_coefficients lagrange )
