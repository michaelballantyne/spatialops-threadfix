include_directories( ${spatialops_SOURCE_DIR}/include . )

add_executable( test_field testField.cpp )
target_link_libraries( test_field ${SpatialOps_LIBRARY} ${TPL_LIBRARIES} )
add_test( test_field test_field )

add_executable( test_memwindow testMemWindow.cpp )
target_link_libraries( test_memwindow ${SpatialOps_LIBRARY} )
add_test( test_memwindow test_memwindow )

add_executable( test_memw_ops testOperator.cpp )
target_link_libraries( test_memw_ops testlib ${SpatialOps_LIBRARY} )
add_test( test_memw_ops test_memw_ops )

add_executable( test_stencil2 test_stencil2.cpp Grid.cpp )
target_link_libraries( test_stencil2 ${SpatialOps_LIBRARY} )
add_test( stencil2_bcxyz   test_stencil2  --bcx --bcy --bcz )
add_test( stencil2_bcxy    test_stencil2  --bcx --bcy )
add_test( stencil2_bcxz    test_stencil2  --bcx --bcz )
add_test( stencil2_bcyz    test_stencil2  --bcz --bcy )
add_test( stencil2_bcx     test_stencil2  --bcx )
add_test( stencil2_bcy     test_stencil2  --bcy )
add_test( stencil2_bcz     test_stencil2  --bcz )
add_test( stencil2_no_bc   test_stencil2  --bcz )
add_test( stencil2_no_bc2  test_stencil2  --nx 5 --ny 6 --nz 7 )
add_test( stencil2_no_bc3  test_stencil2  --nx 6 --ny 7 --nz 5 )

add_executable( test_stencil4 test_stencil4.cpp Grid.cpp )
target_link_libraries( test_stencil4 ${SpatialOps_LIBRARY} )


if( ENABLE_BOOST_SERIALIZATION )
  add_executable( test_serialize testSerialization.cpp )
  target_link_libraries( test_serialize ${SpatialOps_LIBRARY} ${Boost_LIBRARIES} )
  add_test( test_serialization test_serialize )
endif()