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

add_executable( test_stencil test_stencil.cpp Grid.cpp )
target_link_libraries( test_stencil ${SpatialOps_LIBRARY} )

if( ENABLE_BOOST_SERIALIZATION )
  add_executable( test_serialize testSerialization.cpp )
  target_link_libraries( test_serialize ${SpatialOps_LIBRARY} ${Boost_LIBRARIES} )
  add_test( test_serialization test_serialize )
endif()