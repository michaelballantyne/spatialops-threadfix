set( libs spatialops-structured spatialops ${TPL_LIBRARIES} )

include_directories( ${spatialops_SOURCE_DIR}/include . )

add_executable( test_field testField.cpp )
target_link_libraries( test_field ${libs} )
add_test( test_field test_field )

add_executable( test_memwindow testMemWindow.cpp )
target_link_libraries( test_memwindow ${libs} )
add_test( test_memwindow test_memwindow )

add_executable( test_fvtools testFVTools.cpp )
target_link_libraries( test_fvtools ${libs} )
#jcs add this once we get the new SpatialOps stuff done.
# add_test( test_fvtools test_fvtools )

if( ENABLE_BOOST_SERIALIZATION )
  add_executable( test_serialize testSerialization.cpp )
  target_link_libraries( test_serialize ${libs} )
  add_test( test_serialization test_serialize )
endif()