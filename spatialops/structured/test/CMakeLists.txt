include_directories( ${spatialops_SOURCE_DIR}/include )

add_executable( test_field testField.cpp )
target_link_libraries( test_field ${SpatialOps_LIBRARY} ${TPL_LIBRARIES} )

add_test( test_field test_field )