set( libs spatialops ${TPL_LIBRARIES} )

add_executable( fielddim field_dim_test.cpp )
target_link_libraries( fielddim ${libs} )
add_test( field_dims fielddim )

add_executable( simple_field SimpleFieldTest.cpp )
target_link_libraries( simple_field ${libs} )
add_test( simple_field simple_field )

add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )