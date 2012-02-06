set( libs spatialops ${TPL_LIBRARIES} )

add_executable( simple_field SimpleFieldTest.cpp )
target_link_libraries( simple_field ${libs} )
add_test( simple_field simple_field )

add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )

add_executable( assign_test ScalingAssignmentTest.cpp )
target_link_libraries( assign_test ${libs} )

add_executable( stencil_test StencilChainTest.cpp )
target_link_libraries( stencil_test ${libs} )

add_subdirectory( heatEqn )