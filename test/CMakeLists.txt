set( libs spatialops ${TPL_LIBRARIES} )

add_executable( simple_field SimpleFieldTest.cpp )
target_link_libraries( simple_field ${libs} )
add_test( simple_field simple_field )

add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )

add_executable( assign_test ScalingAssignmentTest.cpp )
target_link_libraries( assign_test ${libs} )

add_executable( stencil_test StencilTest.cpp )
target_link_libraries( stencil_test ${libs} )

add_executable( stencil_test_deep StencilTestDeepChain.cpp )
target_link_libraries( stencil_test_deep ${libs} )

add_executable( stencil_chain_test StencilChainTest.cpp )
target_link_libraries( stencil_chain_test ${libs} )

add_executable( till_10_test StencilChainTill10Test.cpp )
target_link_libraries( till_10_test ${libs} )

add_executable( till_20_test StencilChainTill20Test.cpp )
target_link_libraries( till_20_test ${libs} )

add_executable( till_5_test StencilChainTill5Test.cpp )
target_link_libraries( till_5_test ${libs} )

add_executable( cond_test cond_test.cpp )
target_link_libraries( cond_test ${libs} )

add_executable( error_test ErrorTest.cpp )
target_link_libraries( error_test ${libs} )

add_subdirectory( heatEqn )
