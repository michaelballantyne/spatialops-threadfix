set( libs spatialops ${TPL_LIBRARIES} )

add_executable( simple_field SimpleFieldTest.cpp )
target_link_libraries( simple_field ${libs} )
add_test( simple_field simple_field )

add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )

add_executable( assign_test ScalingAssignmentTest.cpp )
target_link_libraries( assign_test ${libs} )

add_executable( stencil_test1 StencilTest1.cpp )
target_link_libraries( stencil_test1 ${libs} )

add_executable( stencil_test2 StencilTest2.cpp )
target_link_libraries( stencil_test2 ${libs} )

add_executable( stencil_test3 StencilTest3.cpp )
target_link_libraries( stencil_test3 ${libs} )

add_executable( stencil_test4 StencilTest4.cpp )
target_link_libraries( stencil_test4 ${libs} )

add_executable( stencil_test5 StencilTest5.cpp )
target_link_libraries( stencil_test5 ${libs} )

add_executable( stencil_test6 StencilTest6.cpp )
target_link_libraries( stencil_test6 ${libs} )

add_executable( stencil_test7 StencilTest7.cpp )
target_link_libraries( stencil_test7 ${libs} )

add_executable( stencil_test8 StencilTest8.cpp )
target_link_libraries( stencil_test8 ${libs} )

add_executable( stencil_test9 StencilTest9.cpp )
target_link_libraries( stencil_test9 ${libs} )

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
