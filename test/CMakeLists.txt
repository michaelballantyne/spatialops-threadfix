nebo_cuda_prep_dir()

set( libs spatialops spatialops-structured spatialops-stencil ${TPL_LIBRARIES} )

nebo_add_executable( gpu_test GPUtest.cpp )
target_link_libraries( gpu_test ${libs} )

nebo_add_executable( simple_field SimpleFieldTest.cpp )
target_link_libraries( simple_field ${libs} )
add_test( simple_field simple_field )

nebo_add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )

nebo_add_executable( assign_test ScalingAssignmentTest.cpp )
target_link_libraries( assign_test ${libs} )

nebo_add_executable( stencil_test StencilTest.cpp )
target_link_libraries( stencil_test ${libs} )

nebo_add_executable( stencil_test_deep StencilTestDeepChain.cpp )
target_link_libraries( stencil_test_deep ${libs} )

nebo_add_executable( stencil_chain_test StencilChainTest.cpp )
target_link_libraries( stencil_chain_test ${libs} )

nebo_add_executable( cond_test cond_test.cpp )
target_link_libraries( cond_test ${libs} )

add_executable( filter_test FilterTest.cpp )
target_link_libraries( filter_test ${libs} )

add_subdirectory( heatEqn )
