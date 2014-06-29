nebo_cuda_prep_dir()

set( libs spatialops spatialops-structured spatialops-stencil ${TPL_LIBRARIES} )

nebo_add_executable( field_exprs testFieldExprs.cpp )
target_link_libraries( field_exprs ${libs} )
add_test( field_exprs field_exprs )

nebo_add_executable( test_nebo NeboTest.cpp )
target_link_libraries( test_nebo ${libs} )
add_test( test_nebo test_nebo )

nebo_add_executable( test_pointwise_simple SimplePointwiseTest.cpp )
target_link_libraries( test_pointwise_simple ${libs} )
add_test( test_pointwise_simple test_pointwise_simple )

add_subdirectory( heatEqn )
