set( libs spatialops-onedim spatialops ${TPL_LIBRARIES} )
if( ENABLE_CUDA )
	set(libs ${libs} CUDAMemoryManager)
endif( ENABLE_CUDA )

add_executable( lagrange testLagrange.cpp )
target_link_libraries( lagrange ${libs} )
add_test( lagrange_coefficients lagrange )
