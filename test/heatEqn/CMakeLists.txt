add_executable( heateqn main.cpp )
set( libs ${TPL_LIBRARIES} spatialops )
if( ENABLE_CUDA )
	set(libs ${libs} spatialops_cuda )
endif( ENABLE_CUDA )

target_link_libraries( heateqn ${libs} )

add_test( heateqn heateqn )