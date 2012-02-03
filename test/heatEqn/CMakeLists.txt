add_executable( heateqn main.cpp )
set( libs ${TPL_LIBRARIES} spatialops )

target_link_libraries( heateqn ${libs} )

add_test( heateqn heateqn )