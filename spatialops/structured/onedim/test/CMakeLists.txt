set( libs spatialops-onedim spatialops ${TPL_LIBRARIES} )

add_executable( lagrange testLagrange.cpp )
target_link_libraries( lagrange ${libs} )
add_test( lagrange_coefficients lagrange )
