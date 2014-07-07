nebo_cuda_prep_dir()

nebo_add_executable( example_1 Example1_field_creation.cpp )
target_link_libraries( example_1 spatialops )
add_test( example_1 example_1 )

nebo_add_executable( example_2 Example1_field_creation.cpp )
target_link_libraries( example_2 spatialops )
add_test( example_2 example_2 )
