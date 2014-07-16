nebo_cuda_prep_dir()

nebo_add_executable( example_1 Example1_field_creation.cpp )
target_link_libraries( example_1 spatialops )
add_test( example_1 example_1 )

nebo_add_executable( example_2 Example2_field_operations.cpp )
target_link_libraries( example_2 spatialops )
add_test( example_2 example_2 )

nebo_add_executable( example_3 Example3_field_reductions.cpp )
target_link_libraries( example_3 spatialops )
add_test( example_3 example_3 )

nebo_add_executable( example_35 Example3_5_field_types.cpp )
target_link_libraries( example_35 spatialops )
add_test( example_35 example_35 )

nebo_add_executable( example_field_type_inference Example4_field_type_inference.cpp )
target_link_libraries( example_field_type_inference spatialops )
add_test( example_field_type_inference example_field_type_inference )

nebo_add_executable( example_stencil_types Example5_stencil_type_inference.cpp )
target_link_libraries( example_stencil_types spatialops )
add_test( example_stencil_types example_stencil_types )

nebo_add_executable( example_stencils Example7_stencils.cpp )
target_link_libraries( example_stencils spatialops )
add_test( example_stencils example_stencils )

nebo_add_executable( example_masks Example8_masks.cpp )
target_link_libraries( example_masks spatialops )
add_test( example_masks example_masks )
