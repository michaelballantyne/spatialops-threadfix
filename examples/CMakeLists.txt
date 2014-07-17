nebo_cuda_prep_dir()

nebo_add_executable( example1_field_creation Example1_field_creation.cpp )
target_link_libraries( example1_field_creation spatialops )
add_test( example1_field_creation example1_field_creation )
install( TARGETS example1_field_creation RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example2_field_operations Example2_field_operations.cpp )
target_link_libraries( example2_field_operations spatialops )
add_test( example2_field_operations example2_field_operations )
install( TARGETS example2_field_operations RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example3_field_reductions Example3_field_reductions.cpp )
target_link_libraries( example3_field_reductions spatialops )
add_test( example3_field_reductions example3_field_reductions )
install( TARGETS example2_field_operations RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example3_5_field_types Example3_5_field_types.cpp )
target_link_libraries( example3_5_field_types spatialops )
add_test( example3_5_field_types example3_5_field_types )
install( TARGETS example3_5_field_types RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example4_field_type_inference Example4_field_type_inference.cpp )
target_link_libraries( example4_field_type_inference spatialops )
add_test( example4_field_type_inference example4_field_type_inference )
install( TARGETS example4_field_type_inference RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example5_stencil_types Example5_stencil_type_inference.cpp )
target_link_libraries( example5_stencil_types spatialops )
add_test( example5_stencil_types example5_stencil_types )
install( TARGETS example5_stencil_types RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example7_stencils Example7_stencils.cpp )
target_link_libraries( example7_stencils spatialops )
add_test( example7_stencils example_stencils )
install( TARGETS example7_stencils RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example8_masks Example8_masks.cpp )
target_link_libraries( example8_masks spatialops )
#add_test( example8_masks example8_masks )
#install( TARGETS example8_stencils RUNTIME DESTINATION bin/spatialops )

nebo_add_executable( example9_heat_equation Example9_heat_equation.cpp )
target_link_libraries( example9_heat_equation spatialops )
add_test( example9_heat_equation example_heat_equation )
install( TARGETS example9_heat_equation RUNTIME DESTINATION bin/spatialops )
