add_executable( heateqn main.cpp )
target_link_libraries( heateqn
  ${TPL_LIBRARIES}
  spatialops
  )

add_test( heateqn heateqn )