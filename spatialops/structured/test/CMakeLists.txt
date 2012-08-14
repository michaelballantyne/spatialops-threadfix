set( libs spatialops-structured spatialops ${TPL_LIBRARIES} )

include_directories( ${spatialops_SOURCE_DIR}/include . )

add_executable( test_field testField.cpp )
target_link_libraries( test_field ${libs} )
add_test( test_field test_field )

add_executable( test_memwindow testMemWindow.cpp )
target_link_libraries( test_memwindow ${libs} )
add_test( test_memwindow test_memwindow )

add_executable( test_grid testGrid.cpp )
target_link_libraries( test_grid ${libs} )
add_test( test_grid test_grid )

add_executable( test_tools testTools.cpp )
target_link_libraries( test_tools ${TPL_LIBRARIES} )
add_test( test_tools test_tools )

if( ENABLE_BOOST_SERIALIZATION )
  add_executable( test_serialize testSerialization.cpp )
  target_link_libraries( test_serialize ${libs} )
  add_test( test_serialization test_serialize )
endif()

macro( build_test tname exec tout )
  add_test( CUDA_${tname} ${exec} ) 

  set_tests_properties ( CUDA_${tname} PROPERTIES PASS_REGULAR_EXPRESSION ${tout} )
endmacro()

if( CUDA_CUDA_LIBRARY )
  include_directories( ${CUDA_INCLUDE_DIRS} )
  cuda_add_library( cudaVecLib STATIC CudaVecLib.cu )
  target_link_libraries( cudaVecLib ${CUDA_LIBRARIES} )

  add_executable( ctest ctest.cpp )
  set( driver_link_libs cudaVecLib ${CUDA_LIBRARIES} )
  target_link_libraries( ctest ${driver_link_libs} )
  
  set( cxxtemp ${CMAKE_CXX_FLAGS} )
  set( CMAKE_CXX_FLAGS "-DDEBUG_CUDA_VERBOSE ${cxxtemp}" )
  add_executable( debug_memtest_flags memtest1.cpp )
  target_link_libraries( debug_memtest_flags spatialops )
  set( CMAKE_CXX_FLAGS ${cxxtemp} )
  
  add_executable( simplememtest memtest1.cpp )
  target_link_libraries( simplememtest spatialops )
  
  add_executable( sharedmemtest memtest2.cpp )
  target_link_libraries( sharedmemtest spatialops )
  
  #add_executable( extfieldtest testExternalField.cpp )
  #target_link_libraries( extfieldtest ${libs} )
  
  cuda_add_executable( extstenciltest ExternalStencilTest.cu )
  set( driver_link_libs ${libs} ${CUDA_LIBRARIES} )
  target_link_libraries( extstenciltest ${driver_link_libs} )

  build_test( CompileTest ctest "Success" )
  build_test( MemoryDebugFlags debug_memtest_flags "Success" )
  build_test( MemoryTest simplememtest "Success" )
  build_test( SharedMemoryTest sharedmemtest "Success" )
  #build_test( ExternalFieldTest extfieldtest "Success" )
  build_test( ExternalStencilTest extstenciltest "Success" )
  
else( CUDA_CUDA_LIBRARY )
  message( STATUS "Must include CUDA libraries to build cuda test base" )
endif( CUDA_CUDA_LIBRARY )
