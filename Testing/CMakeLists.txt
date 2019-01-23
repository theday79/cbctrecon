set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(test_target
  CbctRecon_test
  )

set(${test_target}_SRC
  cbctrecon_test.cpp
  cbctrecon_test.hpp
  )

set(${test_target}_DATA
  DATA{Input/Alderson/CT/,REGEX:.*}
  DATA{Input/Alderson/CBCT/,RECURSE:,REGEX:.*}
  )

add_cbctrecon_test(
    TARGET "${test_target}"
    SRC_FILES "${${test_target}_SRC}"
    DATA_ARGS "${${test_target}_DATA}"
  )

set_tests_properties(${test_target} PROPERTIES TIMEOUT 6000)

# Because ctest runs in CMAKE_BINARY_DIR we'll apparently also need the opencl file there
add_custom_command(TARGET ${test_target}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../Library/CbctReconLib/fdk_opencl.cl" "${CMAKE_BINARY_DIR}"
  )
install(TARGETS ${test_target}
  RUNTIME DESTINATION bin
  )

