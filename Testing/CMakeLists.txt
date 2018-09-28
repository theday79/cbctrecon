set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_AUTOMOC ON)

set(test_target
  CbctRecon_test
  )

set(${test_target}_SRC
  cbctrecon_test.cpp
  cbctrecon_test.hpp
  )

set(${test_target}_DATA
  DATA{Input/Alderson/CT/,REGEX:.*}
  )
  
add_cbctrecon_test(
    TARGET "${test_target}"
    SRC_FILES "${${test_target}_SRC}"
    DATA_ARGS "${${test_target}_DATA}"
  )

install(TARGETS ${test_target}
  RUNTIME DESTINATION bin
  )

