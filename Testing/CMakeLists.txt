
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_executable(CbctRecon_test
  cbctrecon_test.cpp
  cbctrecon_test.hpp
  )

add_cbctrecon_test(CbctRecon_test
  DATA{Input/Alderson/CT/,REGEX:.*}
  )

target_link_libraries(CbctRecon_test
  PRIVATE CbctReconLib
  )

target_include_directories(CbctRecon_test
  PUBLIC ${CBCTRECON_INCLUDE_DIRS}
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  )

install(TARGETS CbctRecon_test
  RUNTIME DESTINATION bin
  )
