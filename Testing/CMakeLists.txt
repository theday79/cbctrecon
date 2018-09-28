set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_cbctrecon_test(
    TNAME CbctRecon_test
    SRC_FILES cbctrecon_test.cpp cbctrecon_test.hpp
    DATA_ARGS DATA{Input/Alderson/CT/,REGEX:.*}
  )

install(TARGETS CbctRecon_test
  RUNTIME DESTINATION bin
  )

