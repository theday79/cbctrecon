
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_cbctrecon_test(CbctRecon_test
  cbctrecon_test.cpp
  cbctrecon_test.hpp
  )

set(APP_DIR ${CbctRecon_SOURCE_DIR}/Applications/CbctRecon)
set(APP_FOLDERS
  ${APP_DIR}
  ${APP_DIR}/qt-utils/
  ${APP_DIR}/rtkExtension/
  ${APP_DIR}/ui/
)
set(GEN_UI_INCLUDE_DIR ${CMAKE_BINARY_DIR}/Applications/CbctRecon/CbctReconLib_autogen/include)

target_link_libraries(CbctRecon_test
  PRIVATE CbctReconLib
  )
target_include_directories(CbctRecon_test
  PRIVATE ${APP_FOLDERS}
  PRIVATE ${GEN_UI_INCLUDE_DIR}
  )
