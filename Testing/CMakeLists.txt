
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

target_link_libraries(CbctRecon_test
  PRIVATE CbctReconLib
  PRIVATE ${PLM_TARGETS}
  PRIVATE ${CBCTRECON_DEP_LIBS}
  )

target_include_directories(CbctRecon_test
  PRIVATE ${APP_FOLDERS}
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  )
