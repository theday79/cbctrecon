if (USE_TINYREFL)
  include(${TINYREFL_SOURCE_DIR}/tool/driver.cmake)
endif()

add_cbctrecon_test(CbctRecon_test
  cbctrecon_test.cpp
  )

set(APP_DIR ${CbctRecon_SOURCE_DIR}/Applications/CbctRecon)
set(APP_FOLDERS
  ${APP_DIR}
  ${APP_DIR}/qt-utils/
  ${APP_DIR}/rtkExtension/
  ${APP_DIR}/ui/
)

target_link_libraries(CbctRecon_test PRIVATE CbctReconLib)
target_include_directories(CbctRecon_test
  PRIVATE ${APP_FOLDERS}
  )

if(USE_TINYREFL)
  tinyrefl_tool(TARGET CbctRecon_test
    HEADERS
    cbctrecon_test.h
    ${APP_DIR}/cbctrecon.h
    )
endif()
