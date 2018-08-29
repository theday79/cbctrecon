if (USE_TINYREFL)
  include(${TINYREFL_SOURCE_DIR}/tool/driver.cmake)
endif()

add_cbctrecon_test(CbctRecon_test
  cbctrecon_test.cpp
  )

set(APP_DIR ${CbctRecon_SOURCE_DIR}/Applications/CbctRecon)

target_link_libraries(CbctRecon_test PRIVATE CbctReconLib)
target_include_directories(CbctRecon_test PUBLIC ${APP_DIR})

if(USE_TINYREFL)
  tinyrefl_tool(TARGET CbctRecon_test
    HEADERS
    cbctrecon_test.h
    ${APP_DIR}/cbctrecon.h
    )
endif()
