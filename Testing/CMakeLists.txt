if (USE_TINYREFL)
  include(${TINYREFL_SOURCE_DIR}/tool/driver.cmake)
endif()

add_cbctrecon_test(CbctRecon_test
  cbctrecon_test.cpp
  )

target_link_libraries(CbctRecon_test PRIVATE CbctRecon)

tinyrefl_tool(TARGET CbctRecon_test
  HEADERS cbctrecon_test.h ${CbctRecon_HEADERS}
  )
