get_filename_component(_CBCTRECONExternalData_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${_CBCTRECONExternalData_DIR}/ExternalData.cmake)


set(ExternalData_BINARY_ROOT ${CMAKE_BINARY_DIR}/ExternalData)

set(ExternalData_URL_TEMPLATES "" CACHE STRING
  "Additional URL templates for the ExternalData CMake script to look for testing data. E.g.
file:///var/bigharddrive/%(algo)/%(hash)")
mark_as_advanced(ExternalData_URL_TEMPLATES)
list(APPEND ExternalData_URL_TEMPLATES
  # Data published by MIDAS
  
"https://midas3.kitware.com/midas/api/rest?method=midas.bitstream.download&checksum=%(hash)&algorithm=%(algo)"
  )

# Tell ExternalData commands to transform raw files to content links.
# TODO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT MD5)
