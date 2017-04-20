SET(DEPENDENCIES DCMTK ITKv4 RTK Plastimatch gPMC)
mark_as_superbuild(DEPENDENCIES:STRING)
include(ExternalProjectAddSource)

include(ListToString)

ExternalProject_Include_Dependencies(CbctRecon DEPENDS_VAR DEPENDENCIES)


#------------------------------------------------------------------------------
# Configure and build CbctRecon
#------------------------------------------------------------------------------
set(proj CbctRecon)

ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  DEPENDS ${DEPENDENCIES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${CbctRecon_BINARY_INNER_SUBDIR}
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
    -DADDITIONAL_C_FLAGS:STRING=${ADDITIONAL_C_FLAGS}
    -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
    -DCbctRecon_REQUIRED_C_FLAGS:STRING=${CbctRecon_REQUIRED_C_FLAGS}
    -DCbctRecon_REQUIRED_CXX_FLAGS:STRING=${CbctRecon_REQUIRED_CXX_FLAGS}
    -DCbctRecon_SUPERBUILD:BOOL=OFF
    -DCbctRecon_SUPERBUILD_DIR:PATH=${CbctRecon_BINARY_DIR}
    -D${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME:STRING=${${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_MAJOR:STRING=${${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_MAJOR}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_MINOR:STRING=${${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_MINOR}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_PATCH:STRING=${${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_PATCH}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_TWEAK:STRING=${${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_TWEAK}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_RC:STRING=${${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_RC}
    -DQt5_DIR:STRING=${Qt5_DIR}
    -DITK_DIR:STRING=${ITK_DIR}
    -DRTK_DIR:STRING=${RTK_DIR}
    -DDCMTK_DIR:STRING=${DCMTK_DIR}
    -DPlastimatch_DIR:STRING=${Plastimatch_DIR}
    -DPlastimatch_BUILD_DIR:STRING=${Plastimatch_BUILD_DIR}
    -DPlastimatch_SRC:STRING=${Plastimatch_SRC}
    -DCbctRecon_APPLICATIONS_DIR:PATH=${CbctRecon_APPLICATIONS_DIR}
    -DCbctRecon_EXTENSION_SOURCE_DIRS:STRING=${CbctRecon_EXTENSION_SOURCE_DIRS}
    -DExternalData_OBJECT_STORES:PATH=${ExternalData_OBJECT_STORES}
    ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
  INSTALL_COMMAND ""
  )

# This custom external project step forces the build and later
# steps to run whenever a top level build is done...
#
# BUILD_ALWAYS flag is available in CMake 3.1 that allows force build
# of external projects without this workaround. Remove this workaround
# and use the CMake flag instead, when CbctRecon's required minimum CMake
# version will be at least 3.1.
#
if(CMAKE_CONFIGURATION_TYPES)
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${CMAKE_CFG_INTDIR}/${proj}-build")
else()
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build")
endif()
ExternalProject_Add_Step(${proj} forcebuild
  COMMAND ${CMAKE_COMMAND} -E remove ${BUILD_STAMP_FILE}
  COMMENT "Forcing build step for '${proj}'"
  DEPENDEES build
  ALWAYS 1
  )
