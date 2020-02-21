include(${CbctRecon_SOURCE_DIR}/cmake/externals.cmake)

set(git_protocol https://git)

macro(external_proj NAME)
  if(USE_SYSTEM_${NAME})
    find_package(${NAME} REQUIRED)
  else()
    include(External_${NAME}.cmake)
    external_dependency(${NAME} ${${NAME}_GIT_REPOSITORY} ${${NAME}_GIT_TAG})
  endif()
endmacro()

# REPOSITORIES
# set(ZLIB_GIT_REPOSITORY "${git_protocol}hub.com/madler/zlib.git")
# set(ZLIB_GIT_TAG master)

set(GDCM_GIT_REPOSITORY "${git_protocol}hub.com/malaterre/gdcm.git")
set(GDCM_GIT_TAG master)

set(DCMTK_GIT_REPOSITORY "${git_protocol}hub.com/DCMTK/DCMTK.git")
set(DCMTK_GIT_TAG master)

# RTK is deprecating ITKv4 support (may still work), but 5.0.0 and forward should be fine
set(ITK_GIT_REPOSITORY "${git_protocol}hub.com/InsightSoftwareConsortium/ITK.git")
set(ITK_GIT_TAG master)

# It looks like version 19.3 and forward would work with C++17
set(dlib_GIT_REPOSITORY "${git_protocol}hub.com/davisking/dlib.git")
set(dlib_GIT_TAG master)

set(plm Plastimatch)
set(${plm}_GIT_REPOSITORY "${git_protocol}lab.com/plastimatch/plastimatch.git")
set(${plm}_GIT_TAG master)

set(TINYREFL_GIT_REPOSITORY "${git_protocol}lab.com/Manu343726/tinyrefl.git")
set(TINYREFL_GIT_TAG master)

set(ABSEIL_GIT_REPOSITORY "${git_protocol}hub.com/abseil/abseil-cpp.git")
set(ABSEIL_GIT_TAG master)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  external_proj(ABSEIL)
endif()

if(WIN32 AND NOT (USE_ITK_GDCM OR USE_SYSTEM_ITK))
  external_proj(GDCM)
  if (NOT USE_SYSTEM_GDCM)
    set(GDCM_DIR ${GDCM_BUILD_DIR} CACHE PATH "")
  endif()
endif()

# DCMTK
if(NOT USE_ITK_DCMTK)
  external_proj(DCMTK)
  if(NOT USE_SYSTEM_DCMTK)
    set(DCMTK_DIR "${CMAKE_BINARY_DIR}")
    set(DCMTK_TARGETS_IMPORTED 1)
    include_directories("${DCMTK_BINARY_DIR}/config/include/")
  endif()
endif()

# ITK/RTK
external_proj(ITK)
if (NOT USE_SYSTEM_ITK)
    if (${ITK_GIT_TAG} STREQUAL "release-4.13")
    configure_file(
      ${CbctRecon_SOURCE_DIR}/External/patches/ExternalRTK.remote.cmake
      ${ITK_SOURCE_DIR}/Modules/Remote/
      COPYONLY
    )
  endif()
else()
  include(${ITK_USE_FILE})
endif()

if(USE_ITK_DCMTK)
  set(DCMTK_DIR "${ITK_BINARY_DIR}/Modules/ThirdParty/DCMTK/ITKDCMTK_ExtProject-build")
  set(DCMTK_TARGETS_IMPORTED 1)
endif()

# dlib because the one in Plastimatch is two years old
external_proj(dlib)

# Plastimatch
external_proj(${plm})

set(${plm}_DIR ${${plm}_BUILD_DIR})
set(${plm}_SRC ${${plm}_SOURCE_DIR})

if(USE_TINYREFL)
  external_proj(TINYREFL)
endif()

