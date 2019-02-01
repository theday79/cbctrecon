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
set(ZLIB_GIT_REPOSITORY "${git_protocol}hub.com/madler/zlib.git")
set(ZLIB_GIT_TAG master)

set(GDCM_GIT_REPOSITORY "${git_protocol}hub.com/malaterre/gdcm.git")
set(GDCM_GIT_TAG master)

set(DCMTK_GIT_REPOSITORY "${git_protocol}hub.com/DCMTK/DCMTK.git")
set(DCMTK_GIT_TAG master)

set(ITK_GIT_REPOSITORY "${git_protocol}hub.com/InsightSoftwareConsortium/ITK.git")
set(ITK_GIT_TAG release-4.13) #"d92873e33e8a54e933e445b92151191f02feab42") #"v4.13.0"

# set(DLIB_GIT_REPOSITORY "${git_protocol}hub.com/davisking/dlib.git")
# set(DLIB_GIT_TAG master)

set(plm Plastimatch)
set(${plm}_GIT_REPOSITORY "${git_protocol}lab.com/plastimatch/plastimatch.git")
set(${plm}_GIT_TAG master) #"b9b6c7f59b98d52a93d58288af9d5bc25d00f24b")

set(TINYREFL_GIT_REPOSITORY "${git_protocol}lab.com/Manu343726/tinyrefl.git")
set(TINYREFL_GIT_TAG master)

# ZLIB
external_proj(ZLIB)

# GDCM
option(USE_SYSTEM_GDCM "ON=Use system GDCM, OFF=Download GDCM" OFF)
if (NOT USE_SYSTEM_ITK) # Because ITK also defines the gdcm symbols
  external_proj(GDCM)
endif()

# DCMTK
external_proj(DCMTK)

# ITK/RTK
external_proj(ITK)
if (NOT USE_SYSTEM_ITK)
  configure_file(
    ${CbctRecon_SOURCE_DIR}/External/patches/ExternalRTK.remote.cmake
    ${ITK_SOURCE_DIR}/Modules/Remote/
    COPYONLY
  )
else()
  include(${ITK_USE_FILE})
endif()

# dlib because the one in Plastimatch is two years old
# external_proj(DLIB)

# Plastimatch
external_proj(${plm})

set(${plm}_DIR ${${plm}_BUILD_DIR})
set(${plm}_SRC ${${plm}_SOURCE_DIR})

if(USE_TINYREFL)
  external_proj(TINYREFL)
endif()
