if(HAVE_AVX2)
  set(USE_AVX_INSTRUCTIONS
      ON
      CACHE BOOL "" FORCE)
endif()

# CUDA
if(USE_CUDA)
  set(DLIB_USE_CUDA
      ON
      CACHE BOOL "" FORCE)
else()
  set(DLIB_USE_CUDA
      OFF
      CACHE BOOL "" FORCE)
endif()
