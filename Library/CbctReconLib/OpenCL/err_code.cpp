#include "OpenCL/err_code.h"

const char *ocl_err_code(const cl_int err_in) {
  switch (err_in) {
  case CL_SUCCESS:
    return static_cast<char *>("CL_SUCCESS");
  case CL_DEVICE_NOT_FOUND:
    return static_cast<char *>("CL_DEVICE_NOT_FOUND");
  case CL_DEVICE_NOT_AVAILABLE:
    return static_cast<char *>("CL_DEVICE_NOT_AVAILABLE");
  case CL_COMPILER_NOT_AVAILABLE:
    return static_cast<char *>("CL_COMPILER_NOT_AVAILABLE");
  case CL_MEM_OBJECT_ALLOCATION_FAILURE:
    return static_cast<char *>("CL_MEM_OBJECT_ALLOCATION_FAILURE");
  case CL_OUT_OF_RESOURCES:
    return static_cast<char *>("CL_OUT_OF_RESOURCES");
  case CL_OUT_OF_HOST_MEMORY:
    return static_cast<char *>("CL_OUT_OF_HOST_MEMORY");
  case CL_PROFILING_INFO_NOT_AVAILABLE:
    return static_cast<char *>("CL_PROFILING_INFO_NOT_AVAILABLE");
  case CL_MEM_COPY_OVERLAP:
    return static_cast<char *>("CL_MEM_COPY_OVERLAP");
  case CL_IMAGE_FORMAT_MISMATCH:
    return static_cast<char *>("CL_IMAGE_FORMAT_MISMATCH");
  case CL_IMAGE_FORMAT_NOT_SUPPORTED:
    return static_cast<char *>("CL_IMAGE_FORMAT_NOT_SUPPORTED");
  case CL_BUILD_PROGRAM_FAILURE:
    return static_cast<char *>("CL_BUILD_PROGRAM_FAILURE");
  case CL_MAP_FAILURE:
    return static_cast<char *>("CL_MAP_FAILURE");
  case CL_MISALIGNED_SUB_BUFFER_OFFSET:
    return static_cast<char *>("CL_MISALIGNED_SUB_BUFFER_OFFSET");
  case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
    return static_cast<char *>("CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST");
  case CL_INVALID_VALUE:
    return static_cast<char *>("CL_INVALID_VALUE");
  case CL_INVALID_DEVICE_TYPE:
    return static_cast<char *>("CL_INVALID_DEVICE_TYPE");
  case CL_INVALID_PLATFORM:
    return static_cast<char *>("CL_INVALID_PLATFORM");
  case CL_INVALID_DEVICE:
    return static_cast<char *>("CL_INVALID_DEVICE");
  case CL_INVALID_CONTEXT:
    return static_cast<char *>("CL_INVALID_CONTEXT");
  case CL_INVALID_QUEUE_PROPERTIES:
    return static_cast<char *>("CL_INVALID_QUEUE_PROPERTIES");
  case CL_INVALID_COMMAND_QUEUE:
    return static_cast<char *>("CL_INVALID_COMMAND_QUEUE");
  case CL_INVALID_HOST_PTR:
    return static_cast<char *>("CL_INVALID_HOST_PTR");
  case CL_INVALID_MEM_OBJECT:
    return static_cast<char *>("CL_INVALID_MEM_OBJECT");
  case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
    return static_cast<char *>("CL_INVALID_IMAGE_FORMAT_DESCRIPTOR");
  case CL_INVALID_IMAGE_SIZE:
    return static_cast<char *>("CL_INVALID_IMAGE_SIZE");
  case CL_INVALID_SAMPLER:
    return static_cast<char *>("CL_INVALID_SAMPLER");
  case CL_INVALID_BINARY:
    return static_cast<char *>("CL_INVALID_BINARY");
  case CL_INVALID_BUILD_OPTIONS:
    return static_cast<char *>("CL_INVALID_BUILD_OPTIONS");
  case CL_INVALID_PROGRAM:
    return static_cast<char *>("CL_INVALID_PROGRAM");
  case CL_INVALID_PROGRAM_EXECUTABLE:
    return static_cast<char *>("CL_INVALID_PROGRAM_EXECUTABLE");
  case CL_INVALID_KERNEL_NAME:
    return static_cast<char *>("CL_INVALID_KERNEL_NAME");
  case CL_INVALID_KERNEL_DEFINITION:
    return static_cast<char *>("CL_INVALID_KERNEL_DEFINITION");
  case CL_INVALID_KERNEL:
    return static_cast<char *>("CL_INVALID_KERNEL");
  case CL_INVALID_ARG_INDEX:
    return static_cast<char *>("CL_INVALID_ARG_INDEX");
  case CL_INVALID_ARG_VALUE:
    return static_cast<char *>("CL_INVALID_ARG_VALUE");
  case CL_INVALID_ARG_SIZE:
    return static_cast<char *>("CL_INVALID_ARG_SIZE");
  case CL_INVALID_KERNEL_ARGS:
    return static_cast<char *>("CL_INVALID_KERNEL_ARGS");
  case CL_INVALID_WORK_DIMENSION:
    return static_cast<char *>("CL_INVALID_WORK_DIMENSION");
  case CL_INVALID_WORK_GROUP_SIZE:
    return static_cast<char *>("CL_INVALID_WORK_GROUP_SIZE");
  case CL_INVALID_WORK_ITEM_SIZE:
    return static_cast<char *>("CL_INVALID_WORK_ITEM_SIZE");
  case CL_INVALID_GLOBAL_OFFSET:
    return static_cast<char *>("CL_INVALID_GLOBAL_OFFSET");
  case CL_INVALID_EVENT_WAIT_LIST:
    return static_cast<char *>("CL_INVALID_EVENT_WAIT_LIST");
  case CL_INVALID_EVENT:
    return static_cast<char *>("CL_INVALID_EVENT");
  case CL_INVALID_OPERATION:
    return static_cast<char *>("CL_INVALID_OPERATION");
  case CL_INVALID_GL_OBJECT:
    return static_cast<char *>("CL_INVALID_GL_OBJECT");
  case CL_INVALID_BUFFER_SIZE:
    return static_cast<char *>("CL_INVALID_BUFFER_SIZE");
  case CL_INVALID_MIP_LEVEL:
    return static_cast<char *>("CL_INVALID_MIP_LEVEL");
  case CL_INVALID_GLOBAL_WORK_SIZE:
    return static_cast<char *>("CL_INVALID_GLOBAL_WORK_SIZE");
  case CL_INVALID_PROPERTY:
    return static_cast<char *>("CL_INVALID_PROPERTY");

  default:
    return static_cast<char *>("UNKNOWN ERROR");
  }
}

void check_ocl_error(const cl_int err, const char *operation,
                     const char *filename, const int line) {
  if (err != CL_SUCCESS) {
    fprintf(stderr, "Error during operation '%s', ", operation);
    fprintf(stderr, "in '%s' on line %d\n", filename, line);
    fprintf(stderr, "Error code was \"%s\" (%d)\n", ocl_err_code(err), err);
#if defined(_WIN32) && !defined(__MINGW32__)
    system("pause");
#endif
    exit(EXIT_FAILURE);
  }
}
