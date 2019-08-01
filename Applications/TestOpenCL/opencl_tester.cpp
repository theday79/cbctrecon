#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "OpenCLFFTFilter.h"

int main(int argc, char** argv){
  auto plat_dev = getPlatformAndDeviceID(0);

  {
    // Print some platform info:
    const auto plat = std::get<0>(plat_dev);
    char plat_name[128];
    auto cl_err = clGetPlatformInfo(plat, CL_PLATFORM_NAME,
            sizeof(plat_name), &plat_name, nullptr);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get platform name, CL_ERR: " << cl_err << "\n";
      return -1;
    }
    std::cerr << "Using platform: " << plat_name << "\n";
  }

  {
    // Print some device info:
    const auto dev = std::get<1>(plat_dev);
    char device_name[128];
    auto cl_err = clGetDeviceInfo(dev, CL_DEVICE_NAME,
            sizeof(device_name), &device_name, nullptr);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get device name, CL_ERR: " << cl_err << "\n";
      return -2;
    }
    std::cerr << "Using device: " << device_name << "\n";

    // Image support?
    cl_bool device_image_support = false;
    cl_err = clGetDeviceInfo(dev, CL_DEVICE_IMAGE_SUPPORT,
            sizeof(cl_bool), &device_image_support, nullptr);
    if (cl_err != CL_SUCCESS){
        std::cerr << "Could not get device image support, CL_ERR: " << cl_err << "\n";
    }
    if (device_image_support){
      std::cerr << "Device has image support :)\n";
    } else {
      std::cerr << "Device does NOT support opencl image types! :(\n";
    }
  }

  return 0;
}
