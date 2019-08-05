// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "OpenCL/ImageFilters.h"
#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"

int main(int argc, char** argv){

  // Get list of platforms
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);


  auto i = 0U;
  for (auto &plat : platforms) {
    // Print some platform info:
    std::string plat_name;
    auto cl_err = plat.getInfo(CL_PLATFORM_NAME, &plat_name);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get platform name, CL_ERR: " << cl_err << "\n";
      return -1;
    }
    std::cerr << "Platform #" << ++i << ": " << plat_name << "\n";
  }

  auto devices = std::vector<cl::Device>();
  getDeviceList(devices);

  i = 0U;
  for (auto &dev : devices) {
    // Print some device info:
    std::string device_name;
    auto cl_err = dev.getInfo(CL_DEVICE_NAME, &device_name);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get device name, CL_ERR: " << cl_err << "\n";
      return -2;
    }
    std::cerr << "Device #" << ++i << ": " << device_name << "\n";

    // Image support?
    cl_bool device_image_support = false;
    cl_err = dev.getInfo(CL_DEVICE_IMAGE_SUPPORT, &device_image_support);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get device image support, CL_ERR: " << cl_err << "\n";
      return -3;
    }
    if (device_image_support){
      std::cerr << "Device has image support :)\n";
    } else {
      std::cerr << "Device does NOT support opencl image types! :(\n";
    }

    // Memory?:
    cl_ulong max_mem_alloc;
    cl_err = dev.getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &max_mem_alloc);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get device memory max alloc, CL_ERR: " << cl_err << "\n";
      return -4;
    }
    std::cerr << "Max memory alloc: " << max_mem_alloc / (1024 * 1024) << " MB\n";

    cl_ulong global_mem;
    cl_err = dev.getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &global_mem);
    if (cl_err != CL_SUCCESS){
      std::cerr << "Could not get device global memory, CL_ERR: " << cl_err << "\n";
      return -5;
    }
    std::cerr << "Global memory: " << global_mem / (1024 * 1024) << " MB\n";
  }

  return 0;
}
