// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "OpenCL/device_picker.hpp"

#include <vector>
#include <string>

namespace crl {
namespace opencl {

std::vector<cl::Device> getDeviceList() {
  // Get list of platforms
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  std::vector<cl::Device> devices;
  // Enumerate devices
  for (auto &platform : platforms) {
    std::vector<cl::Device> plat_devices;
    const auto err = platform.getDevices(CL_DEVICE_TYPE_ALL, &plat_devices);
    checkError(err, "Get devices from platform");
    std::copy(plat_devices.begin(), plat_devices.end(),
              std::back_inserter(devices));
  }

  return devices;
}

std::string getDeviceName(const cl::Device &device) {
  std::string name;
  cl_device_info info = CL_DEVICE_NAME;

  // Special case for AMD
#ifdef CL_DEVICE_BOARD_NAME_AMD
  name = device.getInfo<CL_DEVICE_VENDOR>();
  if (strstr(name.c_str(), "Advanced Micro Devices"))
    info = CL_DEVICE_BOARD_NAME_AMD;
#endif

  const auto err = device.getInfo(info, &name);
  checkError(err, "Get device name");
  return name;
}

} // namespace opencl
} // namespace crl
