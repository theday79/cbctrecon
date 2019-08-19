#include "OpenCL/device_picker.hpp"

unsigned OpenCL_getDeviceList(std::vector<cl::Device> &devices) {
  // Get list of platforms
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  // Enumerate devices
  for (auto &platform : platforms) {
    std::vector<cl::Device> plat_devices;
    const auto err = platform.getDevices(CL_DEVICE_TYPE_ALL, &plat_devices);
    checkError(err, "Get devices from platform");
    devices.insert(devices.end(), plat_devices.begin(), plat_devices.end());
  }

  return devices.size();
}

std::string OpenCL_getDeviceName(const cl::Device &device) {
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

template <typename T> T string_to(const char *str, char *next);

template <> cl_uint string_to<cl_uint>(const char *str, char *next) {
  const auto output = std::strtoul(str, &next, 10);
  return static_cast<cl_uint>(output);
}

template <typename T> int OpenCL_parse(const char *str, T *output) {
  char *next = nullptr;
  *output = string_to<T>(str, next);
  return !strlen(next);
}

void OpenCL_parseArguments(const int argc, char *argv[], cl_uint *deviceIndex) {
  for (auto i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--list")) {
      // Get list of devices
      std::vector<cl::Device> devices;
      const auto numDevices = OpenCL_getDeviceList(devices);

      // Print device names
      if (numDevices == 0) {
        std::cout << "No devices found.\n";
      } else {
        std::cout << "\nDevices:\n";
        for (unsigned int j = 0; j < numDevices; j++) {
          std::cout << j << ": " << OpenCL_getDeviceName(devices[j]) << "\n";
        }
        std::cout << "\n";
      }
      exit(0);
    }
    if (!strcmp(argv[i], "--device")) {
      if (++i >= argc || !OpenCL_parse<cl_uint>(argv[i], deviceIndex)) {
        std::cout << "Invalid device index\n";
        exit(1);
      }
    } else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      std::cout << "\n";
      std::cout << "Usage: ./program [OPTIONS]\n\n";
      std::cout << "Options:\n";
      std::cout << "  -h  --help               Print the message\n";
      std::cout << "      --list               List available devices\n";
      std::cout << "      --device     INDEX   Select device at INDEX\n";
      std::cout << "\n";
      exit(0);
    }
  }
}
