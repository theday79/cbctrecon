// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "rtkOpenCLUtilities.h"
#include <fstream>
#include <vector>

std::vector<cl::Platform> GetListOfOpenCLPlatforms() {
  // Get list of platforms
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  return platforms;
}

std::vector<cl::Device> GetListOfOpenCLDevices(cl::Platform platform) {
  // If no devices support images, returns empty vector

  // Enumerate devices
  auto devices = std::vector<cl::Device>();
  platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

  cl_bool bImageSupport = false;
  // If found, check if supports image.
  auto dev_it =
      std::remove_if(devices.begin(), devices.end(), [](cl::Device device) {
        return !(device.getInfo<CL_DEVICE_IMAGE_SUPPORT>());
      });
  devices.erase(dev_it, devices.end());

  bImageSupport = dev_it != devices.end();

  // If not a good device, switch to Accelerator.
  if (!bImageSupport) {
    devices.clear();
    platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);

    // If found, check if supports image.
    auto dev_it =
        std::remove_if(devices.begin(), devices.end(), [](cl::Device device) {
          return !(device.getInfo<CL_DEVICE_IMAGE_SUPPORT>());
        });
    devices.erase(dev_it, devices.end());

    bImageSupport = dev_it != devices.end();
  }

  // If still not a good device, switch to CPU.
  if (!bImageSupport) {
    devices.clear();
    platform.getDevices(CL_DEVICE_TYPE_CPU, &devices);

    // If found, check if supports image.
    auto dev_it =
        std::remove_if(devices.begin(), devices.end(), [](cl::Device device) {
          return !(device.getInfo<CL_DEVICE_IMAGE_SUPPORT>());
        });
    devices.erase(dev_it, devices.end());

    bImageSupport = dev_it != devices.end();
  }

  return devices;
}

void CreateAndBuildOpenCLProgramFromSourceFile(std::string &fileName,
                                               const cl::Context &context,
                                               cl::Program &program) {
  char *oclSource;
  size_t size;
  cl_int error;

  // Open file stream
  // const std::string& completeFileName(fileName);
  // //std::string(RTK_BINARY_DIR) + std::string("/") +
  std::fstream f(fileName.c_str(), (std::fstream::in | std::fstream::binary));

  // Check if we have opened file stream
  if (f.is_open()) {
    // Find the stream size
    f.seekg(0, std::fstream::end);
    size = static_cast<size_t>(f.tellg());
    f.seekg(0, std::fstream::beg);

    oclSource = new char[size + 1];

    // Read file
    f.read(oclSource, size);
    f.close();
    oclSource[size] = '\0';
  } else {
    itkGenericExceptionMacro(<< "Could not read OpenCL source file "
                             << fileName);
  }

  program = cl::Program(context, oclSource, false, &error);
  // (    context, 1, const_cast<const char **>(&oclSource), &size, &error);
  if (error != CL_SUCCESS)
    itkGenericExceptionMacro(
        << "Could not create OpenCL sampler object, error code: " << error);

  error = program.build();
  // error = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (error != CL_SUCCESS) {
    auto id = context.getInfo<CL_CONTEXT_DEVICES>();

    auto log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(&error);

    itkGenericExceptionMacro(
        << "OPENCL ERROR with clBuildProgram. The log is:\n"
        << &log[0]);
  }
  delete[] oclSource;
}
