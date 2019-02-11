// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

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
#include <vector>
#include <fstream>

std::vector<cl_platform_id> GetListOfOpenCLPlatforms() {
  cl_uint numberOfPlatforms;
  OPENCL_CHECK_ERROR(clGetPlatformIDs(0, nullptr, &numberOfPlatforms));

  std::vector<cl_platform_id> platformList(numberOfPlatforms);
  OPENCL_CHECK_ERROR(
      clGetPlatformIDs(numberOfPlatforms, &platformList[0], nullptr));

  return platformList;
}

std::vector<cl_device_id> GetListOfOpenCLDevices(cl_platform_id platform) {
  cl_uint numberOfDevices;
  cl_int error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr,
                                &numberOfDevices);

  std::vector<cl_device_id> deviceList(numberOfDevices);
  if (error != -1){ // -1 means we didn't find a GPU
    OPENCL_CHECK_ERROR(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU,
                                      numberOfDevices, &deviceList[0], nullptr));
  }

  cl_bool bImageSupport = 0u;
  // If found, check if supports image.
  if (numberOfDevices > 0)
    OPENCL_CHECK_ERROR(clGetDeviceInfo(deviceList[0], CL_DEVICE_IMAGE_SUPPORT,
                                       sizeof(cl_bool), &bImageSupport,
                                       nullptr));

  // If not a good device, switch to CPU.
  if (bImageSupport == 0u) {
    OPENCL_CHECK_ERROR(clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, nullptr,
                                      &numberOfDevices));
    deviceList.resize(numberOfDevices);
    OPENCL_CHECK_ERROR(clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU,
                                      numberOfDevices, &deviceList[0],
                                      nullptr));
  }
  /*
    char buf[1024];
    OPENCL_CHECK_ERROR( clGetDeviceInfo (deviceList[0], CL_DEVICE_NAME,
    sizeof(buf), buf, NULL) ); std::cout << "Buf=" << buf << std::endl;
  */
  return deviceList;
}

void CreateAndBuildOpenCLProgramFromSourceFile(std::string &fileName,
                                               const cl_context &context,
                                               cl_program &program) {
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

  program = clCreateProgramWithSource(
      context, 1, const_cast<const char **>(&oclSource), &size, &error);
  if (error != CL_SUCCESS)
    itkGenericExceptionMacro(
        << "Could not create OpenCL sampler object, error code: " << error);

  error = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (error != CL_SUCCESS) {
    cl_device_id id;
    clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(cl_device_id), &id,
                     nullptr);

    size_t logSize;
    OPENCL_CHECK_ERROR(clGetProgramBuildInfo(program, id, CL_PROGRAM_BUILD_LOG,
                                             0, nullptr, &logSize));
    // auto *log = new char[logSize + 1];
    auto log = std::vector<char>(logSize + 1);
    OPENCL_CHECK_ERROR(clGetProgramBuildInfo(program, id, CL_PROGRAM_BUILD_LOG,
                                             logSize, &log[0], &logSize));
    log[logSize] = '\0';
    itkGenericExceptionMacro(<< "OPENCL ERROR with clBuildProgram. The log is:\n"
                             << &log[0]);
  }
  delete[] oclSource;
}
