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

#ifndef RTKOPENCLUTILITIES_H
#define RTKOPENCLUTILITIES_H

#ifdef RTK_USE_OPENCL

#ifndef CBCTRECON_OPENCL_VERSION
#define CBCTRECON_OPENCL_VERSION 120
#endif

#if CBCTRECON_OPENCL_VERSION >= 210
#define CL_HPP_MINIMUM_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 210
#else
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#endif

#include <OpenCL/cl2.hpp>

#include <string>
#include <vector>

#include "itkMacro.h" // itkGenericExceptionMacro

/** \brief Macro to check errors when running an OpenCL command.
 *
 * \author Simon Rit
 *
 * \ingroup Macro
 */
#define OPENCL_CHECK_ERROR(cmd)                                                \
  {                                                                            \
    cl_int rc = cmd;                                                           \
    if (rc != CL_SUCCESS)                                                      \
      itkGenericExceptionMacro(<< "OPENCL ERROR with " << #cmd                 \
                               << ". Returned value is " << rc);               \
  }

/** \brief Get the list of OpenCL compatible platforms
 *
 * \author Simon Rit
 *
 * \ingroup Functions
 */
std::vector<cl::Platform> GetListOfOpenCLPlatforms();

/** \brief Get the list of OpenCL compatible devices
 *
 * \author Simon Rit
 *
 * \ingroup Functions
 */
std::vector<cl::Device> GetListOfOpenCLDevices(cl::Platform platform);

/** \brief Builds an OpenCL program in a given filename given a context
 *
 * \author Simon Rit
 *
 * \ingroup Functions
 */
void CreateAndBuildOpenCLProgramFromSourceFile(std::string &fileName,
                                               const cl::Context &context,
                                               cl::Program &program);

#endif
#endif
