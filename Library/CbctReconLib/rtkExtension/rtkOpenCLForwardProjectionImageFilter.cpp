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
/*
 * Completely rewritten from Cuda version to OpenCL
 * by A. Gravgaard
 *
 */

/*****************
 *  rtk #includes *
 *****************/
#include "rtkConfiguration.h"

#include <string>

#if CBCTRECON_OPENCL_VERSION >= 210
#define CL_HPP_MINIMUM_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 210
#else
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#endif

// #define CL_HPP_ENABLE_EXCEPTIONS

#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"
#include "OpenCL/util.hpp"

#include "rtkOpenCLForwardProjectionImageFilter.h"

template <typename T> std::string stringify_array(const std::vector<T> &arr) {
  std::string out;
  for (unsigned int i = 0; i < arr.size() - 1; ++i) {
    out += std::to_string(arr.at(i)) + ",";
  }
  out += std::to_string(arr.at(arr.size() - 1)) + " ";

  return out;
}

template <typename T, unsigned int N>
std::string stringify_array(const std::array<T, N> &arr) {
  std::string out;
  for (unsigned int i = 0; i < N - 1; ++i) {
    out += std::to_string(arr.at(i)) + ",";
  }
  out += std::to_string(arr.at(N - 1)) + " ";

  return out;
}

std::string make_OpenCL_defines_str(OpenCL_forwardProject_options &fwd_opts) {
  auto cl_defines = std::string(" ");

  cl_defines += std::string("-DC_PROJ_SIZE=") +
                stringify_array<unsigned int, 3>(fwd_opts.projSize);

  cl_defines +=
      std::string("-DC_BOX_MIN=") + stringify_array<float, 3>(fwd_opts.box_min);

  cl_defines +=
      std::string("-DC_BOX_MAX=") + stringify_array<float, 3>(fwd_opts.box_max);

  cl_defines +=
      std::string("-DC_SPACING=") + stringify_array<float, 3>(fwd_opts.spacing);

  cl_defines += std::string("-DC_VOL_SIZE=") +
                stringify_array<unsigned int, 3>(fwd_opts.volSize);

  cl_defines +=
      std::string("-DC_TSTEP=") + std::to_string(fwd_opts.t_step) + " ";

  cl_defines += std::string("-DSLAB_SIZE=") +
                std::to_string(fwd_opts.projSize.at(2)) + " ";

  /*
  cl_defines += std::string("-DC_TRANS_PROJ_IDX_TRN_MAT=") +
                stringify_array<float>(
                    fwd_opts.translatedProjectionIndexTransformMatrices);
                    */

  if (fwd_opts.radiusCylindricalDetector > 0.00001) {
    cl_defines += std::string("-DC_RADIUS=") +
                  std::to_string(fwd_opts.radiusCylindricalDetector) + " ";

    cl_defines +=
        std::string("-DC_TRANS_VOL_TRN_MAT=") +
        stringify_array<float>(fwd_opts.translatedVolumeTransformMatrices);
  } else {
    cl_defines += std::string("-DC_RADIUS_ZERO=true ");
  }

  /*
  cl_defines += std::string("-DC_SOURCE_POS=") +
                stringify_array<float>(fwd_opts.source_positions);
                */

  switch (fwd_opts.vectorLength) {
  case 1:
    cl_defines += "-DVECTOR_LENGTH=1 ";
    break;
  case 3:
    cl_defines += "-DVECTOR_LENGTH=3 ";
    break;
  default:
    itkGenericExceptionMacro("Vector length " << fwd_opts.vectorLength
                                              << " is not supported.")
  }

  return cl_defines;
}

void OpenCL_forward_project(float *h_proj_in, float *h_proj_out, float *h_vol,
                            OpenCL_forwardProject_options &fwd_opts) {

  std::vector<cl::Device> devices;
  OpenCL_getDeviceList(devices);

  auto err = CL_SUCCESS;
  // Attempt first device if none with image_support
  cl::Device device = devices.at(5);
  /*
  for (auto &dev : devices) {
    auto device_image_support = dev.getInfo<CL_DEVICE_IMAGE_SUPPORT>(&err);
    checkError(err, "Get device image support");
    if (device_image_support) {
      device = dev;
      break;
    }
  }*/
  auto ctx = cl::Context(device);
  auto queue = cl::CommandQueue(ctx);

  auto tot_proj_size = fwd_opts.projSize.at(0) * fwd_opts.projSize.at(1) *
                       fwd_opts.projSize.at(2);

  // Maybe do the pinned memory trick?

  auto tot_vol_size =
      fwd_opts.volSize.at(0) * fwd_opts.volSize.at(1) * fwd_opts.volSize.at(2);
  if (!h_vol) {
    std::cerr << "h_vol ptr is invalid\n";
  }
  if (!(h_vol + tot_vol_size - 1)) {
    std::cerr << "last index of h_vol is invalid\n";
  }

  auto program = cl::Program(util::loadProgram("forward_proj.cl"), false, &err);
  checkError(err, "Load forward_proj.cl");

  // Constant memory
  auto cl_defines = make_OpenCL_defines_str(fwd_opts);
  std::cerr << "DEBUG: OpenCL defines: " << cl_defines << "\n";

  err = program.build(cl_defines.c_str());
  checkError(err, "Build forward_proj.cl");

  /*
   * __kernel void kernel_forwardProject(
   *  __global float *dev_proj_out,
   *  __global float *dev_proj_in,
   *  __global float *dev_vol,
   *  __read_only image3d_t tex_vol,
   *  __constant float *c_translatedProjectionIndexTransformMatrices,
   *  __constant float *c_sourcePos)
   */
  auto kernel_forwardProject =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer, cl::Image3D,
                        cl::Buffer, cl::Buffer>(program,
                                                "kernel_forwardProject", &err);
  checkError(err, "Create kernel  functor");

  auto req_dev_alloc = tot_proj_size * sizeof(float);
  auto avail_dev_alloc = device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
  if (avail_dev_alloc < req_dev_alloc) {
    std::cerr << "Oh no, device doesn't have enough memory!\n"
              << "It has: " << avail_dev_alloc / 1024 / 1024 << "MB\n"
              << "But: " << req_dev_alloc / 1024 / 1024 << "MB was required!\n";
  }

  // Output projection = input + forward proj
  cl::Buffer dev_proj_out(ctx, CL_MEM_WRITE_ONLY, sizeof(float) * tot_proj_size,
                          nullptr, &err);
  checkError(err, "Alloc proj_out on device");

  // Input projection (add to this)
  cl::Buffer dev_proj_in =
      cl::Buffer(ctx, h_proj_in, h_proj_in + tot_proj_size, true, true, &err);
  checkError(err, "Alloc proj_in on device");

  // Volume to forward project
  cl::Buffer dev_vol(ctx, h_vol, h_vol + tot_vol_size, true, true, &err);
  checkError(err, "Alloc vol on device");

  // Create an array of textures
  auto dev_tex_vol = cl::Image3D(
      ctx, CL_MEM_READ_ONLY, cl::ImageFormat(CL_INTENSITY, CL_FLOAT),
      fwd_opts.volSize[0], fwd_opts.volSize[1], fwd_opts.volSize[2],
      /*row_pitch*/ 0, /*slice_pitch*/ 0, nullptr, &err);
  checkError(err, "Create 3D Image texture");

  std::array<cl::size_type, 3> origin = {{0, 0, 0}};
  std::array<cl::size_type, 3> region = {
      {static_cast<cl::size_type>(fwd_opts.volSize.at(0)),
       static_cast<cl::size_type>(fwd_opts.volSize.at(1)),
       static_cast<cl::size_type>(fwd_opts.volSize.at(2))}};
  err = queue.enqueueWriteImage(dev_tex_vol, CL_TRUE, origin, region, 0, 0,
                                h_vol);
  checkError(err, "Copy 3D Image to device");

  auto dev_trn_prj_idx_trf_mats = cl::Buffer(
      ctx, fwd_opts.translatedProjectionIndexTransformMatrices.begin(),
      fwd_opts.translatedProjectionIndexTransformMatrices.end(), true, true,
      &err);
  checkError(err, "Copy trnslProjIndexMats to device");

  auto dev_source_positions =
      cl::Buffer(ctx, fwd_opts.source_positions.begin(),
                 fwd_opts.source_positions.end(), true, true, &err);
  checkError(err, "Copy source_positions to device");

  auto local_workgroup_size = cl::NDRange(16, 16);
  auto global_workgroup_size =
      cl::NDRange(fwd_opts.projSize[0], fwd_opts.projSize[1]);

  kernel_forwardProject(
      cl::EnqueueArgs(queue, global_workgroup_size, local_workgroup_size),
      dev_proj_out, dev_proj_in, dev_vol, dev_tex_vol, dev_trn_prj_idx_trf_mats,
      dev_source_positions, err);
  checkError(err, "Forward kernel, enqueue");

  err = queue.finish();
  checkError(err, "Forward kernel, finish queue");

  // Read back data from device to host
  err = cl::copy(queue, dev_proj_out, h_proj_out, h_proj_out + tot_proj_size);
  checkError(err, "Copy proj back from device");

  auto out_sum = 0.0;
  for (auto i = 0U; i < tot_proj_size; ++i) {
    out_sum += h_proj_out[i];
  }
  std::cerr << "DEBUG: h_proj_out test sum over " << tot_proj_size
            << " elements: " << out_sum << "\n";
}
