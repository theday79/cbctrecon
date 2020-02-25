/*
 * Completely rewritten from RTK Cuda forward projection filter to OpenCL
 * version by A. Gravgaard
 *
 */

#include <cassert>
#include <memory>
#include <string>

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

//#define CL_HPP_ENABLE_EXCEPTIONS

#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"
#include "OpenCL/util.hpp"

#include "OpenCL/ForwardProjectionImageFilter.hpp"

// On Linux with intel opencl runtime, you can debug opencl kernels:
// #define DEBUG_OPENCL

enum enKernelName { en_kernel_forwardProject_test, en_kernel_forwardProject };

class FwdKernelMan {

private:
  cl::Program m_program;
  bool m_initialized = false;
  bool m_defaults_is_set = false;
  std::vector<cl::Kernel> m_kernel_list;
  cl::Context m_ctx;
  cl::Device m_dev;

  static cl::Program build_ocl_program(const cl::Context &ctx,
                                       std::string &defines) {
    const auto cl_source = util::loadProgram("forward_proj.cl");

    auto err = CL_SUCCESS;
    auto program = cl::Program(ctx, cl_source, /* build? */ false, &err);
#ifdef DEBUG_OPENCL
    // For use with GDB, Requires intel OpenCL runtime
    defines += " -g -s "
               "/home/andreas/Projects/build-cbct/clang8-itk5-Debug/"
               "bin/forward_proj.cl ";
#endif
    err = program.build(defines.c_str());
#ifdef DEBUG_OPENCL
    std::cerr << "DEBUG: OpenCL defines: " << defines << "\n";
#endif
    checkError(err, "Build program");

    return program;
  }

  void make_defaults(cl_int *err) {
    if (this->m_defaults_is_set) {
      return;
    }

    auto devices = OpenCL_getDeviceList();

    // Attempt first device if none with image_support
    auto device = devices.at(0);
    auto good_device_found = false;
    auto image_devices = std::vector<cl::Device>();
    for (auto &dev : devices) {
      const auto device_image_support =
          dev.getInfo<CL_DEVICE_IMAGE_SUPPORT>(err);
      checkError(*err, "Get device image support");
      if (device_image_support) {
        image_devices.push_back(dev);
      }

      const auto device_type = dev.getInfo<CL_DEVICE_TYPE>(err);
      checkError(*err, "Get device type");

      if (device_type != CL_DEVICE_TYPE_CPU && device_image_support) {
        device = dev;
        good_device_found = true;
      }
    }
    if (!good_device_found && !image_devices.empty()) {
      device = image_devices.at(0);
    }
    auto device_name = device.getInfo<CL_DEVICE_NAME>(err);
    checkError(*err, "Get device name");
#ifdef DEBUG_OPENCL
    std::cerr << "Using " << device_name << " for fwd proj.\n";
#endif

    const auto default_ctx = cl::Context(device);
#if CL_HPP_MINIMUM_OPENCL_VERSION >= 200
    auto deviceQueue =
        cl::DeviceCommandQueue::makeDefault(default_ctx, device, err);
#else
    cl::Device::setDefault(device);
    cl::Context::setDefault(default_ctx);
    cl::CommandQueue::setDefault(cl::CommandQueue(default_ctx, device, 0, err));
#endif
    checkError(*err, "Set default queue");

    this->m_defaults_is_set = true;
  }

public:
  FwdKernelMan() = default;
  void initialize(std::string &defines) {
    auto err = CL_SUCCESS;
    this->make_defaults(&err);
    m_dev = cl::Device::getDefault(&err);
    checkError(err, "Get default device");
    m_ctx = cl::Context::getDefault(&err);
    checkError(err, "Get default context");
    m_program = build_ocl_program(m_ctx, defines);
    m_program.createKernels(&m_kernel_list);
    m_initialized = true;
  }

  cl::Kernel getKernel(const enKernelName kernel) {
    assert(m_initialized);

    auto &out_kernel = m_kernel_list.at(kernel);
    const auto kernel_name = out_kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
#ifdef DEBUG_OPENCL
    std::cerr << kernel_name << "\n";
#endif
    return out_kernel;
  }
};

// static KernelMan kernel_man;

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

std::string
make_OpenCL_defines_str(const OpenCL_forwardProject_options &fwd_opts) {
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
    // TODO: Check vectorlength at a higher level so this is impossible to reach
    std::cerr << "This vectorlength is not implemented, you'll get "
                 "CL_PROGRAM_BUILD_FAILURE:\n";
    break; // Will fail with -11 CL_PROGRAM_BUILD_FAILURE
  }

  return cl_defines;
}

void OpenCL_forward_project(float *h_proj_in, float *h_proj_out, float *h_vol,
                            OpenCL_forwardProject_options &fwd_opts) {

  auto err = CL_SUCCESS;
  // Constant memory
  auto cl_defines = make_OpenCL_defines_str(fwd_opts);

  auto kernel_man = std::make_unique<FwdKernelMan>();
  kernel_man->initialize(cl_defines);

  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  const auto tot_proj_size = fwd_opts.projSize.at(0) * fwd_opts.projSize.at(1) *
                             fwd_opts.projSize.at(2);

  const auto tot_vol_size =
      fwd_opts.volSize.at(0) * fwd_opts.volSize.at(1) * fwd_opts.volSize.at(2);

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
                        cl::Buffer, cl::Buffer>(
          kernel_man->getKernel(en_kernel_forwardProject));
  checkError(err, "Create kernel functor");

  const auto req_dev_alloc = tot_proj_size * sizeof(float);
  const auto avail_dev_alloc =
      cl::Device::getDefault().getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>(&err);
  checkError(err, "Get device max alloc size");
  if (avail_dev_alloc < req_dev_alloc) {
    std::cerr << "Oh no, device doesn't have enough memory!\n"
              << "It has: " << avail_dev_alloc / 1024 / 1024 << "MB\n"
              << "But: " << req_dev_alloc / 1024 / 1024 << "MB was required!\n";
  }

  // Create an array of textures
  auto channel_order = CL_INTENSITY;
  if (fwd_opts.vectorLength == 3) {
    channel_order = CL_RGB;
  }
  const auto dev_tex_vol = cl::Image3D(
      ctx, CL_MEM_READ_ONLY, cl::ImageFormat(channel_order, CL_FLOAT),
      fwd_opts.volSize[0], fwd_opts.volSize[1], fwd_opts.volSize[2],
      /*row_pitch*/ 0, /*slice_pitch*/ 0, nullptr, &err);
  checkError(err, "Create 3D Image texture");

  // Output projection = input + forward proj
  const cl::Buffer dev_proj_out(ctx, CL_MEM_READ_WRITE, req_dev_alloc);
  checkError(err, "Alloc proj_out on device");

  // Input projection (add to this)
  const auto dev_proj_in =
      cl::Buffer(queue, h_proj_in, h_proj_in + tot_proj_size, true, true, &err);
  checkError(err, "Alloc proj_in on device");

  // Volume to forward project
  const auto dev_vol =
      cl::Buffer(queue, h_vol, h_vol + tot_vol_size, true, true, &err);
  checkError(err, "Alloc vol on device");

  // Write image data to texture object
  const std::array<cl::size_type, 3> origin = {{0, 0, 0}};
  const std::array<cl::size_type, 3> region = {
      {static_cast<cl::size_type>(fwd_opts.volSize.at(0)),
       static_cast<cl::size_type>(fwd_opts.volSize.at(1)),
       static_cast<cl::size_type>(fwd_opts.volSize.at(2))}};
  err = queue.enqueueWriteImage(dev_tex_vol, CL_TRUE, origin, region, 0, 0,
                                h_vol);
  checkError(err, "Copy 3D Image to device");

  // Write geometries to constant memory
  const auto dev_trn_prj_idx_trf_mats = cl::Buffer(
      queue, fwd_opts.translatedProjectionIndexTransformMatrices.begin(),
      fwd_opts.translatedProjectionIndexTransformMatrices.end(), true, true,
      &err);
  checkError(err, "Copy trnslProjIndexMats to device");

  const auto dev_source_positions =
      cl::Buffer(queue, fwd_opts.source_positions.begin(),
                 fwd_opts.source_positions.end(), true, true, &err);
  checkError(err, "Copy source_positions to device");

  const auto local_workgroup_size = cl::NDRange(16, 16);
  const auto global_workgroup_size =
      cl::NDRange(fwd_opts.projSize[0], fwd_opts.projSize[1]);
  // Run kernel
  kernel_forwardProject(
      cl::EnqueueArgs(queue, global_workgroup_size, local_workgroup_size),
      dev_proj_out, dev_proj_in, dev_vol, dev_tex_vol, dev_trn_prj_idx_trf_mats,
      dev_source_positions, err);
  checkError(err, "Forward kernel, enqueue");

  err = queue.finish();
  checkError(err, "Forward kernel, finish queue");

  // Read back data from device to host
  err = queue.enqueueReadBuffer(dev_proj_out, CL_TRUE, 0, req_dev_alloc,
                                h_proj_out);
  checkError(err, "Read out buffer from device");

  err = queue.finish();
  checkError(err, "Forward kernel, finish queue");
}
