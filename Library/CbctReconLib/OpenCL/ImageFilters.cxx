// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <cassert>

#include "OpenCL/ImageFilters.h"

#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"
#include "OpenCL/util.hpp"

// For crop_by_structure
#include "StructureSet.h"

// On Linux with intel opencl runtime, you can debug opencl kernels:
// #define DEBUG_OPENCL

enum enKernel {
  en_padding_kernel,
  en_multiply_kernel,
  en_multiply_kernel2D,
  en_subtract_kernel2D,
  en_divide_kernel3D_ushort,
  en_add_const_kernel,
  en_add_const_with_thresh_kernel,
  en_add_mul_const_kernel,
  en_add_mul_const_with_thresh_kernel,
  en_min_max_kernel,
  en_min_max_kernel2,
  en_crop_by_struct_kernel
};

class KernelMan {

private:
  cl::Program m_program;
  bool m_initialized = false;
  bool m_defaults_is_set = false;
  std::vector<cl::Kernel> m_kernel_list;
  cl::Context m_ctx;
  cl::Device m_dev;

  static cl::Program build_ocl_program(const cl::Context &ctx,
                                       std::string &defines) {
#if CL_HPP_TARGET_OPENCL_VERSION >= 210 ||                                     \
    (CL_HPP_TARGET_OPENCL_VERSION == 200 && defined(CL_HPP_USE_IL_KHR))
    // const auto cl_source = util::loadProgramIL("kernels_cl20.spv");
    const auto cl_source = util::loadProgram("filters.cl");
#else
    const auto cl_source = util::loadProgram("filters.cl");
#endif

    auto err = CL_SUCCESS;
    auto program = cl::Program(ctx, cl_source, /* build? */ false, &err);
#if defined(DEBUG_OPENCL) && !defined(_WIN32)
    // For use with GDB, Requires intel OpenCL runtime
    defines += " -g -s "
               "/home/andreas/Projects/build-cbct/clang8-itk5-Debug/"
               "bin/filters.cl ";
#endif
    err = program.build(defines.c_str());
#ifdef DEBUG_OPENCL
    std::cerr << "DEBUG: OpenCL defines: " << defines << "\n";
#endif
    checkError(err, "Build program");

    return program;
  }

  void make_defaults(cl_int *err, const size_t req_dev_alloc) {
    if (this->m_defaults_is_set) {
      const auto dev_max_alloc =
          cl::Device::getDefault().getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>(err);
      checkError(*err, "Get dev max alloc pre-check");
      if (dev_max_alloc > req_dev_alloc) {
        return;
      }
    }

    std::vector<cl::Device> devices;
    OpenCL_getDeviceList(devices);

    // Attempt first device if none with image_support
    auto device = devices.at(0);
    auto good_device_found = false;
    auto image_devices = std::vector<cl::Device>();
    for (auto &dev : devices) {
      const auto avail_dev_alloc =
          dev.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>(err);
      checkError(*err, "Get device max alloc");
      if (avail_dev_alloc > req_dev_alloc) {
        image_devices.push_back(dev);
      }

      const auto device_type = dev.getInfo<CL_DEVICE_TYPE>(err);
      checkError(*err, "Get device type");

      if (device_type != CL_DEVICE_TYPE_CPU &&
          avail_dev_alloc > req_dev_alloc) {
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
    std::cerr << "Using " << device_name << "\n";
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
  KernelMan() = default;
  void initialize(std::string &defines, const size_t req_dev_alloc = 0) {
    auto err = CL_SUCCESS;
    this->make_defaults(&err, req_dev_alloc);
    m_dev = cl::Device::getDefault(&err);
    checkError(err, "Get default device");
    m_ctx = cl::Context::getDefault(&err);
    checkError(err, "Get default context");
    m_program = build_ocl_program(m_ctx, defines);
    m_program.createKernels(&m_kernel_list);
    m_initialized = true;
  }

  cl::Kernel getKernel(const enKernel kernel) const {
    assert(m_initialized);

#ifdef _WIN32
    // Maybe only necessary for xeon phi
    const auto get_kernel_name = [](enKernel en_kernel) {
      switch (en_kernel) {
      case en_padding_kernel:
        return "padding_kernel";
      case en_multiply_kernel:
        return "multiply_kernel";
      case en_multiply_kernel2D:
        return "multiply_kernel2D";
      case en_subtract_kernel2D:
        return "subtract_kernel2D";
      case en_divide_kernel3D_ushort:
        return "divide_kernel3D_Ushort";
      case en_add_const_kernel:
        return "add_const_kernel";
      case en_add_const_with_thresh_kernel:
        return "add_const_with_thresh_kernel";
      case en_add_mul_const_kernel:
        return "add_mul_const_kernel";
      case en_add_mul_const_with_thresh_kernel:
        return "add_mul_const_with_thresh_kernel";
      case en_min_max_kernel:
        return "min_max_kernel";
      case en_min_max_kernel2:
        return "min_max_kernel2";
      case en_crop_by_struct_kernel:
        return "crop_by_struct_kernel";
      };
    };

    const auto get_kernel_by_name =
        [](const std::string &kernel_name,
           const std::vector<cl::Kernel> &kernel_list) {
          for (auto &cur_kernel : kernel_list) {
            const auto cur_kernel_name =
                cur_kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
            if (cur_kernel_name == kernel_name) {
              return cur_kernel;
            }
          }
          std::cerr << "OH NO, SOMETHING BAD WILL HAPPEN SOON!\n";
          return kernel_list.at(0);
        };

    const auto &kernel_name = get_kernel_name(kernel);
    auto out_kernel = get_kernel_by_name(kernel_name, m_kernel_list);
#else
    // It seems linux always has the kernels in order, probably it's just
    // windows which doesn't
    auto &out_kernel = m_kernel_list.at(kernel);
    const auto &kernel_name = out_kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
#endif

#ifdef DEBUG_OPENCL
    std::cerr << kernel_name << "\n";
#endif
    return out_kernel;
  }
};

static KernelMan kernel_man;

void OpenCL_initialize(const size_t required_mem_alloc_size,
                       std::string &defines) {
  kernel_man.initialize(defines, required_mem_alloc_size);
}

void OpenCL_initialize(std::string &defines) { kernel_man.initialize(defines); }

cl::NDRange get_local_work_size_small(const cl::Device &device) {
  const auto device_name = device.getInfo<CL_DEVICE_NAME>();
  if (device_name == "Intel(R) Iris(TM) Pro Graphics 5200") {
    return {16};
  }
  return {128};
}
cl::NDRange get_local_work_size_large(const cl::Device &device) {
  const auto device_name = device.getInfo<CL_DEVICE_NAME>();
  if (device_name == "Intel(R) Iris(TM) Pro Graphics 5200") {
    return {32};
  }
  return {128};
}

void OpenCL_padding(const cl_int4 &paddingIndex, const cl_uint4 &paddingSize,
                    const cl_uint4 &inputSize,
                    const float *hostVolume, // input
                    float *hostPaddedVolume, // output
                    std::vector<float> &mirrorWeights) {
  auto err = CL_SUCCESS;

  const auto pv_size =
      static_cast<size_t>(paddingSize.x * paddingSize.y * paddingSize.z);
  const auto pv_buffer_size = pv_size * sizeof(float);

  auto defines = std::string("");
  kernel_man.initialize(defines, pv_buffer_size);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */

  // just because intellisense couldn't understand it below...
  const auto w_buf_sizeof = static_cast<cl_uint>(mirrorWeights.size());

  const auto weights_d =
      cl::Buffer(ctx, mirrorWeights.begin(), mirrorWeights.end(),
                 /*read only*/ true, /*use host ptr*/ false, &err);
  checkError(err, "Pad, allocate weights on device");

  const auto devicePaddedVolume =
      cl::Buffer(ctx, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, pv_buffer_size,
                 (void *)&hostPaddedVolume[0], &err);
  checkError(err, "Pad, Allocate device padded volume");

  const auto v_buffer_size = sizeof(hostVolume);
  const auto deviceVolume =
      cl::Buffer(ctx, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, v_buffer_size,
                 (void *)&hostVolume[0], &err);
  checkError(err, "Pad, allocate volume device buffer");

  auto padding_kernel =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_int4, cl_uint4, cl_uint4,
                        cl::Buffer, cl_uint>(
          kernel_man.getKernel(en_padding_kernel));
  // padding_kernel(program, "padding_kernel");

  const auto local_work_size =
      get_local_work_size_small(cl::Device::getDefault());

  const cl::NDRange global_work_size(pv_size);

  // Call to kernel
  padding_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                 deviceVolume, devicePaddedVolume, paddingIndex, paddingSize,
                 inputSize, weights_d, w_buf_sizeof);

  // Execute kernel
  err = queue.finish();
  checkError(err, "Padding finish queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, devicePaddedVolume, hostPaddedVolume,
                 hostPaddedVolume + pv_size);
  checkError(err, "Padding copy back from device");
}

void OpenCL_subtract2Dfrom3DbySlice_InPlace(
    FloatImageType::Pointer &projections,
    const FloatImage2DType::Pointer &filter) {

  const auto inputSize = projections->GetBufferedRegion().GetSize();
  const auto subSize = filter->GetBufferedRegion().GetSize();

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput =
      memorySizeInput * sizeof(FloatImageType::ValueType);

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  kernel_man.initialize(defines, memoryByteSizeInput);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  auto *buffer = projections->GetBufferPointer();
  auto *sub_buffer = filter->GetBufferPointer();
  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  const auto deviceSubBuffer =
      cl::Buffer(ctx, sub_buffer, sub_buffer + (subSize[0] * subSize[1]), true,
                 true, nullptr);

  // Create program
  auto subtract_kernel2D = cl::KernelFunctor<cl::Buffer, cl_uint4, cl::Buffer>(
      kernel_man.getKernel(en_subtract_kernel2D));

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  const cl_uint4 inputDim = {{static_cast<cl_uint>(inputSize[0]),
                              static_cast<cl_uint>(inputSize[1]),
                              static_cast<cl_uint>(inputSize[2]), 0}};

  subtract_kernel2D(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                    deviceBuffer, inputDim, deviceSubBuffer);

  err = queue.finish();
  checkError(err, "Finish subtract2d queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  checkError(err, "Sub, copy back from device");
}

itk::Image<float, 3U>::Pointer OpenCL_divide3Dby3D_OutOfPlace(
    const itk::Image<unsigned short, 3U>::Pointer &Num3D,
    const itk::Image<unsigned short, 3U>::Pointer &Denum3D) {

  const auto region = Num3D->GetLargestPossibleRegion();
  const auto buffer = Num3D->GetBufferPointer();
  auto inputSize = region.GetSize();
  const auto sub_buffer = Denum3D->GetBufferPointer();
  auto subSize = Denum3D->GetLargestPossibleRegion().GetSize();

  auto outImage = itk::Image<float, 3U>::New();
  const auto projCT_size = region.GetSize();
  const auto projCT_idxStart = region.GetIndex();
  const auto projCT_spacing = Num3D->GetSpacing();
  const auto projCT_origin = Num3D->GetOrigin();

  itk::Image<float, 3U>::RegionType projCT_region;
  projCT_region.SetSize(projCT_size);
  projCT_region.SetIndex(projCT_idxStart);

  outImage->SetRegions(projCT_region);
  outImage->SetSpacing(projCT_spacing);
  outImage->SetOrigin(projCT_origin);

  outImage->Allocate();
  auto *out_buffer = outImage->GetBufferPointer();

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_ushort);

  const auto memorySizeSub = subSize[0] * subSize[1] * subSize[2];

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  kernel_man.initialize(defines, memoryByteSizeInput);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  const auto deviceSubBuffer = cl::Buffer(
      ctx, sub_buffer, sub_buffer + memorySizeSub, true, true, nullptr);

  const auto deviceOutBuffer = cl::Buffer(
      ctx, out_buffer, out_buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto divide_kernel3D_Ushort =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer>(
          kernel_man.getKernel(en_divide_kernel3D_ushort));

  const auto local_work_size =
      get_local_work_size_small(cl::Device::getDefault());
  const cl::NDRange global_work_size(memorySizeInput);

  divide_kernel3D_Ushort(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      deviceSubBuffer, deviceOutBuffer);

  err = queue.finish();
  checkError(err, "Finish divide3d queue");

  /* Fetch results of calculations. */

  err = cl::copy(queue, deviceOutBuffer, out_buffer,
                 out_buffer + memorySizeInput);

  checkError(err, "Sub3D, could not copy out_buffer back from device");

  return outImage;
}

void OpenCL_AddConst_InPlace(cl_float *buffer,
                             const itk::Image<float, 3U>::SizeType &inputSize,
                             const cl_float constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  kernel_man.initialize(defines, memoryByteSizeInput);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float>(
      kernel_man.getKernel(en_add_const_kernel));

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  add_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                   deviceBuffer, constant);

  err = queue.finish();
  checkError(err, "Finish subtract2d queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  checkError(err, "Add, copy data back from device");
}

void OpenCL_AddConst_MulConst_InPlace(
    cl_float *buffer, const itk::Image<float, 3U>::SizeType &inputSize,
    const cl_float add_constant, const cl_float mul_constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  kernel_man.initialize(defines, memoryByteSizeInput);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_mul_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float, cl_float>(
      kernel_man.getKernel(en_add_mul_const_kernel));

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  add_mul_const_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      add_constant, mul_constant);

  err = queue.finish();
  checkError(err, "Finish add_mul_const queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  checkError(err, "addMul, copy data back from device");
}

void OpenCL_AddConst_InPlace_2D(
    cl_float *buffer, const itk::Image<float, 2U>::SizeType &inputSize,
    const cl_float constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  kernel_man.initialize(defines, memoryByteSizeInput);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float>(
      kernel_man.getKernel(en_add_const_kernel));

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  add_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                   deviceBuffer, constant);

  err = queue.finish();
  checkError(err, "Finish add_const queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  checkError(err, "add2d could not copy back from device");
}

cl_float2 OpenCL_min_max_3D(cl_float *buffer,
                            const itk::Image<float, 3U>::SizeType &inputSize) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];

  return OpenCL_min_max_1D(buffer, memorySizeInput);
}

cl_float2 OpenCL_min_max_2D(cl_float *buffer,
                            const itk::Image<float, 2U>::SizeType &inputSize) {

  const auto memorySizeInput = inputSize[0] * inputSize[1];

  return OpenCL_min_max_1D(buffer, memorySizeInput);
}

cl_float2 OpenCL_min_max_recurse_unified(const size_t inputSize,
                                         cl::Buffer &deviceBuffer) {
  auto err = CL_SUCCESS;
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  const auto recurse_local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  if (inputSize < 257) {
    // read device buffer
    auto *buffer = (cl_float2 *)queue.enqueueMapBuffer(
        deviceBuffer, CL_TRUE, CL_MAP_READ, 0, inputSize * sizeof(cl_float2));
    err = queue.finish();
    checkError(err, "Finish min_max_recurese queue");

    cl_float2 out;
    out.x = std::numeric_limits<float>::max();
    out.y = std::numeric_limits<float>::min();
    for (cl_uint i = 0; i < inputSize; i++) {
      if (out.x > buffer[i].x) {
        out.x = buffer[i].x;
      }
      if (out.y < buffer[i].y) {
        out.y = buffer[i].y;
      }
    }

    err = queue.enqueueUnmapMemObject(deviceBuffer, buffer);
    checkError(err, "OpenCL_min_max_recurse_unified unmap");

    return out;
  }

  // Create program
  auto min_max_kernel =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint, cl_uint>(
          kernel_man.getKernel(en_min_max_kernel2));

  /* Prepare OpenCL memory objects and place data inside them. */

  cl_uint divider = 128;
  while (true) {
    if (inputSize % divider == 0) {
      break;
    }
    divider /= 2;
    if (divider == 2) {
      break;
    }
  }

  const auto outputDim = inputSize / divider;
  auto local_work_size = recurse_local_work_size.get()[0];

  while (local_work_size > 32) {
    const auto n_work_groups =
        (outputDim + outputDim % local_work_size) / local_work_size;
    if (n_work_groups % 2 !=
            0 // global must be evenly divisable with local work group
        || n_work_groups * local_work_size < local_work_size) {
      local_work_size /= 2;
    } else {
      break;
    }
  }
  const auto global_work_size = (outputDim + outputDim % local_work_size) /
                                local_work_size *
                                local_work_size; // seems redundant, but isn't

  const auto memoryByteSizeSub = global_work_size * sizeof(cl_float2);
  // to avoid access violation in kernel

  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR, memoryByteSizeSub,
                             nullptr, &err);
  checkError(err, "Allocate device sub buffer in unified memory");

  min_max_kernel(cl::EnqueueArgs(queue, cl::NDRange(global_work_size),
                                 cl::NDRange(local_work_size)),
                 deviceBuffer, deviceSubBuffer, divider, inputSize, err);
  checkError(err, "min_max_kernel call recurse unified");

  err = queue.finish();
  checkError(err, "Finish min_max_recurse queue");

  const auto result =
      OpenCL_min_max_recurse_unified(outputDim, deviceSubBuffer);

  err = queue.finish();
  checkError(err, "Finish unmap min_max_recurse queue");

  return result;
}

cl_float2 OpenCL_min_max_unified(cl_float *buffer, const size_t memorySizeInput,
                                 const size_t memorySizeSub,
                                 const cl::NDRange local_work_size,
                                 const size_t divider) {
  auto err = CL_SUCCESS;
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  // Create program
  auto min_max_kernel = cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint>(
      kernel_man.getKernel(en_min_max_kernel));

  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR, memoryByteSizeSub);

  const cl::NDRange global_work_size(memorySizeSub);

  min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                 deviceBuffer, deviceSubBuffer, divider);

  err = queue.finish();
  checkError(err, "Finish add_const queue");

  const auto result =
      OpenCL_min_max_recurse_unified(memorySizeSub, deviceSubBuffer);

  err = queue.finish();
  checkError(err, "add_const finish");

  return result;
}

cl_float2 OpenCL_min_max_pinned(cl_float *buffer, const size_t memorySizeInput,
                                const size_t memorySizeSub,
                                const cl::NDRange local_work_size,
                                const size_t divider) {

  auto err = CL_SUCCESS;
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  // Create program
  auto min_max_kernel = cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint>(
      kernel_man.getKernel(en_min_max_kernel));

  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  const cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR,
                                         memoryByteSizeSub);
  auto *sub_buffer = (cl_float2 *)queue.enqueueMapBuffer(
      devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  const cl::NDRange global_work_size(memorySizeSub);

  min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                 deviceBuffer, deviceSubBuffer, divider);

  err = queue.finish();
  checkError(err, "Finish add_const queue");

  const auto result =
      OpenCL_min_max_recurse(sub_buffer, memorySizeSub, deviceSubBuffer);

  err = queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  checkError(err, "add_const unmap");
  err = queue.finish();
  checkError(err, "add_const finish");

  return result;
}

cl_float2 OpenCL_min_max_recurse(cl_float2 *buffer, const size_t inputSize,
                                 cl::Buffer &deviceBuffer) {

  const auto recurse_local_work_size =
      get_local_work_size_large(cl::Device::getDefault());
  auto queue = cl::CommandQueue::getDefault();

  if (inputSize < 257) {
    auto err = queue.enqueueReadBuffer(deviceBuffer, CL_TRUE, 0,
                                       inputSize * sizeof(cl_float2), buffer);
    checkError(err, "min_max_recurse read buffer");
    err = queue.finish();
    checkError(err, "Finish min_max_recurese queue");

    cl_float2 out;
    out.x = std::numeric_limits<float>::max();
    out.y = std::numeric_limits<float>::min();
    for (cl_uint i = 0; i < inputSize; i++) {
      if (out.x > buffer[i].x) {
        out.x = buffer[i].x;
      }
      if (out.y < buffer[i].y) {
        out.y = buffer[i].y;
      }
    }
    return out;
  }

  // Create program
  auto min_max_kernel =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint, cl_uint>(
          kernel_man.getKernel(en_min_max_kernel2));

  /* Prepare OpenCL memory objects and place data inside them. */

  cl_uint divider = 128;
  while (true) {
    if (inputSize % divider == 0) {
      break;
    }
    divider /= 2;
    if (divider == 2) {
      break;
    }
  }

  const auto outputDim = inputSize / divider;
  auto local_work_size = recurse_local_work_size.get()[0];

  while (local_work_size > 32) {
    const auto n_work_groups =
        (outputDim + outputDim % local_work_size) / local_work_size;
    if (n_work_groups % 2 !=
            0 // global must be evenly divisable with local work group
        || n_work_groups * local_work_size < local_work_size) {
      local_work_size /= 2;
    } else {
      break;
    }
  }
  const auto global_work_size = (outputDim + outputDim % local_work_size) /
                                local_work_size *
                                local_work_size; // seems redundant, but isn't

  const auto memoryByteSizeSub = global_work_size * sizeof(cl_float2);
  // to avoid access violation in kernel

  const auto ctx = cl::Context::getDefault();
  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  const cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR,
                                         memoryByteSizeSub);
  auto *sub_buffer = (cl_float2 *)queue.enqueueMapBuffer(
      devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  auto err = CL_SUCCESS;
  min_max_kernel(cl::EnqueueArgs(queue, cl::NDRange(global_work_size),
                                 cl::NDRange(local_work_size)),
                 deviceBuffer, deviceSubBuffer, divider, inputSize, err);
  checkError(err, "min_max_kernel call recurse pinned");

  err = queue.finish();
  checkError(err, "Finish min_max_recurse queue");

  /* Fetch results of calculations. */
  err = queue.enqueueReadBuffer(deviceSubBuffer, CL_TRUE, 0, memoryByteSizeSub,
                                sub_buffer);
  checkError(err, "read sub buffer min_max_recurse queue");
  err = queue.finish();
  checkError(err, "Finish read sub buffer min_max_recurse queue");

  const auto result =
      OpenCL_min_max_recurse(sub_buffer, outputDim, deviceSubBuffer);

  err = queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  checkError(err, "Unmap sub pinned min_max_recurse queue");
  err = queue.finish();
  checkError(err, "Finish unmap min_max_recurse queue");

  return result;
}

cl_float2 OpenCL_min_max_1D(cl_float *buffer, const size_t memorySizeInput) {

  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);
  auto defines = std::string("");
  kernel_man.initialize(defines, memoryByteSizeInput);

  const auto local_work_size =
      get_local_work_size_small(cl::Device::getDefault());

  size_t divider = 128;
  while (true) {
    if (memorySizeInput % divider == 0) {
      break;
    }
    divider /= 2;
    if (divider == 2) {
      break;
    }
  }

  const auto memorySizeSub = memorySizeInput / divider;

#if CL_HPP_TARGET_OPENCL_VERSION < 200
  auto dev = cl::Device::getDefault();
  const auto unified_memory = dev.getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>();
  if (unified_memory) {
    return OpenCL_min_max_unified(buffer, memorySizeInput, memorySizeSub,
                                  local_work_size, divider);
  }
#endif
  return OpenCL_min_max_pinned(buffer, memorySizeInput, memorySizeSub,
                               local_work_size, divider);
}

void OpenCL_crop_by_struct_InPlace(UShortImageType::Pointer &ct_image,
                                   const Rtss_roi_modern &voi) {

  const auto inputSize = ct_image->GetBufferedRegion().GetSize();
  const cl_ulong4 in_size = {{inputSize[0], inputSize[1], inputSize[2], 0}};

  const auto inputOrigin = ct_image->GetOrigin();
  const cl_float2 in_orig = {{static_cast<cl_float>(inputOrigin[0]),
                              static_cast<cl_float>(inputOrigin[1])}};

  const auto inputSpacing = ct_image->GetSpacing();
  const cl_float2 in_spacing = {{static_cast<cl_float>(inputSpacing[0]),
                                 static_cast<cl_float>(inputSpacing[1])}};

  const auto memorySizeInput = inputSize[0] * inputSize[1];

  auto defines = std::string("");
  kernel_man.initialize(defines, memorySizeInput * sizeof(cl_ushort));

  auto err = CL_SUCCESS;
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /*
   * __kernel void crop_by_struct_kernel(__global ushort *dev_vol,
   *                                 __constant float2 *structure,
   *                                 ulong number_of_verti, ulong4 vol_dim,
   *                                 float2 vol_offset, float2 vol_spacing)
   */
  // Create program
  auto crop_kernel = cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_ulong,
                                       cl_ulong4, cl_float2, cl_float2>(
      kernel_man.getKernel(en_crop_by_struct_kernel));

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());
  const cl::NDRange global_work_size(memorySizeInput);

  for (auto slice = 0ul; slice < inputSize[2]; ++slice) {
    auto *buffer = ct_image->GetBufferPointer() + slice * memorySizeInput;

    /* Prepare OpenCL memory objects and place data inside them. */
    const auto deviceBuffer =
        cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, &err);
    checkError(err, "Create device buffer ushort image");

    std::vector<cl_float2> struct_buffer;
    auto prev_z_pos = voi.pslist.at(1).coordinates.at(0).z;
    for (auto &contour : voi.pslist) {
      const auto slice_pos = inputOrigin[2] + inputSpacing[2] * slice;
      const auto z_pos = contour.coordinates.at(0).z;
      if (fabs(z_pos - slice_pos) < fabs(z_pos - prev_z_pos)) {
        for (const auto coord : contour.coordinates) {
          struct_buffer.push_back({coord.x, coord.y});
        }
        break;
      }
      prev_z_pos = z_pos;
    }

    cl::Buffer d_struct_buffer;
    if (!struct_buffer.empty()) {
      d_struct_buffer = cl::Buffer(struct_buffer.begin(), struct_buffer.end(),
                                   true, false, &err);
      checkError(err, "Create structure buffer");
    }

    const cl_ulong n_coords = struct_buffer.size();

    crop_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                deviceBuffer, d_struct_buffer, n_coords, in_size, in_orig,
                in_spacing);

    err = queue.finish();
    checkError(err, "Finish crop by struct queue");

    /* Fetch results of calculations. */
    err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

    checkError(err, "Add, copy data back from device");
  }
}
