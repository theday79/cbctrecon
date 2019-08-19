// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <cassert>

#include "ImageFilters.h"
#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"
#include "OpenCL/util.hpp"

// On Linux with intel opencl runtime, you can debug opencl kernels:
// #define DEBUG_OPENCL

enum enKernel {
  en_fdk_kernel_nn,
  en_OpenCLFDKBackProjectionImageFilterKernel,
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
  en_min_max_kernel2
};

class KernelMan {

private:
  cl::Program m_program;
  bool m_initialized = false;
  std::vector<cl::Kernel> m_kernel_list;

  static cl::Program build_ocl_program(cl::Context &ctx) {
#if CL_HPP_TARGET_OPENCL_VERSION >= 210 ||                                     \
    (CL_HPP_TARGET_OPENCL_VERSION == 200 && defined(CL_HPP_USE_IL_KHR))
    // const auto cl_source = util::loadProgramIL("kernels_cl20.spv");
    const auto cl_source = util::loadProgram("fdk_opencl.cl");
#else
    const auto cl_source = util::loadProgram("fdk_opencl.cl");
#endif
    auto err = CL_SUCCESS;
#ifdef DEBUG_OPENCL
    auto program = cl::Program(ctx, cl_source, /* build? */ false, &err);
    // For use with GDB, Requires intel OpenCL runtime
    err = program.build(" -g -s "
                        "/home/andreas/Projects/build-cbct/clang8-itk5-Debug/"
                        "bin/fdk_opencl.cl");
#else
    auto program = cl::Program(ctx, cl_source, /* build? */ true, &err);
#endif

    if (err != CL_SUCCESS) {
      std::cerr << "Could not build OpenCL program! error code: " << err
                << "\n";
    }
    return program;
  }

public:
  cl::Context m_ctx;
  cl::Device m_dev;
  KernelMan() = default;
  void initialize(const cl::Device &dev) {
    if (m_initialized && dev == m_dev) {
      return;
    }
    m_ctx = cl::Context(dev);
    m_program = build_ocl_program(m_ctx);
    m_program.createKernels(&m_kernel_list);
    m_initialized = true;
  }

  cl::Kernel getKernel(const enKernel kernel) {
    assert(m_initialized);

    auto &out_kernel = m_kernel_list.at(kernel);
    const auto kernel_name = out_kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
    std::cerr << kernel_name << "\n";
    return out_kernel;
  }
};

static KernelMan kernel_man;

// If this function fails, something is probably wrong with the local OpenCL ICD
// setup. Therefore we will just let it run and not throw, so the user gets as
// much info as possible before the crash
auto getDeviceByReqAllocSize(const size_t required_mem_alloc_size) {

  auto devices = std::vector<cl::Device>();
  OpenCL_getDeviceList(devices);

  auto default_id = 0;
  // Maybe read device number from file?

  for (auto &dev : devices) {
    const auto dev_available = dev.getInfo<CL_DEVICE_AVAILABLE>();
    if (!dev_available) {
      ++default_id;
      continue;
    }
    const auto max_mem_alloc = dev.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
    if (max_mem_alloc > required_mem_alloc_size) {
      return dev;
    }
  }
  // Let's hope that max mem alloc is not strictly necessary
  return devices.at(default_id);
}

void OpenCL_initialize(const size_t required_mem_alloc_size) {
  const auto dev = getDeviceByReqAllocSize(required_mem_alloc_size);
  kernel_man.initialize(dev);
}

void OpenCL_initialize(cl::Device &dev) { kernel_man.initialize(dev); }

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

  const auto pv_size = paddingSize.x * paddingSize.y * paddingSize.z;
  const auto pv_buffer_size = pv_size * sizeof(float);

  const auto dev = getDeviceByReqAllocSize(pv_buffer_size);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

  /* Prepare OpenCL memory objects and place data inside them. */

  // just because intellisense couldn't understand it below...
  const auto w_buf_sizeof = static_cast<cl_uint>(mirrorWeights.size());

  const auto weights_d =
      cl::Buffer(ctx, mirrorWeights.begin(), mirrorWeights.end(),
                 /*read only*/ true, /*use host ptr*/ false, &err);
  if (err != CL_SUCCESS) {
    std::cout << "PAD::Could not create weigths device buffer, error code: "
              << err << std::endl;
  }

  const auto devicePaddedVolume =
      cl::Buffer(ctx, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, pv_buffer_size,
                 (void *)&hostPaddedVolume[0], &err);

  if (err != CL_SUCCESS) {
    std::cout
        << "PAD::Could not write OpenCL devicePaddedVolume buffer, error code: "
        << err << std::endl;
  }

  const auto v_buffer_size = sizeof(hostVolume);
  const auto deviceVolume =
      cl::Buffer(ctx, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, v_buffer_size,
                 (void *)&hostVolume[0], &err);

  if (err != CL_SUCCESS) {
    std::cout << "PAD::Could not write OpenCL deviceVolume buffer, error code: "
              << err << std::endl;
  }

  auto padding_kernel =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_int4, cl_uint4, cl_uint4,
                        cl::Buffer, cl_uint>(
          kernel_man.getKernel(en_padding_kernel));
  // padding_kernel(program, "padding_kernel");

  const auto local_work_size = get_local_work_size_small(dev);

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
  if (err != CL_SUCCESS) {
    std::cerr << "PAD::Something went wrong, error code: " << err << "\n";
  }

  // std::cout << "Did it work? we just don't know..." << std::endl;
}

void OpenCL_subtract2Dfrom3DbySlice_InPlace(
    FloatImageType::Pointer &projections,
    const FloatImage2DType::Pointer &filter) {

  const auto inputSize = projections->GetBufferedRegion().GetSize();
  const auto subSize = filter->GetBufferedRegion().GetSize();

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput =
      memorySizeInput * sizeof(FloatImageType::ValueType);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

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

  const auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  const cl_uint4 inputDim = {{static_cast<cl_uint>(inputSize[0]),
                              static_cast<cl_uint>(inputSize[1]),
                              static_cast<cl_uint>(inputSize[2]), 0}};

  subtract_kernel2D(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                    deviceBuffer, inputDim, deviceSubBuffer);

  auto err = queue.finish();
  checkError(err, "Finish subtract2d queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "SUB::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }
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

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

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

  const auto local_work_size = get_local_work_size_small(dev);
  const cl::NDRange global_work_size(memorySizeInput);

  divide_kernel3D_Ushort(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      deviceSubBuffer, deviceOutBuffer);

  auto err = queue.finish();
  checkError(err, "Finish divide3d queue");

  /* Fetch results of calculations. */

  err = cl::copy(queue, deviceOutBuffer, out_buffer,
                 out_buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "SUB3D::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }

  return outImage;
}

void OpenCL_AddConst_InPlace(cl_float *buffer,
                             const itk::Image<float, 3U>::SizeType &inputSize,
                             const cl_float constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float>(
      kernel_man.getKernel(en_add_const_kernel));

  const auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  add_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                   deviceBuffer, constant);

  auto err = queue.finish();
  checkError(err, "Finish subtract2d queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "ADD::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }
}

void OpenCL_AddConst_MulConst_InPlace(
    cl_float *buffer, const itk::Image<float, 3U>::SizeType &inputSize,
    const cl_float add_constant, const cl_float mul_constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_mul_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float, cl_float>(
      kernel_man.getKernel(en_add_mul_const_kernel));

  const auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  add_mul_const_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      add_constant, mul_constant);

  auto err = queue.finish();
  checkError(err, "Finish add_mul_const queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "ADDMUL::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }
}

void OpenCL_AddConst_InPlace_2D(
    cl_float *buffer, const itk::Image<float, 2U>::SizeType &inputSize,
    const cl_float constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float>(
      kernel_man.getKernel(en_add_const_kernel));

  const auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  add_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                   deviceBuffer, constant);

  auto err = queue.finish();
  checkError(err, "Finish add_const queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "ADD::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }
}

cl_float2 OpenCL_min_max_3D(cl_float *buffer,
                            const itk::Image<float, 3U>::SizeType &inputSize) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];

  return OpenCL_min_max_1D(buffer, memorySizeInput);
  ;
}

cl_float2 OpenCL_min_max_2D(cl_float *buffer,
                            const itk::Image<float, 2U>::SizeType &inputSize) {

  const auto memorySizeInput = inputSize[0] * inputSize[1];

  return OpenCL_min_max_1D(buffer, memorySizeInput);
}

cl_float2 OpenCL_min_max_1D(cl_float *buffer, const size_t memorySizeInput) {
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  kernel_man.initialize(dev);
  auto &ctx = kernel_man.m_ctx;
  cl::CommandQueue queue(ctx);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  // Create program
  auto min_max_kernel = cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint>(
      kernel_man.getKernel(en_min_max_kernel));

  const auto local_work_size = get_local_work_size_small(dev);

  cl_uint divider = 128;
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
  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  const cl::NDRange global_work_size(memorySizeSub);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  const cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR,
                                         memoryByteSizeSub);
  auto *sub_buffer = (cl_float2 *)(queue.enqueueMapBuffer(
      devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub));

  min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                 deviceBuffer, deviceSubBuffer, divider);

  auto err = queue.finish();
  checkError(err, "Finish add_const queue");

  const auto recurse_local_work_size = get_local_work_size_large(dev);
  const auto result =
      OpenCL_min_max_recurse(sub_buffer, memorySizeSub, deviceSubBuffer, queue,
                             recurse_local_work_size);

  err = queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  checkError(err, "add_const unmap");
  err = queue.finish();
  checkError(err, "add_const finish");

  return result;
}

cl_float2 OpenCL_min_max_recurse(cl_float2 *buffer, const cl_uint inputSize,
                                 cl::Buffer &deviceBuffer,
                                 cl::CommandQueue &queue,
                                 cl::NDRange nd_local_work_size) {

  auto &ctx = kernel_man.m_ctx;

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
  auto local_work_size = nd_local_work_size.size();

  while (local_work_size > 32) {
    if ((outputDim + outputDim % local_work_size) / local_work_size % 2 !=
        0) { // global must be evenly divisable with local work group
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

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  const cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR,
                                         memoryByteSizeSub);
  auto *sub_buffer = (cl_float2 *)queue.enqueueMapBuffer(
      devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  min_max_kernel(cl::EnqueueArgs(queue, cl::NDRange(global_work_size),
                                 cl::NDRange(local_work_size)),
                 deviceBuffer, deviceSubBuffer, divider, inputSize);

  auto err = queue.finish();
  checkError(err, "Finish min_max_recurse queue");

  /* Fetch results of calculations. */
  err = queue.enqueueReadBuffer(deviceSubBuffer, CL_TRUE, 0, memoryByteSizeSub,
                                sub_buffer);
  checkError(err, "read sub buffer min_max_recurse queue");
  err = queue.finish();
  checkError(err, "Finish read sub buffer min_max_recurse queue");

  const auto result = OpenCL_min_max_recurse(
      sub_buffer, outputDim, deviceSubBuffer, queue, nd_local_work_size);

  err = queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  checkError(err, "Unmap sub pinned min_max_recurse queue");
  err = queue.finish();
  checkError(err, "Finish unmap min_max_recurse queue");

  return result;
}
