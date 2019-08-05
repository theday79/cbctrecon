// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "ImageFilters.h"
#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"
#include "OpenCL/util.hpp"

// If this function fails, something is probably wrong with the local OpenCL ICD
// setup. Therefore we will just let it run and not throw, so the user gets as
// much info as possible before the crash
const auto getDeviceByReqAllocSize(const size_t required_mem_alloc_size) {

  auto devices = std::vector<cl::Device>();
  getDeviceList(devices);

  // Maybe read device number from file?

  for (auto & dev : devices) {
    cl_ulong max_mem_alloc;
    dev.getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &max_mem_alloc);
    if (max_mem_alloc > required_mem_alloc_size){
      return dev;
    }
  }
  return devices.at(0);
}

cl::NDRange get_local_work_size_small(cl::Device device) {
  auto device_name = std::string();
  device.getInfo(CL_DEVICE_NAME, &device_name);
  if (device_name == "Intel(R) Iris(TM) Pro Graphics 5200") {
    return cl::NDRange(16);
  }
  return cl::NDRange(128);
}
cl::NDRange get_local_work_size_large(cl::Device device) {
  auto device_name = std::string();
  device.getInfo(CL_DEVICE_NAME, &device_name);
  if (device_name == "Intel(R) Iris(TM) Pro Graphics 5200") {
    return cl::NDRange(32);
  }
  return cl::NDRange(128);
}

void OpenCL_padding(const cl_int4 &paddingIndex, const cl_uint4 &paddingSize,
                    const cl_uint4 &inputSize,
                    const float *hostVolume, // input
                    float *hostPaddedVolume, // output
                    std::vector<float> &mirrorWeights) {
  auto err = CL_SUCCESS;

  const auto pv_size =
      paddingSize.x * paddingSize.y * paddingSize.z;
  const auto pv_buffer_size = pv_size * sizeof(float);

  const auto dev = getDeviceByReqAllocSize(pv_buffer_size);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);


  /* Prepare OpenCL memory objects and place data inside them. */

  // just because intellisense couldn't understand it below...
  const auto w_buf_sizeof = static_cast<cl_uint>(mirrorWeights.size());

  auto weights_d = cl::Buffer(ctx, mirrorWeights.begin(), mirrorWeights.end(), /*read only*/true, /*use host ptr*/ false, &err);
  if (err != CL_SUCCESS) {
    std::cout
        << "PAD::Could not create weigths device buffer, error code: "
        << err << std::endl;
  }

  auto devicePaddedVolume = cl::Buffer(ctx, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, pv_buffer_size, (void*)&hostPaddedVolume[0], &err);

  if (err != CL_SUCCESS) {
    std::cout
        << "PAD::Could not write OpenCL devicePaddedVolume buffer, error code: "
        << err << std::endl;
  }

  const auto v_buffer_size = sizeof(hostVolume);
  auto deviceVolume = cl::Buffer(ctx, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, v_buffer_size, (void *)&hostVolume[0], &err);

  if (err != CL_SUCCESS) {
    std::cout << "PAD::Could not write OpenCL deviceVolume buffer, error code: "
              << err << std::endl;
  }


  cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_int4, cl_uint4, cl_uint4, cl::Buffer, cl_uint> padding_kernel(program, "padding_kernel");

  auto local_work_size = get_local_work_size_small(dev);

  const cl::NDRange global_work_size(pv_size);

  // Call to kernel
  padding_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
   deviceVolume, devicePaddedVolume, paddingIndex, paddingSize, inputSize, weights_d, w_buf_sizeof);

  // Execute kernel
  queue.finish();

  /* Fetch results of calculations. */
  err = cl::copy(queue, devicePaddedVolume, hostPaddedVolume, hostPaddedVolume + pv_size);
  if (err != CL_SUCCESS){
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
  const auto memoryByteSizeInput = memorySizeInput *
                               sizeof(FloatImageType::ValueType);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  const auto memorySizeSub =
      subSize[0] * subSize[1] * sizeof(FloatImage2DType::ValueType);

  auto *buffer = projections->GetBufferPointer();
  auto *sub_buffer = filter->GetBufferPointer();
  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  auto deviceSubBuffer = cl::Buffer(ctx, sub_buffer, sub_buffer + (subSize[0] * subSize[1]), true, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl_uint4, cl::Buffer> subtract_kernel2D(program, "subtract_kernel2D");

  auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  const cl_uint4 inputDim = {{static_cast<cl_uint>(inputSize[0]),
                              static_cast<cl_uint>(inputSize[1]),
                              static_cast<cl_uint>(inputSize[2]), 0}};

  subtract_kernel2D(cl::EnqueueArgs(queue, global_work_size, local_work_size),
          deviceBuffer, inputDim, deviceSubBuffer);

  queue.finish();

  /* Fetch results of calculations. */
  auto err = CL_SUCCESS;
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
  auto out_buffer = outImage->GetBufferPointer();

  const auto memorySizeInput =
      inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_ushort);
  const auto memoryByteSizeOutput = memorySizeInput * sizeof(cl_float);

  const auto memorySizeSub =
      subSize[0] * subSize[1] * subSize[2];
  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_ushort);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  auto deviceSubBuffer = cl::Buffer(ctx, sub_buffer, sub_buffer + memorySizeSub, true, true, nullptr);

  auto deviceOutBuffer = cl::Buffer(ctx, out_buffer, out_buffer + memorySizeInput, false, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer> divide_kernel3D_ushort(program, "divide_kernel3D_Ushort");

  auto local_work_size = get_local_work_size_small(dev);
  const cl::NDRange global_work_size(memorySizeInput);

  divide_kernel3D_ushort(cl::EnqueueArgs(queue, global_work_size, local_work_size),
          deviceBuffer, deviceSubBuffer, deviceOutBuffer);

  queue.finish();

  /* Fetch results of calculations. */

  auto err = CL_SUCCESS;
  err = cl::copy(queue, deviceOutBuffer, out_buffer, out_buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "SUB3D::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }

  return outImage;
}

void OpenCL_AddConst_InPlace(cl_float *buffer,
                             const itk::Image<float, 3U>::SizeType &inputSize,
                             const cl_float constant) {

  const auto memorySizeInput =
      inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput =
      memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl_float> add_const_kernel(program, "add_const_kernel");

  auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  add_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
          deviceBuffer, constant);

  queue.finish();

  /* Fetch results of calculations. */
  auto err = CL_SUCCESS;
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "ADD::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }

}

void OpenCL_AddConst_MulConst_InPlace(
    cl_float *buffer, const itk::Image<float, 3U>::SizeType &inputSize,
    const cl_float add_constant, const cl_float mul_constant) {

  const auto memorySizeInput =
      inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput =
      memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl_float, cl_float> add_mul_const_kernel(program, "add_mul_const_kernel");

  auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  add_mul_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
          deviceBuffer, add_constant, mul_constant);

  queue.finish();

  /* Fetch results of calculations. */
  auto err = CL_SUCCESS;
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "ADDMUL::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }

}

void OpenCL_AddConst_InPlace_2D(
    cl_float *buffer, const itk::Image<float, 2U>::SizeType &inputSize,
    const cl_float constant) {

  const auto memorySizeInput =
      inputSize[0] * inputSize[1];
  const auto memoryByteSizeInput =
      memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl_float> add_const_kernel(program, "add_const_kernel");

  auto local_work_size = get_local_work_size_large(dev);

  const cl::NDRange global_work_size(memorySizeInput);

  add_const_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer, constant);

  queue.finish();

  /* Fetch results of calculations. */
  auto err = CL_SUCCESS;
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  if (err != CL_SUCCESS) {
    std::cout << "ADD::Could not read OpenCL buffer, error code: " << err
              << std::endl;
  }

}

cl_float2 OpenCL_min_max(cl_float *buffer,
                         const itk::Image<float, 3U>::SizeType &inputSize) {

  const auto memorySizeInput =
      inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput =
      memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint> min_max_kernel(program, "min_max_kernel");

  auto local_work_size = get_local_work_size_small(dev);

  cl_uint divider = 128;
  while (true) {
    if (inputSize[0] % divider == 0) {
      break;
    }
    divider /= 2;
    if (divider == 2) {
      break;
    }
  }

  const cl_uint4 outputDim = {{static_cast<cl_uint>(inputSize[0]) / divider,
                               static_cast<cl_uint>(inputSize[1]),
                               static_cast<cl_uint>(inputSize[2]), 0}};

  const auto memorySizeSub = outputDim.x * outputDim.y * outputDim.z;
  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  const cl::NDRange global_work_size(memorySizeSub);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR, memoryByteSizeSub);
  cl_float2 *sub_buffer = (cl_float2*)queue.enqueueMapBuffer(
          devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
          deviceBuffer, deviceSubBuffer, divider);

  queue.finish();

  queue.enqueueReadBuffer(deviceSubBuffer, CL_TRUE, 0, memoryByteSizeSub, sub_buffer);
  queue.finish();

  const auto recurse_local_work_size = get_local_work_size_large(dev);
  const auto result = OpenCL_min_max_recurse(sub_buffer, memorySizeSub,
          deviceSubBuffer, ctx, queue, program,
          recurse_local_work_size);

  queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  queue.finish();

  return result;
}

cl_float2 OpenCL_min_max_2D(cl_float *buffer,
                            const itk::Image<float, 2U>::SizeType &inputSize) {

  const auto memorySizeInput = inputSize[0] * inputSize[1];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  const auto dev = getDeviceByReqAllocSize(memoryByteSizeInput);
  cl::Context ctx(dev);
  cl::CommandQueue queue(ctx);
  const auto cl_source = util::loadProgram("fdk_opencl.cl");
  auto program = cl::Program(ctx, cl_source, /* build? */ true);

  /* Prepare OpenCL memory objects and place data inside them. */
  auto deviceBuffer = cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  // Create program
  cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint> min_max_kernel(program, "min_max_kernel");

  auto local_work_size = get_local_work_size_small(dev);

  cl_uint divider = 128;
  while (true) {
    if (inputSize[0] % divider == 0) {
      break;
    }
    divider /= 2;
    if (divider == 2) {
      break;
    }
  }

  const cl_uint2 outputDim = {{static_cast<cl_uint>(inputSize[0]) / divider,
                               static_cast<cl_uint>(inputSize[1])}};

  const auto memorySizeSub = outputDim.x * outputDim.y;
  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  const cl::NDRange global_work_size(memorySizeSub);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR, memoryByteSizeSub);
  cl_float2 *sub_buffer = (cl_float2*)queue.enqueueMapBuffer(
          devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
          deviceBuffer, deviceSubBuffer, divider);

  queue.finish();

  queue.enqueueReadBuffer(deviceSubBuffer, CL_TRUE, 0, memoryByteSizeSub, sub_buffer);
  queue.finish();

  auto recurse_local_work_size = get_local_work_size_large(dev);
  const auto result = OpenCL_min_max_recurse(sub_buffer, memorySizeSub,
                                deviceSubBuffer, ctx, queue, program,
                                recurse_local_work_size);

  queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  queue.finish();

  return result;
}

cl_float2 OpenCL_min_max_recurse(const cl_float2 *buffer,
                                 const cl_uint inputSize,
                                 cl::Buffer &deviceBuffer,
                                 cl::Context &ctx,
                                 cl::CommandQueue &queue,
                                 cl::Program &program,
                                 cl::NDRange nd_local_work_size) {

  if (inputSize < 257) {
    cl_float2 out = {{65535.0, -9999.0}};
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
  cl::KernelFunctor<cl::Buffer, cl::Buffer, cl_uint, cl_uint> min_max_kernel(program, "min_max_kernel2");

  const auto memoryByteSizeInput = inputSize * sizeof(cl_float2);

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
  cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR, memoryByteSizeSub);
  cl_float2 *sub_buffer = (cl_float2*)queue.enqueueMapBuffer(
          devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  min_max_kernel(cl::EnqueueArgs(queue, cl::NDRange(global_work_size), cl::NDRange(local_work_size)),
          deviceBuffer, deviceSubBuffer, divider, outputDim);

  queue.finish();

  /* Fetch results of calculations. */
  queue.enqueueReadBuffer(deviceSubBuffer, CL_TRUE, 0, memoryByteSizeSub,
                            sub_buffer);
  queue.finish();

  const auto result = OpenCL_min_max_recurse(sub_buffer, outputDim, deviceSubBuffer,
          ctx, queue, program,
          nd_local_work_size);

  queue.enqueueUnmapMemObject(devicePinnedSubBuffer, sub_buffer);
  queue.finish();

  return result;
}

