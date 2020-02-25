// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <cassert>
#include <unordered_map>

#include "OpenCL/ImageFilters.hpp"

#include "OpenCL/cl2.hpp"
#include "OpenCL/device_picker.hpp"
#include "OpenCL/util.hpp"

// For crop_by_structure
#include "StructureSet.h"

// #define USE_XEON_PHI

// On Linux with intel opencl runtime, you can debug opencl kernels:
// #define DEBUG_OPENCL

#ifdef LOWPASS_FFT
// ITK Low-pass fourier filter
#include <itkFFTShiftImageFilter.h>
#include <itkForwardFFTImageFilter.h>
#include <itkGaussianImageSource.h>
#include <itkInverseFFTImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkWrapPadImageFilter.h>
#else
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#endif // LOWPASS_FFT

// We want to use the opencl types as template args, but don't care that the
// attributes (alignment) is ignored
#ifndef _WIN32
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

enum class enKernel {
  en_padding_kernel,
  en_multiply_kernel,
  en_multiply_kernel2D,
  en_subtract_kernel2D,
  en_divide_kernel3D_loginv_ushort,
  en_add_const_kernel,
  en_add_const_with_thresh_kernel,
  en_add_mul_const_kernel,
  en_add_mul_const_with_thresh_kernel,
  en_min_max_kernel,
  en_min_max_kernel2,
  en_crop_by_struct_kernel,
  en_LogItoI_subtract_median_ItoLogI,
  en_LogItoI_subtract_median_y_x,
  en_ItoLogI
};

template <typename E, typename T, size_t N>
using enum_map = std::array<std::pair<E, T>, N>;

using kernel_map = enum_map<enKernel, const char *, 15>;
constexpr kernel_map kernel_names = {
    {{enKernel::en_padding_kernel, "padding_kernel"},
     {enKernel::en_multiply_kernel, "multiply_kernel"},
     {enKernel::en_multiply_kernel2D, "multiply_kernel2D"},
     {enKernel::en_subtract_kernel2D, "subtract_kernel2D"},
     {enKernel::en_divide_kernel3D_loginv_ushort,
      "divide_kernel3D_loginv_Ushort"},
     {enKernel::en_add_const_kernel, "add_const_kernel"},
     {enKernel::en_add_const_with_thresh_kernel,
      "add_const_with_thresh_kernel"},
     {enKernel::en_add_mul_const_kernel, "add_mul_const_kernel"},
     {enKernel::en_add_mul_const_with_thresh_kernel,
      "add_mul_const_with_thresh_kernel"},
     {enKernel::en_min_max_kernel, "min_max_kernel"},
     {enKernel::en_min_max_kernel2, "min_max_kernel2"},
     {enKernel::en_crop_by_struct_kernel, "crop_by_struct_kernel"},
     {enKernel::en_LogItoI_subtract_median_ItoLogI,
      "log_i_to_i_subtract_median_i_to_log_i"},
     {enKernel::en_LogItoI_subtract_median_y_x,
      "log_i_to_i_subtract_median_y_x"},
     {enKernel::en_ItoLogI, "i_to_log_i_kernel"}}};

template <auto key, typename EM>
constexpr auto ce_find_key(const EM &input_arr) {
  constexpr auto N = std::tuple_size_v<EM>;
  static_assert(N > 0, "Array has 0 elements!");

  // Pair type (a special case of tuple):
  using PT = typename std::tuple_element_t<0, EM>;
  static_assert(std::tuple_size_v<PT> == 2, "Array doesn't contain pair type");

  using E = typename std::tuple_element_t<0, PT>;
  static_assert(std::is_enum_v<E>, "Key must be enum type");
  static_assert(std::is_same_v<decltype(key), E>,
                "Key didn't match enum type of array");

  for (auto key_val : input_arr) {
    if (key_val.first == key) {
      return key_val.second;
    }
  }
}

void print_prof_info(cl::Event &evt) {
  auto status = evt.getInfo<CL_EVENT_COMMAND_EXECUTION_STATUS>();
  if (status != CL_COMPLETE) {
    evt.wait();
  }
#ifdef DEBUG_OPENCL
  std::cerr << "Kernel took: "
            << static_cast<double>(
                   evt.getProfilingInfo<CL_PROFILING_COMMAND_END>() -
                   evt.getProfilingInfo<CL_PROFILING_COMMAND_START>()) *
                   1e6 // ns to ms
            << " ms\n";
#endif
}

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

    std::vector<cl::Device> devices = OpenCL_getDeviceList();

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

      if (device_type == CL_DEVICE_TYPE_GPU && // Xeon Phi is much too slow
                                               // compared with a GPU
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
    if (device != cl::Device::setDefault(device)) {
      std::cerr << "Another default device was already set!\n";
    }
    if (default_ctx != cl::Context::setDefault(default_ctx)) {
      std::cerr << "Another default context was already set!\n";
    }
    auto deviceQueue = cl::CommandQueue::setDefault(
        cl::CommandQueue(default_ctx, device, 0, err));
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

  template <enKernel kernel> cl::Kernel getKernel() const {
    assert(m_initialized);

    constexpr auto kernel_name = ce_find_key<kernel>(kernel_names);

    /*kernel_name i is captured for free in clang but not in gcc, so we'll use
     * auto-capturing, [=] */
    const auto it_kernel = std::find_if(
        m_kernel_list.begin(), m_kernel_list.end(), [=](cl::Kernel cur_kernel) {
          return (cur_kernel.getInfo<CL_KERNEL_FUNCTION_NAME>()) == kernel_name;
        });
    if (it_kernel == m_kernel_list.end()) {
      std::cerr << "OH NO, SOMETHING BAD WILL HAPPEN SOON!\n";
    }
    const auto out_kernel = *it_kernel;

#ifdef DEBUG_OPENCL
    const auto &kernel_name = out_kernel.getInfo<CL_KERNEL_FUNCTION_NAME>();
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
  } else if (device_name == "Intel(R) Many Integrated Core Acceleration Card") {
    return {16};
  }
  return {128};
}
cl::NDRange get_local_work_size_large(const cl::Device &device) {
  const auto device_name = device.getInfo<CL_DEVICE_NAME>();
  if (device_name == "Intel(R) Iris(TM) Pro Graphics 5200") {
    return {32};
  } else if (device_name == "Intel(R) Many Integrated Core Acceleration Card") {
    return {128};
  }
  return {128};
}

auto initialize_opencl(std::string &defines, size_t buffer_size) {
  auto err = CL_SUCCESS;
  kernel_man.initialize(defines, buffer_size);
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");
  return std::make_pair(ctx, queue);
}

void OpenCL_padding(const cl_int4 &paddingIndex, const cl_uint4 &paddingSize,
                    const cl_uint4 &inputSize,
                    const float *hostVolume, // input
                    float *hostPaddedVolume, // output
                    std::vector<float> &mirrorWeights) {
  auto err = CL_SUCCESS;

  const auto pv_size =
      static_cast<size_t>(paddingSize.x) * paddingSize.y * paddingSize.z;
  const auto pv_buffer_size = pv_size * sizeof(float);

  auto defines = std::string("");
  auto [ctx, queue] = initialize_opencl(defines, pv_buffer_size);

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
          kernel_man.getKernel<enKernel::en_padding_kernel>());
  // padding_kernel(program, "padding_kernel");

  const auto local_work_size =
      get_local_work_size_small(cl::Device::getDefault());

  const cl::NDRange global_work_size(pv_size);

  // Call to kernel
  auto evt =
      padding_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                     deviceVolume, devicePaddedVolume, paddingIndex,
                     paddingSize, inputSize, weights_d, w_buf_sizeof);

  // Execute kernel
  err = queue.finish();
  print_prof_info(evt);
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
  auto [ctx, queue] = initialize_opencl(defines, memoryByteSizeInput);

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
      kernel_man.getKernel<enKernel::en_subtract_kernel2D>());

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  const cl_uint4 inputDim = {{static_cast<cl_uint>(inputSize[0]),
                              static_cast<cl_uint>(inputSize[1]),
                              static_cast<cl_uint>(inputSize[2]), 0}};

  auto evt = subtract_kernel2D(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      inputDim, deviceSubBuffer);

  err = queue.finish();
  print_prof_info(evt);
  checkError(err, "Finish subtract2d queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  checkError(err, "Sub, copy back from device");
}

itk::Image<float, 3U>::Pointer OpenCL_divide3Dby3D_loginv_OutOfPlace(
    const itk::Image<unsigned short, 3U>::Pointer &Num3D,
    const itk::Image<unsigned short, 3U>::Pointer &Denum3D) {

  const auto region = Num3D->GetBufferedRegion();
  const auto buffer = Num3D->GetBufferPointer();
  auto inputSize = region.GetSize();

  const auto sub_buffer = Denum3D->GetBufferPointer();
  auto subSize = Denum3D->GetLargestPossibleRegion().GetSize();

  auto outImage = itk::Image<float, 3U>::New();
  outImage->SetRegions(region);
  outImage->Allocate();

  auto out_buffer = outImage->GetBufferPointer();

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_ushort);

  const auto memorySizeSub = subSize[0] * subSize[1] * subSize[2];

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  auto [ctx, queue] = initialize_opencl(defines, memoryByteSizeInput);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, true, true, nullptr);

  const auto deviceSubBuffer = cl::Buffer(
      ctx, sub_buffer, sub_buffer + memorySizeSub, true, true, nullptr);

  const auto deviceOutBuffer = cl::Buffer(
      ctx, out_buffer, out_buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto divide_kernel3D_loginv_Ushort =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer>(
          kernel_man.getKernel<enKernel::en_divide_kernel3D_loginv_ushort>());

  const auto local_work_size =
      get_local_work_size_small(cl::Device::getDefault());
  const cl::NDRange global_work_size(memorySizeInput);

  auto evt = divide_kernel3D_loginv_Ushort(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      deviceSubBuffer, deviceOutBuffer);

  err = queue.finish();
  print_prof_info(evt);
  checkError(err, "Finish divide3d queue");

  /* Fetch results of calculations. */

  err = cl::copy(queue, deviceOutBuffer, out_buffer,
                 out_buffer + memorySizeInput);

  checkError(err, "Sub3D, could not copy out_buffer back from device");

  return outImage;
}

#ifdef USE_XEON_PHI
#include "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.8.275\windows\mkl\include\mkl.h"

void OpenCL_AddConst_InPlace(cl_float *buffer,
                             const itk::Image<float, 3U>::SizeType &inputSize,
                             const cl_float constant) {}

#else

void OpenCL_AddConst_InPlace(cl_float *buffer,
                             const itk::Image<float, 3U>::SizeType &inputSize,
                             const cl_float constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  auto [ctx, queue] = initialize_opencl(defines, memoryByteSizeInput);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float>(
      kernel_man.getKernel<enKernel::en_add_const_kernel>());

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  auto evt = add_const_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      constant);

  err = queue.finish();
  print_prof_info(evt);
  checkError(err, "Finish subtract2d queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer, buffer, buffer + memorySizeInput);

  checkError(err, "Add, copy data back from device");
}
#endif

void OpenCL_AddConst_MulConst_InPlace(
    cl_float *buffer, const itk::Image<float, 3U>::SizeType &inputSize,
    const cl_float add_constant, const cl_float mul_constant) {

  const auto memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2];
  const auto memoryByteSizeInput = memorySizeInput * sizeof(cl_float);

  auto err = CL_SUCCESS;
  auto defines = std::string("");
  auto [ctx, queue] = initialize_opencl(defines, memoryByteSizeInput);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_mul_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float, cl_float>(
      kernel_man.getKernel<enKernel::en_add_mul_const_kernel>());

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  auto evt = add_mul_const_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      add_constant, mul_constant);

  err = queue.finish();
  print_prof_info(evt);
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
  auto [ctx, queue] = initialize_opencl(defines, memoryByteSizeInput);

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer =
      cl::Buffer(ctx, buffer, buffer + memorySizeInput, false, true, nullptr);

  // Create program
  auto add_const_kernel = cl::KernelFunctor<cl::Buffer, cl_float>(
      kernel_man.getKernel<enKernel::en_add_const_kernel>());

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());

  const cl::NDRange global_work_size(memorySizeInput);

  auto evt = add_const_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size), deviceBuffer,
      constant);

  err = queue.finish();
  print_prof_info(evt);
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
          kernel_man.getKernel<enKernel::en_min_max_kernel2>());

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

  auto evt =
      min_max_kernel(cl::EnqueueArgs(queue, cl::NDRange(global_work_size),
                                     cl::NDRange(local_work_size)),
                     deviceBuffer, deviceSubBuffer, divider, inputSize, err);
  checkError(err, "min_max_kernel call recurse unified");

  err = queue.finish();
  print_prof_info(evt);
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
      kernel_man.getKernel<enKernel::en_min_max_kernel>());

  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR, memoryByteSizeSub);

  const cl::NDRange global_work_size(memorySizeSub);

  auto evt =
      min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                     deviceBuffer, deviceSubBuffer, divider);

  err = queue.finish();
  print_prof_info(evt);
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
      kernel_man.getKernel<enKernel::en_min_max_kernel>());

  const auto memoryByteSizeSub = memorySizeSub * sizeof(cl_float2);

  cl::Buffer deviceSubBuffer(ctx, CL_MEM_READ_WRITE, memoryByteSizeSub);
  const cl::Buffer devicePinnedSubBuffer(ctx, CL_MEM_ALLOC_HOST_PTR,
                                         memoryByteSizeSub);
  auto *sub_buffer = (cl_float2 *)queue.enqueueMapBuffer(
      devicePinnedSubBuffer, CL_TRUE, CL_MAP_READ, 0, memoryByteSizeSub);

  const cl::NDRange global_work_size(memorySizeSub);

  auto evt =
      min_max_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                     deviceBuffer, deviceSubBuffer, divider);

  err = queue.finish();
  print_prof_info(evt);
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
          kernel_man.getKernel<enKernel::en_min_max_kernel2>());

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
  auto evt =
      min_max_kernel(cl::EnqueueArgs(queue, cl::NDRange(global_work_size),
                                     cl::NDRange(local_work_size)),
                     deviceBuffer, deviceSubBuffer, divider, inputSize, err);
  checkError(err, "min_max_kernel call recurse pinned");

  err = queue.finish();
  print_prof_info(evt);
  checkError(err, "Finish min_max_recurse queue");

  /* Fetch results of calculations. */
  err = queue.enqueueReadBuffer(deviceSubBuffer, CL_TRUE, 0, memoryByteSizeSub,
                                sub_buffer);
  checkError(err, "read sub buffer min_max_recurse queue");
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
      kernel_man.getKernel<enKernel::en_crop_by_struct_kernel>());

  const auto local_work_size =
      get_local_work_size_large(cl::Device::getDefault());
  const cl::NDRange global_work_size(memorySizeInput);

  std::vector<cl::Buffer> deviceBuffer(inputSize[2]);
  std::vector<cl::Event> events(inputSize[2]);
  for (size_t slice = 0; slice < inputSize[2]; ++slice) {
    auto *buffer = ct_image->GetBufferPointer() + slice * memorySizeInput;

    /* Prepare OpenCL memory objects and place data inside them. */
    deviceBuffer.at(slice) =
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

    events.push_back(
        crop_kernel(cl::EnqueueArgs(queue, global_work_size, local_work_size),
                    deviceBuffer.at(slice), d_struct_buffer, n_coords, in_size,
                    in_orig, in_spacing));
  }

  err = queue.finish();
  std::for_each(events.begin(), events.end(), print_prof_info);
  checkError(err, "Finish crop by struct queue");

  for (auto slice = 0ul; slice < inputSize[2]; ++slice) {
    auto *buffer = ct_image->GetBufferPointer() + slice * memorySizeInput;
    /* Fetch results of calculations. */
    err = cl::copy(queue, deviceBuffer.at(slice), buffer,
                   buffer + memorySizeInput);

    checkError(err, "Add, copy data back from device");
  }
}

FloatImage2DType::Pointer OpenCL_LogItoI_subtract_median_ItoLogI(
    const FloatImage2DType::Pointer &proj_raw,
    const FloatImage2DType::Pointer &proj_scatter_intensity,
    const unsigned int median_radius) {
  const auto inputSize = proj_raw->GetBufferedRegion().GetSize();
  const cl_ulong2 in_size = {{inputSize[0], inputSize[1]}};
  if (inputSize != proj_scatter_intensity->GetBufferedRegion().GetSize()) {
    std::cerr << "Raw proj and scatter map was not the same size!\n";
    return nullptr;
  }
  if (std::pow(median_radius * 2, 2) > 64) {
    std::cerr
        << "Oh no, the kernel doesn't support a median radius this big!\n";
    return nullptr;
  }

  const auto raw_buffer = proj_raw->GetBufferPointer();
  const auto sca_buffer = proj_scatter_intensity->GetBufferPointer();

  const auto memorySizeInput = inputSize[0] * inputSize[1];

  auto defines = std::string("");
  kernel_man.initialize(defines, 2 * memorySizeInput * sizeof(cl_float));

  auto err = CL_SUCCESS;
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer_raw =
      cl::Buffer(raw_buffer, raw_buffer + memorySizeInput, true, true, &err);
  checkError(err,
             "Alloc device raw buffer, log_i_to_i_subtract_median_i_to_log_i");

  const auto deviceBuffer_sca =
      cl::Buffer(sca_buffer, sca_buffer + memorySizeInput, true, true, &err);
  checkError(
      err,
      "Alloc device scatter buffer, log_i_to_i_subtract_median_i_to_log_i");

  auto proj_corr = FloatImage2DType::New();
  proj_corr->SetRegions(proj_raw->GetBufferedRegion());
  proj_corr->Allocate();

  auto out_buffer = proj_corr->GetBufferPointer();
  auto deviceOutBuffer =
      cl::Buffer(out_buffer, out_buffer + memorySizeInput, false, true, &err);
  checkError(
      err, "Alloc device output buffer, log_i_to_i_subtract_median_i_to_log_i");

  auto log_i_to_i_subtract_median_i_to_log_i_kernel =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer, cl_ulong2, cl_int>(
          kernel_man.getKernel<enKernel::en_LogItoI_subtract_median_ItoLogI>());

  const auto local_work_size =
      get_local_work_size_small(cl::Device::getDefault());
  const auto global_work_size = cl::NDRange(memorySizeInput);

  auto evt = log_i_to_i_subtract_median_i_to_log_i_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size),
      deviceBuffer_raw, deviceBuffer_sca, deviceOutBuffer, in_size,
      static_cast<cl_int>(median_radius), err);
  checkError(err,
             "Enqueue kernel and args, log_i_to_i_subtract_median_i_to_log_i");

  err = queue.finish();
  print_prof_info(evt);
  checkError(err, "Finish log_i_to_i_subtract_median_i_to_log_i queue");

  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceOutBuffer, out_buffer,
                 out_buffer + memorySizeInput);
  checkError(
      err, "log_i_to_i_subtract_median_i_to_log_i, copy data back from device");

  return proj_corr;
}

#ifdef LOWPASS_FFT

int divisible_by_235_const(const int size) {
  auto input_size = size;
  auto ok = true;
  while (ok) {
    if (input_size % 2 == 0) {
      // ok so far
      input_size /= 2; // compiler should optimize so division and modulo is
                       // calculated simultaneously
    } else if (input_size % 3 == 0) {
      // ok so far
      input_size /= 3;
    } else if (input_size % 5 == 0) {
      // ok so far
      input_size /= 5;
    } else {
      ok = false;
    }
  }
  return input_size;
}

int get_padding(const int input_size) {
  auto cur_padding = 0;
  while (true) {
    // Padding necessary
    if (divisible_by_235_const(input_size + cur_padding) != 1) {
      // While is broken if divisible, else more padding needed.
      cur_padding += 2;
    } else {
      break;
    }
  }
  return cur_padding;
}

// template<typename T>
FloatImage2DType::Pointer gaussian_filter(FloatImage2DType::Pointer input,
                                          const double sigmaValue) {
  using RealImageType = FloatImage2DType;
  using GaussianSourceType = itk::GaussianImageSource<RealImageType>;
  using PadFilterType = itk::WrapPadImageFilter<RealImageType, RealImageType>;
  using ForwardFFTFilterType = itk::ForwardFFTImageFilter<RealImageType>;
  using ComplexImageType = ForwardFFTFilterType::OutputImageType;
  using FFTShiftFilterType =
      itk::FFTShiftImageFilter<RealImageType, RealImageType>;
  using MultiplyFilterType =
      itk::MultiplyImageFilter<ComplexImageType, RealImageType,
                               ComplexImageType>;
  using InverseFilterType =
      itk::InverseFFTImageFilter<ComplexImageType, RealImageType>;

  const auto Dimension = input->GetImageDimension();

  // Some FFT filter implementations, like VNL's, need the image size to be a
  // multiple of small prime numbers.
  auto padFilter = PadFilterType::New();
  padFilter->SetInput(input);
  PadFilterType::SizeType padding{};
  auto input_size = input->GetLargestPossibleRegion().GetSize();

  for (size_t dim = 0; dim < Dimension; dim++) {
    padding[dim] = get_padding(input_size[dim]);
    // Even though the size is usually 512 or 384 which gives padding 0,
    // it would not really speed up the code much to have the checks for these
    // sizes.
  }
  padFilter->SetPadUpperBound(padding);

  auto forwardFFTFilter = ForwardFFTFilterType::New();
  forwardFFTFilter->SetInput(padFilter->GetOutput());
  forwardFFTFilter->UpdateOutputInformation();

  // A Gaussian is used here to create a low-pass filter.
  auto gaussianSource = GaussianSourceType::New();
  gaussianSource->SetNormalized(false);
  gaussianSource->SetScale(1.0);
  const ComplexImageType::ConstPointer transformedInput =
      forwardFFTFilter->GetOutput();
  const auto inputRegion(transformedInput->GetLargestPossibleRegion());

  const auto inputSize = inputRegion.GetSize();
  gaussianSource->SetSize(inputSize);
  const auto inputSpacing = transformedInput->GetSpacing();
  gaussianSource->SetSpacing(inputSpacing);
  const auto inputOrigin = transformedInput->GetOrigin();
  gaussianSource->SetOrigin(inputOrigin);
  const auto inputDirection = transformedInput->GetDirection();
  gaussianSource->SetDirection(inputDirection);

  GaussianSourceType::ArrayType sigma;
  GaussianSourceType::PointType mean;
  sigma.Fill(sigmaValue);
  for (size_t ii = 0; ii < Dimension; ++ii) {
    const auto halfLength = inputSize[ii] * inputSpacing[ii] / 2.0;
    sigma[ii] *= halfLength;
    mean[ii] = inputOrigin[ii] + halfLength;
  }
  mean = inputDirection * mean;
  gaussianSource->SetSigma(sigma);
  gaussianSource->SetMean(mean);

  auto fftShiftFilter = FFTShiftFilterType::New();
  fftShiftFilter->SetInput(gaussianSource->GetOutput());

  auto multiplyFilter = MultiplyFilterType::New();
  multiplyFilter->SetInput1(forwardFFTFilter->GetOutput());
  multiplyFilter->SetInput2(fftShiftFilter->GetOutput());

  auto inverseFFTFilter = InverseFilterType::New();
  inverseFFTFilter->SetInput(multiplyFilter->GetOutput());
  inverseFFTFilter->Update();
  return inverseFFTFilter->GetOutput();
}
#else
FloatImage2DType::Pointer gaussian_filter(FloatImage2DType::Pointer input,
                                          const double sigmaValue) {
  using SmoothingFilterType =
      itk::SmoothingRecursiveGaussianImageFilter<FloatImage2DType,
                                                 FloatImage2DType>;
  SmoothingFilterType::Pointer gaussianFilter = SmoothingFilterType::New();
  gaussianFilter->SetInput(input);
  const auto input_size = input->GetBufferedRegion().GetSize();
  if (input_size[0] == input_size[1]) {
    gaussianFilter->SetSigma(
        sigmaValue); // filter specific setting for 512x 512 image
  } else {
    SmoothingFilterType::SigmaArrayType gaussianSigmaArray;
    gaussianSigmaArray[0] = sigmaValue;
    gaussianSigmaArray[1] = .75 * sigmaValue;
    gaussianFilter->SetSigmaArray(
        gaussianSigmaArray); // filter specific setting for 512x384 (varian/2)
  }
  gaussianFilter->Update();
  return gaussianFilter->GetOutput();
}
#endif

FloatImage2DType::Pointer OpenCL_LogItoI_subtract_median_gaussian(
    const FloatImage2DType::Pointer &proj_raw,
    const FloatImage2DType::Pointer &proj_prim,
    const unsigned int median_radius, const float gaussian_sigma) {

  const auto inputSize = proj_raw->GetBufferedRegion().GetSize();
  const cl_ulong2 in_size = {{inputSize[0], inputSize[1]}};

  if (inputSize != proj_prim->GetBufferedRegion().GetSize()) {
    std::cerr << "Raw proj and scatter map was not the same size!\n";
    return nullptr;
  }
  if (median_radius * 2 + 64 > 128) {
    std::cerr
        << "Oh no, the kernel doesn't support a median radius this big!\n";
    return nullptr;
  }

  const auto raw_buffer = proj_raw->GetBufferPointer();
  const auto pri_buffer = proj_prim->GetBufferPointer();

  const auto memorySizeInput = inputSize[0] * inputSize[1];

  auto defines = std::string("");
  kernel_man.initialize(defines, 2 * memorySizeInput * sizeof(cl_float));

  auto err = CL_SUCCESS;
  const auto ctx = cl::Context::getDefault(&err);
  checkError(err, "Get default context");
  auto queue = cl::CommandQueue::getDefault(&err);
  checkError(err, "Get default queue");

  /* Prepare OpenCL memory objects and place data inside them. */
  const auto deviceBuffer_raw = cl::Buffer(
      ctx, raw_buffer, raw_buffer + memorySizeInput, true, true, nullptr);

  const auto deviceBuffer_pri = cl::Buffer(
      ctx, pri_buffer, pri_buffer + memorySizeInput, true, true, nullptr);

  auto proj_sca = FloatImage2DType::New();
  proj_sca->SetRegions(proj_raw->GetBufferedRegion());
  proj_sca->Allocate();

  auto sca_buffer = proj_sca->GetBufferPointer();
  auto deviceBuffer_sca = cl::Buffer(
      ctx, sca_buffer, sca_buffer + memorySizeInput, false, true, nullptr);

  auto log_i_to_i_subtract_median_y_x_kernel =
      cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer, cl_ulong2, cl_int>(
          kernel_man.getKernel<enKernel::en_LogItoI_subtract_median_y_x>());

  auto local_work_size = 128U; // bench test what is optimal
  while (memorySizeInput % local_work_size != 0) {
    local_work_size /= 2;
  }
  const auto actual_work_size = local_work_size - 2 * median_radius;
  const auto n_work_groups = memorySizeInput / actual_work_size;
  const auto super_global_work_size =
      cl::NDRange(local_work_size * n_work_groups);
  const auto super_local_work_size = cl::NDRange(local_work_size);

  auto evt = log_i_to_i_subtract_median_y_x_kernel(
      cl::EnqueueArgs(queue, super_global_work_size, super_local_work_size),
      deviceBuffer_raw, deviceBuffer_pri, deviceBuffer_sca, in_size,
      static_cast<cl_int>(median_radius));

  err = queue.finish();
  print_prof_info(evt);
  checkError(err, "Finish log_i_to_i_subtract_median_y_x queue");

  // Gaussian filter
  // We could use pinned/unified memory or SVM here, but
  // the allocations may be more expensive than the copy.
  /* Fetch results of calculations. */
  err = cl::copy(queue, deviceBuffer_sca, sca_buffer,
                 sca_buffer + memorySizeInput);
  checkError(err, "log_i_to_i_subtract_median, copy data back from device");

  proj_sca = gaussian_filter(proj_sca, gaussian_sigma);

  /* No need to convert back, see GenScatterMap_PrioriCT
  auto out_buffer = proj_sca->GetBufferPointer();
  auto deviceBuffer_out = cl::Buffer(
      ctx, out_buffer, out_buffer + memorySizeInput, false, true, nullptr);

  // And finally transform i back to log i
  auto i_to_log_i_kernel = cl::KernelFunctor<cl::Buffer, cl_ulong2>(
      kernel_man.getKernel<enKernel::en_ItoLogI>());

  const auto local_work_size_small =
      get_local_work_size_small(cl::Device::getDefault());
  const auto global_work_size = cl::NDRange(memorySizeInput);
  i_to_log_i_kernel(
      cl::EnqueueArgs(queue, global_work_size, local_work_size_small),
      deviceBuffer_out, in_size);

  err = queue.finish();
print_prof_info(evt);
  checkError(err, "Finish i_to_log_i queue");

  // Fetch results of calculations.
  err = cl::copy(queue, deviceBuffer_out, out_buffer,
                 out_buffer + memorySizeInput);
  checkError(err, "Add, copy data back from device");
  */
  return proj_sca;
}
