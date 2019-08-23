#ifndef IMAGEFILTERS_H
#define IMAGEFILTERS_H

#include <vector>

#if CBCTRECON_OPENCL_VERSION >= 210
#define CL_HPP_MINIMUM_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 210
#else
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#endif

#include "OpenCL/cl2.hpp"

#include "itkImage.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

// If this function fails, something is probably wrong with the local OpenCL ICD
auto getDeviceByReqAllocSize(size_t required_mem_alloc_size);

CBCTRECON_API void OpenCL_initialize(size_t required_mem_alloc_size);
CBCTRECON_API void OpenCL_initialize(cl::Device &dev);

CBCTRECON_API
void OpenCL_padding(const cl_int4 &paddingIndex, const cl_uint4 &paddingSize,
                    const cl_uint4 &inputSize, const float *hostVolume,
                    float *hostPaddedVolume,
                    const std::vector<float> &mirrorWeights);

CBCTRECON_API
void OpenCL_subtract2Dfrom3DbySlice_InPlace(
    FloatImageType::Pointer &projections,
    const FloatImage2DType::Pointer &filter);

// Actually is divide ln(65535/X) by ln(65535/Y)
CBCTRECON_API
itk::Image<float, 3U>::Pointer OpenCL_divide3Dby3D_OutOfPlace(
    const itk::Image<unsigned short, 3U>::Pointer &Num3D,
    const itk::Image<unsigned short, 3U>::Pointer &Denum3D);

CBCTRECON_API
void OpenCL_AddConst_InPlace(cl_float *buffer,
                             const FloatImageType::SizeType &inputSize,
                             cl_float constant);

CBCTRECON_API
void OpenCL_AddConst_MulConst_InPlace(cl_float *buffer,
                                      const FloatImageType::SizeType &inputSize,
                                      cl_float add_constant,
                                      cl_float mul_constant);

CBCTRECON_API
void OpenCL_AddConst_InPlace_2D(cl_float *buffer,
                                const FloatImage2DType::SizeType &inputSize,
                                cl_float constant);

CBCTRECON_API
cl_float2 OpenCL_min_max_1D(cl_float *buffer, size_t memorySizeInput);

CBCTRECON_API
cl_float2 OpenCL_min_max_3D(cl_float *buffer,
                            const FloatImageType::SizeType &inputSize);

CBCTRECON_API
cl_float2 OpenCL_min_max_2D(cl_float *buffer,
                            const FloatImage2DType::SizeType &inputSize);

CBCTRECON_API
cl_float2 OpenCL_min_max_recurse(cl_float2 *buffer, size_t inputSize,
                                 cl::Buffer &deviceBuffer,
                                 cl::CommandQueue &queue,
                                 cl::NDRange nd_local_work_size);

#endif // IMAGEFILTERS_H
