#ifndef IMAGEFILTERS_H
#define IMAGEFILTERS_H

#include <vector>

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

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

#include "itkImage.h"

#include "OpenCL/cl2.hpp"

struct Rtss_roi_modern;

CBCTRECON_API void OpenCL_initialize(size_t required_mem_alloc_size,
                                     std::string &defines);
CBCTRECON_API void OpenCL_initialize(std::string &defines);

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
itk::Image<float, 3U>::Pointer OpenCL_divide3Dby3D_loginv_OutOfPlace(
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
                                 cl::Buffer &deviceBuffer);

CBCTRECON_API
void OpenCL_crop_by_struct_InPlace(UShortImageType::Pointer &ct_image,
                                   const Rtss_roi_modern &voi);

CBCTRECON_API
FloatImage2DType::Pointer OpenCL_LogItoI_subtract_median_ItoLogI(
    const FloatImage2DType::Pointer &proj_raw,
    const FloatImage2DType::Pointer &proj_scatter_intensity,
    const unsigned int median_radius);

CBCTRECON_API
FloatImage2DType::Pointer OpenCL_LogItoI_subtract_median_gaussian(
    const FloatImage2DType::Pointer &proj_raw,
    const FloatImage2DType::Pointer &proj_prim,
    const unsigned int median_radius, const float gaussian_sigma);

#endif // IMAGEFILTERS_H
