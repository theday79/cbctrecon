#include <itkImage.h>
#include <vector_types.h>
// #include <complex>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#ifndef OpenCLFFTFilter_H_
#define OpenCLFFTFilter_H_

void OpenCL_padding(
	const cl_int4 paddingIndex,
	const cl_uint4 paddingSize,
	const cl_uint4 inputSize,
	const cl_float *deviceVolume,
	cl_float *devicePaddedVolume,
	const std::vector<cl_float> mirrorWeights);

void OpenCL_fft_convolution(const cl_int4 inputDimension,
	const cl_int2 kernelDimension,
	float *hostProjection,
	cl_float2 *deviceKernelFFT);

void OpenCL_subtract3Dfrom2DbySlice_InPlace(
	cl_float* input, 
	const cl_float* subimage, 
	itk::Image<float, 3U>::SizeType inputSize, 
	itk::Image<float, 2U>::SizeType subSize);

// Actually is divide ln(65535/X) by ln(65535/Y)
itk::Image<float, 3U>::Pointer OpenCL_divide3Dby3D_OutOfPlace(
	itk::Image<unsigned short, 3U>::Pointer Num3D,
	itk::Image<unsigned short, 3U>::Pointer Denum3D);

void OpenCL_AddConst_InPlace(
	cl_float* input,
	itk::Image<float, 3U>::SizeType inputSize,
	const cl_float constant);

void OpenCL_AddConst_InPlace_2D(
	cl_float* input,
	itk::Image<float, 2U>::SizeType inputSize,
	const cl_float constant);

cl_float2 OpenCL_min_max(
	const cl_float* input,
	itk::Image<float, 3U>::SizeType inputSize);

cl_float2 OpenCL_min_max_2D(
	const cl_float* input,
	itk::Image<float, 2U>::SizeType inputSize);

cl_float2 OpenCL_min_max_recurse(
	const cl_float2* input,
	const cl_uint inputSize);

#endif OpenCLFFTFilter_H_
