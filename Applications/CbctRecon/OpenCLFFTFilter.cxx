#include "OpenCLFFTFilter.h"
#include <itkImage.h>
#include <complex>
#include <iostream>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <clFFT.h>
#include "rtkOpenCLUtilities.h"

std::tuple<cl_platform_id, cl_device_id>
getPlatformAndDeviceID(const size_t required_mem_alloc_size) {
	cl_platform_id platform;
	cl_device_id   device;
	cl_int         err = CL_SUCCESS;

	/* Setup OpenCL environment. */
	err = clGetPlatformIDs(1, &platform, NULL);

	err |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, NULL);

	cl_ulong max_mem_alloc;
	err |= clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &max_mem_alloc, NULL);

	if (required_mem_alloc_size > max_mem_alloc) {
		std::cout << "Memory required: " << required_mem_alloc_size / (1024 * 1024) << "MB"
			<< " was larger than DEV_MAX_MEM: " << max_mem_alloc / (1024 * 1024) << "MB\n"
			<< "Changing to CPU mode..." << std::endl;
		err |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, NULL);
	}

	return std::make_tuple(platform, device);
}



void
OpenCL_padding(
	const cl_int4 paddingIndex,
	const cl_uint4 paddingSize,
	const cl_uint4 inputSize,
	const cl_float *hostVolume, //input
	cl_float *hostPaddedVolume, //output
	const std::vector<cl_float> mirrorWeights)
{
	cl_int           err = CL_SUCCESS;
	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_program       m_Program;
	cl_kernel        m_Kernel;

	size_t pv_buffer_size = paddingSize.x * paddingSize.y * paddingSize.z * sizeof(cl_float);
	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(pv_buffer_size);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	// std::cout << "Padding OpenCL style" << std::endl;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not create OpenCL context, error code: " << err << std::endl;
	queue = clCreateCommandQueue(ctx, device, 0, &err);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not create OpenCL queue, error code: " << err << std::endl;


	/* Prepare OpenCL memory objects and place data inside them. */
	size_t w_buffer_size = mirrorWeights.size() * sizeof(cl_float);
	const cl_uint w_buf_sizeof = mirrorWeights.size(); // just because intellisense couldn't understand it below...
	cl_mem weights_d = clCreateBuffer(ctx, CL_MEM_READ_ONLY, w_buffer_size, NULL, &err);
	err = clEnqueueWriteBuffer(queue, weights_d, CL_TRUE, 0, w_buffer_size, &(mirrorWeights[0]), 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not write OpenCL weigths_d buffer, error code: " << err << std::endl;


	cl_mem devicePaddedVolume = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, pv_buffer_size, hostPaddedVolume, &err);
	// err = clEnqueueWriteBuffer(queue, devicePaddedVolume, CL_TRUE, 0, pv_buffer_size, hostPaddedVolume, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not write OpenCL devicePaddedVolume buffer, error code: " << err << std::endl;

	size_t v_buffer_size = sizeof(hostVolume);
	cl_mem deviceVolume = clCreateBuffer(ctx, CL_MEM_READ_ONLY, v_buffer_size, NULL, &err);
	err = clEnqueueWriteBuffer(queue, deviceVolume, CL_TRUE, 0, v_buffer_size, hostVolume, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not write OpenCL deviceVolume buffer, error code: " << err << std::endl;


	// some text
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl",
		ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "padding_kernel", &err);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not create OpenCL kernel, error code: " << err << std::endl;

	const size_t local_work_size = 128;

	const size_t global_work_size = paddingSize.x * paddingSize.y * paddingSize.z;

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceVolume);
	err |= clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&devicePaddedVolume);
	err |= clSetKernelArg(m_Kernel, 2, sizeof(cl_int4), &paddingIndex);
	err |= clSetKernelArg(m_Kernel, 3, sizeof(cl_uint4), &paddingSize);
	err |= clSetKernelArg(m_Kernel, 4, sizeof(cl_uint4), &inputSize);
	err |= clSetKernelArg(m_Kernel, 5, sizeof(cl_mem), (void *)&weights_d);
	err |= clSetKernelArg(m_Kernel, 6, sizeof(cl_uint), &w_buf_sizeof);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not parse OpenCL arguments, error code: " << err << std::endl;

	// Call to kernel

	// Execute kernel
	cl_event events[2];
	err = clEnqueueNDRangeKernel(queue,	m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not wait for OpenCL event!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, devicePaddedVolume, CL_TRUE, 0, pv_buffer_size, hostPaddedVolume, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "PAD::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL memory objects. */
	clReleaseMemObject(weights_d);
	clReleaseMemObject(deviceVolume);
	clReleaseMemObject(devicePaddedVolume);

	/* Release OpenCL working objects. */
	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);
	// std::cout << "Did it work? we just don't know..." << std::endl;

}

void clFFT_Forward_Multiply_Backward(float *hostProjection, const cl_float2* hostKernelFFT, const cl_int4 inputDimension, const cl_int2 kernelDimension) {
	cl_int err;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
	cl_context ctx = 0;
	cl_command_queue queue = 0;


	/* FFT library realted declarations */
	clfftPlanHandle planHandle;
	clfftDim dim = CLFFT_3D;

	const size_t memorySizeProjection = inputDimension.x  * inputDimension.y  * inputDimension.z * sizeof(float);
	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeProjection);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	/* Setup clFFT. */
	clfftSetupData fftSetupFwd;
	err = clfftInitSetupData(&fftSetupFwd);
	err = clfftSetup(&fftSetupFwd);

	// device pointers
	cl_mem  deviceProjectionFFT; // outputpointer
	cl_mem  deviceProjection; // copy of host
	cl_int4 fftDimension = inputDimension;
	fftDimension.x = inputDimension.x / 2 + 1;

	const size_t memorySizeProjectionFFT = fftDimension.x    * fftDimension.y    * fftDimension.z * sizeof(cl_float2);
	//std::cout << "Creating buffers..." << std::endl;
	/* Prepare OpenCL memory objects and place data inside them. */
	deviceProjection = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeProjection, NULL, &err);
	deviceProjectionFFT = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, memorySizeProjectionFFT, NULL, &err);

	// we don't have any data to put into deviceProjectionFFT, so we fill it with zeros:
	cl_float2* hostProjectionFFT = new cl_float2[memorySizeProjectionFFT];

	err = clEnqueueWriteBuffer(queue, deviceProjection, CL_TRUE, 0, memorySizeProjection, hostProjection, 0, NULL, NULL);
	err |= clEnqueueWriteBuffer(queue, deviceProjectionFFT, CL_TRUE, 0, memorySizeProjectionFFT, hostProjectionFFT, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "FWD::Could not create OpenCL buffers, error code: " << err << std::endl;

	size_t clLengths[3] = { (size_t)inputDimension.z, (size_t)inputDimension.y, (size_t)inputDimension.x };

	/* Create a default plan for a complex FFT. */
	err = clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);
	if (err != CL_SUCCESS)
		std::cout << "FWD::Could not create default plan, error code: " << err << std::endl;

	/* Set plan parameters. */
	err = clfftSetPlanPrecision(planHandle, CLFFT_SINGLE);
	err |= clfftSetLayout(planHandle, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
	err |= clfftSetResultLocation(planHandle, CLFFT_OUTOFPLACE);
	if (err != CL_SUCCESS)
		std::cout << "FWD::Could not set up plan, error code: " << err << std::endl;

	/* Bake the plan. */
	err |= clfftBakePlan(planHandle, 1, &queue, NULL, NULL);

	/* Execute the plan. */
	err |= clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &queue, 0, NULL, NULL, &deviceProjection, &deviceProjectionFFT, NULL);

	if (err != CL_SUCCESS)
		std::cout << "FWD::Could not perform forward fourier transform, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clFinish(queue);

	if (err != CL_SUCCESS)
		std::cout << "FWD::Could not finish OpenCL queue, error code: " << err << std::endl;

	const size_t memorySizeKernelFFT = kernelDimension.x * kernelDimension.y * sizeof(cl_float2);

	/* Prepare OpenCL memory objects and place data inside them. */
	cl_mem deviceKernelFFT = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeKernelFFT, NULL, &err);

	err = clEnqueueWriteBuffer(queue, deviceProjectionFFT, CL_TRUE, 0, memorySizeProjectionFFT, hostProjectionFFT, 0, NULL, NULL);
	err |= clEnqueueWriteBuffer(queue, deviceKernelFFT, CL_TRUE, 0, memorySizeKernelFFT, hostKernelFFT, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "MUL::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	cl_program m_Program;
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	cl_kernel m_Kernel;
	if (kernelDimension.y == 1)
		m_Kernel = clCreateKernel(m_Program, "multiply_kernel", &err);
	else
		m_Kernel = clCreateKernel(m_Program, "multiply_kernel2D", &err);

	if (err != CL_SUCCESS)
		std::cout << "MUL::Could not create OpenCL kernel, error code: " << err << std::endl;

	cl_event events[2];
	const size_t local_work_size = 128;

	const size_t global_work_size = (size_t)(fftDimension.x * fftDimension.y * fftDimension.z);

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceProjectionFFT);
	err |= clSetKernelArg(m_Kernel, 1, sizeof(cl_int4), &fftDimension);
	err |= clSetKernelArg(m_Kernel, 2, sizeof(cl_mem), (void *)&deviceKernelFFT);
	if (err != CL_SUCCESS)
		std::cout << "MUL::Could not parse args to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	err |= clWaitForEvents(1, &events[0]);
	err |= clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MUL::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "MUL::Could not finish OpenCL queue, error code: " << err << std::endl;

	//DONE resetting queue

	/* Setup clFFT. */
	clfftSetupData fftSetupBck;
	err = clfftInitSetupData(&fftSetupBck);
	err |= clfftSetup(&fftSetupBck);

	/* Create a default plan for a complex FFT. */
	err |= clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);

	/* Set plan parameters. */
	err |= clfftSetPlanPrecision(planHandle, CLFFT_SINGLE);
	err |= clfftSetLayout(planHandle, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
	err |= clfftSetResultLocation(planHandle, CLFFT_OUTOFPLACE);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not set OpenCL plan parameters, error code: " << err << std::endl;

	/* Bake the plan. */
	err = clfftBakePlan(planHandle, 1, &queue, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not bake OpenCL kernel, error code: " << err << std::endl;

	/* Execute the plan. */
	err = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &queue, 0, NULL, NULL, &deviceProjectionFFT, &deviceProjection, NULL);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not create OpenCL kernel, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceProjection, CL_TRUE, 0, memorySizeProjection, hostProjection, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not read OpenCL output, error code: " << err << std::endl;

	/* Release OpenCL memory objects. */
	clReleaseMemObject(deviceProjection);
	clReleaseMemObject(deviceProjectionFFT);
	clReleaseMemObject(deviceKernelFFT);

	/* Release the plan. */
	err = clfftDestroyPlan(&planHandle);

	/* Release clFFT library. */
	clfftTeardown();

	/* Release OpenCL working objects. */
	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not release OpenCL kernel, error code: " << err << std::endl;
}


void clFFT_backwards(float* hostProjection, const cl_float2* hostProjectionFFT, const cl_int4 inputDimension) {

	// INV FFT START
	cl_int err;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
	cl_context ctx = 0;
	cl_command_queue queue = 0;

	/* FFT library realted declarations */
	clfftPlanHandle planHandle;
	clfftDim dim = CLFFT_3D; 

	const size_t memorySizeProjection = inputDimension.x  * inputDimension.y  * inputDimension.z * sizeof(float);
	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeProjection);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	//std::cout << "Creating buffers..." << std::endl;

	cl_int4 fftDimension = inputDimension;
	fftDimension.x = inputDimension.x / 2 + 1;

	const size_t memorySizeProjectionFFT = fftDimension.x    * fftDimension.y    * fftDimension.z   * sizeof(cl_float2);

	/* Prepare OpenCL memory objects and place data inside them. */
	cl_mem deviceProjectionFFT = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  memorySizeProjectionFFT, NULL, &err);
	cl_mem deviceProjection    = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeProjection, NULL, &err);

	err = clEnqueueWriteBuffer(queue, deviceProjectionFFT, CL_TRUE, 0, memorySizeProjectionFFT, hostProjectionFFT, 0, NULL, NULL);
	err = clEnqueueWriteBuffer(queue, deviceProjection,    CL_TRUE, 0, memorySizeProjection,    hostProjection, 0, NULL, NULL);

	//std::cout << "Buffers created..." << std::endl;

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could create OpenCL buffers, error code: " << err << std::endl;

	size_t clLengths[3] = { (size_t)inputDimension.z, (size_t)inputDimension.y, (size_t)inputDimension.x };

	//DONE resetting queue

	/* Setup clFFT. */
	clfftSetupData fftSetupBck;
	err = clfftInitSetupData(&fftSetupBck);
	err |= clfftSetup(&fftSetupBck);


	/* Create a default plan for a complex FFT. */
	err |= clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);

	/* Set plan parameters. */
	err |= clfftSetPlanPrecision(planHandle, CLFFT_SINGLE);
	err |= clfftSetLayout(planHandle, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
	err |= clfftSetResultLocation(planHandle, CLFFT_OUTOFPLACE);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not set OpenCL plan parameters, error code: " << err << std::endl;

	/* Bake the plan. */
	err = clfftBakePlan(planHandle, 1, &queue, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not bake OpenCL kernel, error code: " << err << std::endl;

	/* Execute the plan. */
	err = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &queue, 0, NULL, NULL, &deviceProjectionFFT, &deviceProjection, NULL);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not create OpenCL kernel, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceProjection, CL_TRUE, 0, memorySizeProjection, hostProjection, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not read OpenCL output, error code: " << err << std::endl;
	
	/* Release OpenCL memory objects. */
	clReleaseMemObject(deviceProjection);
	clReleaseMemObject(deviceProjectionFFT);

	/* Release the plan. */
	err = clfftDestroyPlan(&planHandle);

	/* Release clFFT library. */
	clfftTeardown();

	/* Release OpenCL working objects. */
	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "BCK::Could not release OpenCL kernel, error code: " << err << std::endl;

}

void
OpenCL_fft_convolution(const cl_int4 inputDimension,
	const cl_int2 kernelDimension,
	float *hostProjection,
	cl_float2 *hostKernelFFT)
{
	//std::cout << "Forward projecting..." << std::endl;
	clFFT_Forward_Multiply_Backward(hostProjection, hostKernelFFT, inputDimension, kernelDimension);

	//OpenCL_multiply_by_kernel(hostProjectionFFT, hostKernelFFT, inputDimension, kernelDimension);

	//clFFT_backwards(hostProjection, hostProjectionFFT, inputDimension);

	printf(".");
}



void OpenCL_subtract3Dfrom2DbySlice_InPlace(cl_float* buffer, const cl_float* sub_buffer, 
	itk::Image<float, 3U>::SizeType inputSize, itk::Image<float, 2U>::SizeType subSize) {
	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;
	cl_mem  deviceSubBuffer;

	const size_t memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2] * sizeof(cl_float);
	const size_t memorySizeSub = subSize[0] * subSize[1] * sizeof(cl_float);

	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not setup OpenCL, error code: " << err << std::endl;

	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeInput, NULL, &err);
	deviceSubBuffer = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeSub, NULL, &err);

	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);
	err = clEnqueueWriteBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "SUB::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "subtract_kernel2D", &err);

	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not create OpenCL kernel, error code: " << err << std::endl;

	cl_event events[2];
	const size_t local_work_size = 128;

	const size_t global_work_size = inputSize[0]* inputSize[1] * inputSize[2]; 

	const cl_uint4 inputDim = {
		(cl_uint) inputSize[0],
		(cl_uint) inputSize[1],
		(cl_uint) inputSize[2],
		0 };

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_int4), &inputDim);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 2, sizeof(cl_mem), (void *)&deviceSubBuffer);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not parse arg 3 to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceSubBuffer);
	clReleaseMemObject(deviceBuffer);

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	// std::cout << "Did it work? we just don't know..." << std::endl; IT WORKS
	if (err != CL_SUCCESS)
		std::cout << "SUB::Could not create OpenCL kernel, error code: " << err << std::endl;

}


itk::Image<float, 3U>::Pointer OpenCL_divide3Dby3D_OutOfPlace(
	itk::Image<unsigned short, 3U>::Pointer Num3D,
	itk::Image<unsigned short, 3U>::Pointer Denum3D) {

	unsigned short* buffer = Num3D->GetBufferPointer();
	itk::Image<unsigned short, 3U>::SizeType inputSize = Num3D->GetLargestPossibleRegion().GetSize();
	unsigned short* sub_buffer = Denum3D->GetBufferPointer();
	itk::Image<unsigned short, 3U>::SizeType subSize = Denum3D->GetLargestPossibleRegion().GetSize();


	itk::Image<float, 3U>::Pointer outImage = itk::Image<float, 3U>::New();
	itk::Image<float, 3U>::SizeType projCT_size = Num3D->GetLargestPossibleRegion().GetSize();
	itk::Image<float, 3U>::IndexType projCT_idxStart = Num3D->GetLargestPossibleRegion().GetIndex();
	itk::Image<float, 3U>::SpacingType projCT_spacing = Num3D->GetSpacing();
	itk::Image<float, 3U>::PointType  projCT_origin = Num3D->GetOrigin();

	itk::Image<float, 3U>::RegionType projCT_region;
	projCT_region.SetSize(projCT_size);
	projCT_region.SetIndex(projCT_idxStart);

	outImage->SetRegions(projCT_region);
	outImage->SetSpacing(projCT_spacing);
	outImage->SetOrigin(projCT_origin);

	outImage->Allocate();

	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;
	cl_mem  deviceSubBuffer;
	cl_mem  deviceOutBuffer;

	const size_t memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2] * sizeof(cl_ushort);
	const size_t memorySizeOutput = inputSize[0] * inputSize[1] * inputSize[2] * sizeof(cl_float);
	const size_t memorySizeSub = subSize[0] * subSize[1] * subSize[2] * sizeof(cl_ushort);

	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not setup OpenCL, error code: " << err << std::endl;


	// std::cout << "SUB3D::Input xyz size: " << inputSize[0] * inputSize[1] * inputSize[2] << std::endl;
	// std::cout << "SUB3D::subIm xyz size: " << subSize[0] * subSize[1] * subSize[2] << std::endl;

	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeInput, NULL, &err);
	deviceSubBuffer = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeSub, NULL, &err);
	deviceOutBuffer = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, memorySizeOutput, NULL, &err);

	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);
	err = clEnqueueWriteBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);
	err = clEnqueueWriteBuffer(queue, deviceOutBuffer, CL_TRUE, 0, memorySizeOutput, outImage->GetBufferPointer(), 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "divide_kernel3D_Ushort", &err);

	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not create OpenCL kernel, error code: " << err << std::endl;


	cl_event events[2];
	const size_t local_work_size = 128;

	const size_t global_work_size = inputSize[0] * inputSize[1] * inputSize[2];
	const cl_uint4 inputDim = {
		(cl_uint)inputSize[0],
		(cl_uint)inputSize[1],
		(cl_uint)inputSize[2],
		0 };

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&deviceSubBuffer);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 2, sizeof(cl_mem), (void *)&deviceOutBuffer);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not parse arg 3 to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */

	err = clEnqueueReadBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, outImage->GetBufferPointer(), 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not read OpenCL buffer, error code: " << err << std::endl;
	
	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceSubBuffer);
	clReleaseMemObject(deviceOutBuffer);
	clReleaseMemObject(deviceBuffer);

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	// std::cout << "Did it work? we just don't know..." << std::endl; IT WORKS
	if (err != CL_SUCCESS)
		std::cout << "SUB3D::Could not create OpenCL kernel, error code: " << err << std::endl;

	
	return outImage;
}


void OpenCL_AddConst_InPlace(cl_float* buffer, 
	itk::Image<float, 3U>::SizeType inputSize,
	const cl_float constant) {

	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;

	const size_t memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2] * sizeof(cl_float);

	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not setup OpenCL, error code: " << err << std::endl;

	//itk::Image<float, 3U>::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();

	// std::cout << "ADD::Input xy size: " << inputSize[0] * inputSize[1] << std::endl;

	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeInput, NULL, &err);

	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	// std::cout << "build subtraction kernel..." << std::endl;
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "add_const_kernel", &err);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not create OpenCL kernel, error code: " << err << std::endl;


	std::cout << "Kernel created." << std::endl;
	cl_event events[2];
	const size_t local_work_size = 128;

	const size_t global_work_size = inputSize[0] * inputSize[1] * inputSize[2]; 

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_float), &constant);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceBuffer);

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not create OpenCL kernel, error code: " << err << std::endl;

}

void OpenCL_AddConst_InPlace_2D(cl_float* buffer,
	itk::Image<float, 2U>::SizeType inputSize,
	const cl_float constant) {

	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;

	const size_t memorySizeInput = inputSize[0] * inputSize[1] * sizeof(cl_float);
	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not setup OpenCL, error code: " << err << std::endl;

	//itk::Image<float, 3U>::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();

	// std::cout << "ADD::Input xy size: " << inputSize[0] * inputSize[1] << std::endl;

	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeInput, NULL, &err);

	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "add_const_kernel", &err);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not create OpenCL kernel, error code: " << err << std::endl;


	cl_event events[2];
	const size_t local_work_size = 128;

	const size_t global_work_size = inputSize[0] * inputSize[1];

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_float), &constant);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceBuffer);

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "ADD::Could not create OpenCL kernel, error code: " << err << std::endl;

}

cl_float2 OpenCL_min_max(
	const cl_float* buffer,
	itk::Image<float, 3U>::SizeType inputSize) {

	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;

	const size_t memorySizeInput = inputSize[0] * inputSize[1] * inputSize[2] * sizeof(cl_float);
	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not setup OpenCL, error code: " << err << std::endl;

	//std::cout << "MMM::Input xy size: " << inputSize[0] * inputSize[1] * inputSize[2] << std::endl;

	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeInput, NULL, &err);
	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);
	
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	//std::cout << "build min max kernel..." << std::endl;
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "min_max_kernel", &err);

	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not create OpenCL kernel, error code: " << err << std::endl;


	//std::cout << "Kernel created." << std::endl;
	cl_event events[2];
	const size_t local_work_size = 128;

	cl_uint divider = 128;
	while (true) {
		if (inputSize[0] % divider == 0)
			break;
		divider /= 2;
		if (divider == 2)
			break;
	}

	const cl_uint4 outputDim = {
		(cl_uint)inputSize[0] / divider,
		(cl_uint)inputSize[1],
		(cl_uint)inputSize[2],
		0 };

	const size_t global_work_size = outputDim.x * outputDim.y * outputDim.z;

	const cl_uint4 inputDim = {
		(cl_uint)inputSize[0],
		(cl_uint)inputSize[1],
		(cl_uint)inputSize[2],
		0 };


	const size_t memorySizeSub = outputDim.x * outputDim.y * outputDim.z * sizeof(cl_float2);
	//std::cout << "MMM::output xy size: " << outputDim.x * outputDim.y * outputDim.z << std::endl;
	cl_float2* sub_buffer = new cl_float2[outputDim.x * outputDim.y * outputDim.z];
	cl_mem  deviceSubBuffer;
	deviceSubBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeSub, NULL, &err);
	err = clEnqueueWriteBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);


	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&deviceSubBuffer);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 2, sizeof(cl_uint), &divider);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not parse arg 3 to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceSubBuffer);
	clReleaseMemObject(deviceBuffer);

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "MMM::Could not create OpenCL kernel, error code: " << err << std::endl;

	return OpenCL_min_max_recurse(sub_buffer, (cl_uint) global_work_size);
}

cl_float2 OpenCL_min_max_2D(
	const cl_float* buffer,
	itk::Image<float, 2U>::SizeType inputSize) {

	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;

	const size_t memorySizeInput = inputSize[0] * inputSize[1] * sizeof(cl_float);
	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not setup OpenCL, error code: " << err << std::endl;


	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeInput, NULL, &err);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not create OpenCL mem object, error code: " << err << std::endl;

	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not create OpenCL buffers, error code: " << err << std::endl; //-38

	// Create program
	//std::cout << "build min max kernel..." << std::endl;
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "min_max_kernel", &err);

	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not create OpenCL kernel, error code: " << err << std::endl;


	//std::cout << "Kernel created." << std::endl;
	cl_event events[2];
	const size_t local_work_size = 128;

	cl_uint divider = 128;
	while (true) {
		if (inputSize[0] % divider == 0)
			break;
		divider /= 2;
		if (divider == 2)
			break;
	}

	const cl_uint2 outputDim = {
		(cl_uint)inputSize[0] / divider,
		(cl_uint)inputSize[1]
	};

	const size_t global_work_size = outputDim.x * outputDim.y;

	const cl_uint2 inputDim = {
		(cl_uint)inputSize[0],
		(cl_uint)inputSize[1]
	};


	const size_t memorySizeSub = outputDim.x * outputDim.y * sizeof(cl_float2);
	//std::cout << "MMM::output xy size: " << outputDim.x * outputDim.y * outputDim.z << std::endl;
	cl_float2* sub_buffer = new cl_float2[outputDim.x * outputDim.y];
	cl_mem  deviceSubBuffer;
	deviceSubBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeSub, NULL, &err);
	err = clEnqueueWriteBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);


	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&deviceSubBuffer);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 2, sizeof(cl_uint), &divider);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not parse arg 3 to OpenCL kernel, error code: " << err << std::endl;

	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceSubBuffer);
	clReleaseMemObject(deviceBuffer);

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "MM2D::Could not create OpenCL kernel, error code: " << err << std::endl;

	return OpenCL_min_max_recurse(sub_buffer, (cl_uint)global_work_size);
}

cl_float2 OpenCL_min_max_recurse(
	const cl_float2* buffer,
	const cl_uint inputSize) {

	if (inputSize < 257) {
		cl_float2 out = { 65535.0, -9999.0 };
		for (cl_uint i = 0; i < inputSize; i++) {
			if (out.x > buffer[i].x)
				out.x = buffer[i].x;
			if (out.y < buffer[i].y)
				out.y = buffer[i].y;
		}
		return out;
	}

	cl_context       ctx = 0;
	cl_command_queue queue = 0;
	cl_int           err;

	cl_program            m_Program;
	cl_kernel             m_Kernel;
	cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };

	cl_mem  deviceBuffer;

	const size_t memorySizeInput = inputSize * sizeof(cl_float2);

	std::tuple<cl_platform_id, cl_device_id> dev_tuple =
		getPlatformAndDeviceID(memorySizeInput);
	cl_platform_id platform = std::get<0>(dev_tuple);
	cl_device_id device = std::get<1>(dev_tuple);

	props[1] = (cl_context_properties)platform;
	ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
	queue = clCreateCommandQueue(ctx, device, 0, &err);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not setup OpenCL, error code: " << err << std::endl;

	// std::cout << "MMR::Input xy size: " << inputSize << std::endl;

	/* Prepare OpenCL memory objects and place data inside them. */
	deviceBuffer = clCreateBuffer(ctx, CL_MEM_READ_ONLY, memorySizeInput, NULL, &err);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not create OpenCL mem object, error code: " << err << std::endl;

	err = clEnqueueWriteBuffer(queue, deviceBuffer, CL_TRUE, 0, memorySizeInput, buffer, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not create OpenCL buffers, error code: " << err << std::endl;


	cl_uint divider = 128;
	while (true) {
		if (inputSize % divider == 0)
			break;
		divider /= 2;
		if (divider == 2)
			break;
	}


	const cl_uint outputDim = (cl_uint)inputSize / divider;
	size_t local_work_size = 128;
	while (local_work_size > 32) {
		if (((outputDim + (outputDim % local_work_size)) / local_work_size) % 2 != 0) // global must be evenly divisable with local work group
			local_work_size /= 2;
		else
			break;
	}
	const size_t global_work_size = ((outputDim + (outputDim % local_work_size)) / local_work_size) * local_work_size; // seems redundant, but isn't
	//std::cout << "global/local = " << (float)global_work_size / (float)local_work_size << std::endl;
	const size_t memorySizeSub = global_work_size * sizeof(cl_float2); // to avoid access violation in kernel
	//std::cout << "MMM::output xy size: " << outputDim << std::endl;
	cl_float2* sub_buffer = new cl_float2[global_work_size];
	cl_mem  deviceSubBuffer;
	deviceSubBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, memorySizeSub, NULL, &err);
	err = clEnqueueWriteBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);

	if (err != CL_SUCCESS)
		std::cout << "MMR::Could create OpenCL buffers, error code: " << err << std::endl;

	// Create program
	//std::cout << "build min max kernel 2..." << std::endl;
	CreateAndBuildOpenCLProgramFromSourceFile("fdk_opencl.cl", ctx, m_Program);

	m_Kernel = clCreateKernel(m_Program, "min_max_kernel2", &err);

	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not create OpenCL kernel, error code: " << err << std::endl;


	//std::cout << "Kernel created." << std::endl;

	err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&deviceBuffer);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not parse arg 1 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&deviceSubBuffer);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not parse arg 2 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 2, sizeof(cl_uint), &divider);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not parse arg 3 to OpenCL kernel, error code: " << err << std::endl;
	err = clSetKernelArg(m_Kernel, 3, sizeof(cl_uint), &outputDim);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not parse arg 4 to OpenCL kernel, error code: " << err << std::endl;

	cl_event events[2];
	// Execute kernel
	err = clEnqueueNDRangeKernel(queue, m_Kernel, 1, NULL,
		&global_work_size, &local_work_size,
		0, NULL, &events[0]);

	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not run OpenCL kernel, error code: " << err << std::endl;

	err = clWaitForEvents(1, &events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not wait for OpenCL events!?, error code: " << err << std::endl;
	err = clReleaseEvent(events[0]);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not release OpenCL event, error code: " << err << std::endl;

	/* Wait for calculations to be finished. */
	err = clReleaseKernel(m_Kernel);
	err |= clReleaseProgram(m_Program);
	err |= clFinish(queue);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not finish OpenCL queue, error code: " << err << std::endl;

	/* Fetch results of calculations. */
	err = clEnqueueReadBuffer(queue, deviceSubBuffer, CL_TRUE, 0, memorySizeSub, sub_buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not read OpenCL buffer, error code: " << err << std::endl;

	/* Release OpenCL working objects. */
	clReleaseMemObject(deviceSubBuffer);
	clReleaseMemObject(deviceBuffer);
	delete[] buffer;

	clReleaseCommandQueue(queue);
	clReleaseContext(ctx);

	if (err != CL_SUCCESS)
		std::cout << "MMR::Could not create OpenCL kernel, error code: " << err << std::endl;

	return OpenCL_min_max_recurse(sub_buffer, outputDim);
}