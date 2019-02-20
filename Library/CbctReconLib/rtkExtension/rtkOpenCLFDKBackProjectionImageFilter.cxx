// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
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

// #define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#include "rtkOpenCLFDKBackProjectionImageFilter.h"

#include "rtkOpenCLUtilities.h"

namespace rtk {

OpenCLFDKBackProjectionImageFilter ::OpenCLFDKBackProjectionImageFilter() =
    default;

void OpenCLFDKBackProjectionImageFilter ::InitDevice() {
  // OpenCL init (platform, device, context and command queue)
  auto platforms = GetListOfOpenCLPlatforms();
  auto devices = GetListOfOpenCLDevices(platforms[0]);
  cl_context_properties properties[] = {
      CL_CONTEXT_PLATFORM,
      reinterpret_cast<cl_context_properties>(platforms[0]), 0};

  cl_int error;
  m_Context = clCreateContext(properties, devices.size(), &devices[0], nullptr,
                              nullptr, &error);
  if (error != CL_SUCCESS)
    itkExceptionMacro(<< "Could not create OpenCL context, error code: "
                      << error);

  m_CommandQueue = clCreateCommandQueue(m_Context, devices[0],
                                        CL_QUEUE_PROFILING_ENABLE, &error);
  if (error != CL_SUCCESS)
    itkExceptionMacro(<< "Could not create OpenCL command queue, error code: "
                      << error);

  // OpenCL memory allocation
  m_DeviceMatrix = clCreateBuffer(m_Context, CL_MEM_READ_ONLY,
                                  sizeof(cl_float) * 12, nullptr, &error);
  if (error != CL_SUCCESS)
    itkExceptionMacro(<< "Could not allocate OpenCL matrix buffer, error code: "
                      << error);

  auto volBytes = this->GetOutput()->GetRequestedRegion().GetNumberOfPixels() *
                  sizeof(float);
  m_DeviceVolume = clCreateBuffer(
      m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, volBytes,
      (void *)this->GetInput()->GetBufferPointer(), &error);
  if (error != CL_SUCCESS)
    itkExceptionMacro(<< "Could not allocate OpenCL volume buffer, error code: "
                      << error);

  cl_image_format projFormat{};
  projFormat.image_channel_order = CL_INTENSITY;
  projFormat.image_channel_data_type = CL_FLOAT;
  cl_image_desc im_desc{};
  im_desc.image_type = CL_MEM_OBJECT_IMAGE2D;
  im_desc.image_width =
      this->GetInput(1)->GetLargestPossibleRegion().GetSize()[0];
  im_desc.image_height =
      this->GetInput(1)->GetLargestPossibleRegion().GetSize()[1];
  im_desc.image_row_pitch = 0;
  im_desc.num_mip_levels = 0;
  im_desc.num_samples = 0;

  im_desc.image_depth = 0;
  im_desc.image_array_size = 1; // If you ever change the above from Image2D to
                                // Image2D Array performance will be worsened.
  im_desc.image_slice_pitch = 0;
#ifdef CL_VERSION_2_0
  im_desc.mem_object = NULL;
#else
  im_desc.buffer = nullptr;
#endif

  m_DeviceProjection = clCreateImage(m_Context, CL_MEM_READ_ONLY, &projFormat,
                                     &im_desc, nullptr, &error);

  if (error != CL_SUCCESS)
    itkExceptionMacro(
        << "Could not allocate OpenCL projection image, error code: " << error);

  std::string cl_file("fdk_opencl.cl");
  CreateAndBuildOpenCLProgramFromSourceFile(cl_file, m_Context, m_Program);

  m_Kernel = clCreateKernel(m_Program,
                            "OpenCLFDKBackProjectionImageFilterKernel", &error);
  if (error != CL_SUCCESS)
    itkExceptionMacro(<< "Could not create OpenCL kernel, error code: "
                      << error);

  // Set kernel parameters
  const auto img_size = this->GetOutput()->GetRequestedRegion().GetSize();
  cl_uint4 volumeDim{};
  volumeDim.s[0] = img_size[0];
  volumeDim.s[1] = img_size[1];
  volumeDim.s[2] = img_size[2];
  volumeDim.s[3] = 1;
  OPENCL_CHECK_ERROR(
      clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), &m_DeviceVolume));
  OPENCL_CHECK_ERROR(
      clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), &m_DeviceMatrix));
  OPENCL_CHECK_ERROR(
      clSetKernelArg(m_Kernel, 2, sizeof(cl_mem), &m_DeviceProjection));
  OPENCL_CHECK_ERROR(clSetKernelArg(m_Kernel, 3, sizeof(cl_uint4), &volumeDim));
}

void OpenCLFDKBackProjectionImageFilter ::CleanUpDevice() {
  auto volBytes = this->GetOutput()->GetRequestedRegion().GetNumberOfPixels() *
                  sizeof(float);

  OPENCL_CHECK_ERROR(clReleaseKernel(m_Kernel));
  OPENCL_CHECK_ERROR(clReleaseProgram(m_Program));
  OPENCL_CHECK_ERROR(clFinish(m_CommandQueue));
  OPENCL_CHECK_ERROR(clEnqueueReadBuffer(
      m_CommandQueue, m_DeviceVolume, CL_TRUE, 0, volBytes,
      this->GetOutput()->GetBufferPointer(), 0, nullptr, nullptr));
  OPENCL_CHECK_ERROR(clReleaseMemObject(m_DeviceProjection));
  OPENCL_CHECK_ERROR(clReleaseMemObject(m_DeviceVolume));
  OPENCL_CHECK_ERROR(clReleaseMemObject(m_DeviceMatrix));
  OPENCL_CHECK_ERROR(clReleaseCommandQueue(m_CommandQueue));
  OPENCL_CHECK_ERROR(clReleaseContext(m_Context));
}

void OpenCLFDKBackProjectionImageFilter ::GenerateData() {
  this->AllocateOutputs();

  const auto Dimension = ImageType::ImageDimension;
  const auto input1_largestregion =
      this->GetInput(1)->GetLargestPossibleRegion();
  const unsigned int nProj =
      input1_largestregion.GetSize(Dimension - 1);
  const unsigned int iFirstProj =
      input1_largestregion.GetIndex(Dimension - 1);

  // Ramp factor is the correction for ramp filter which did not account for the
  // divergence of the beam
  // const GeometryPointer geometry =
  //    dynamic_cast<GeometryType *>(this->GetGeometry().GetPointer());

  // Rotation center (assumed to be at 0 yet)
  ImageType::PointType rotCenterPoint;
  rotCenterPoint.Fill(0.0);
  itk::ContinuousIndex<double, Dimension> rotCenterIndex;
  if (!this->GetInput(0)->TransformPhysicalPointToContinuousIndex(rotCenterPoint,
                                                             rotCenterIndex))
  {
    std::cerr << "Center was outside of cube!??" << std::endl;
    return;
  }

  // Include non-zero index in matrix
  itk::Matrix<double, 4, 4> matrixIdxVol;
  matrixIdxVol.SetIdentity();
  for (unsigned int i = 0; i < 3; i++) {
    matrixIdxVol[i][3] = this->GetOutput()->GetRequestedRegion().GetIndex()[i];
    rotCenterIndex[i] -= matrixIdxVol[i][3];
  }
  // Go over each projection
  for (auto iProj = iFirstProj; iProj < iFirstProj + nProj; iProj++) {
    // Extract the current slice
    auto projection = this->GetProjection<ProjectionImageType>(iProj);

    // Index to index matrix normalized to have a correct backprojection weight
    // (1 at the isocenter)
    auto matrix = GetIndexToIndexProjectionMatrix(iProj);

    // We correct the matrix for non zero indexes
    itk::Matrix<double, 3, 3> matrixIdxProj;
    matrixIdxProj.SetIdentity();
    for (unsigned int i = 0; i < 2; i++) { // SR: 0.5 for 2D texture
      matrixIdxProj[i][2] =
          -1 * projection->GetBufferedRegion().GetIndex()[i] + 0.5;
    }

    matrix = matrixIdxProj.GetVnlMatrix() * matrix.GetVnlMatrix() *
             matrixIdxVol.GetVnlMatrix();

    auto perspFactor = matrix[Dimension - 1][Dimension];
    for (unsigned int j = 0; j < Dimension; j++) {
      perspFactor += matrix[Dimension - 1][j] * rotCenterIndex[j];
    }

    matrix /= perspFactor;

    float fMatrix[12];
    for (auto j = 0; j < 3;
         j++) { // Manual unroll 3*(one less +, 4 less %, 4 / is now * )
      fMatrix[j * 4] = matrix[j][0];
      fMatrix[j * 4 + 1] = matrix[j][1];
      fMatrix[j * 4 + 2] = matrix[j][2];
      fMatrix[j * 4 + 3] = matrix[j][3];
    }

    OPENCL_CHECK_ERROR(clEnqueueWriteBuffer(
        m_CommandQueue, m_DeviceMatrix, CL_TRUE, 0, 12 * sizeof(float),
        (void *)&fMatrix[0], 0, nullptr, nullptr));

    const size_t origin[] = {0, 0, 0};
    const size_t region[] = {
        this->GetInput(1)->GetRequestedRegion().GetSize()[0],
        this->GetInput(1)->GetRequestedRegion().GetSize()[1], 1};
    OPENCL_CHECK_ERROR(clEnqueueWriteImage(
        m_CommandQueue, m_DeviceProjection, CL_TRUE, &origin[0], &region[0], 0,
        0, (void *)&projection->GetBufferPointer()[0], 0, nullptr, nullptr));

    // Execute kernel
    cl_event events[2];
    size_t local_work_size = 128;
    auto global_work_size =
        this->GetOutput()->GetRequestedRegion().GetNumberOfPixels();
    OPENCL_CHECK_ERROR(clEnqueueNDRangeKernel(
        m_CommandQueue, m_Kernel, 1, nullptr, &global_work_size,
        &local_work_size, 0, nullptr, &events[0]));
    OPENCL_CHECK_ERROR(clWaitForEvents(1, &events[0]));
    OPENCL_CHECK_ERROR(clReleaseEvent(events[0]));
  }
}

} // end namespace rtk
