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

#ifndef rtkOpenCLForwardProjectionImageFilter_h
#define rtkOpenCLForwardProjectionImageFilter_h

#include "RTKExport.h"
#include "cbctrecon_config.h"

// Conditional definition of the class to pass ITKHeaderTest
#ifdef RTK_USE_OPENCL

#include "rtkForwardProjectionImageFilter.h"

/** \class OpenCLForwardProjectionImageFilter
 * \brief Trilinear interpolation forward projection implemented in OpenCL
 *
 * OpenCLForwardProjectionImageFilter is similar to
 * JosephForwardProjectionImageFilter, except it uses a
 * fixed step between sampling points instead of placing these
 * sampling points only on the main direction slices.
 *
 * The code was developed based on the file tt_project_ray_gpu_kernels.cu of
 * NiftyRec (http://sourceforge.net/projects/niftyrec/) which is distributed
 * under a BSD license. See COPYRIGHT.TXT.
 *
 * \author Marc Vila, updated by Simon Rit and Cyril Mory
 *
 * \ingroup RTK Projector OpenCLImageToImageFilter
 */

namespace rtk {

template <class TInputImage = itk::Image<float, 3>,
          class TOutputImage = itk::Image<float, 3>>
class CBCTRECON_API OpenCLForwardProjectionImageFilter
    : public ForwardProjectionImageFilter<TInputImage, TOutputImage> {
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(OpenCLForwardProjectionImageFilter);

  /** Standard class type alias. */
  using Self = OpenCLForwardProjectionImageFilter;
  using Superclass = ForwardProjectionImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  using VectorType = itk::Vector<float, 3>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OpenCLForwardProjectionImageFilter,
               ForwardProjectionImageFilter);

  /** Set step size along ray (in mm). Default is 1 mm. */
  itkGetConstMacro(StepSize, double);
  itkSetMacro(StepSize, double);

protected:
  OpenCLForwardProjectionImageFilter();
  ~OpenCLForwardProjectionImageFilter(){};

  void GenerateData() override;

private:
  double m_StepSize;
}; // end of class

} // end namespace rtk

#ifndef ITK_MANUAL_INSTANTIATION
#include "rtkOpenCLForwardProjectionImageFilter.hxx"
#endif

#endif // end conditional definition of the class

#endif
