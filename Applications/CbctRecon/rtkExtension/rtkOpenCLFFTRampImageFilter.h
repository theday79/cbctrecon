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

#ifndef rtkOpenCLFFTRampImageFilter_h
#define rtkOpenCLFFTRampImageFilter_h

#include "rtkConfiguration.h"
#include "rtkFFTRampImageFilter.h"
#include "rtkMacro.h"
#include "rtkOpenCLFFTConvolutionImageFilter.h"
#include <itkConceptChecking.h>

namespace rtk {

/** \class FFTRampImageFilter
 * \brief Implements the ramp image filter of the filtered backprojection
 * algorithm.
 *
 * The filter code is based on FFTConvolutionImageFilter by Gaetan Lehmann
 * (see http://hdl.handle.net/10380/3154)
 *
 * \test rtkrampfiltertest.cxx
 *
 * \author Simon Rit
 *
 * \ingroup ImageToImageFilter
 */

class OpenCLFFTRampImageFilter
    : public OpenCLFFTConvolutionImageFilter<FFTRampImageFilter<
          itk::Image<float, 3>, itk::Image<float, 3>, float>> {
public:
  /** Standard class typedefs. */
  typedef OpenCLFFTRampImageFilter Self;
  typedef FFTRampImageFilter<FloatImageType, FloatImageType, float> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(OpenCLFFTRampImageFilter, FFTRampImageFilter);

protected:
  rtkcuda_EXPORT OpenCLFFTRampImageFilter() {}
  ~OpenCLFFTRampImageFilter() {}

  // virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  /** Creates and return a pointer to one line of the ramp kernel in Fourier
   * space. Used in generate data functions.  */
  // void UpdateFFTConvolutionKernel(const SizeType size) ITK_OVERRIDE;

private:
  OpenCLFFTRampImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);           // purposely not implemented
};                                        // end of class

} // end namespace rtk

//#ifndef ITK_MANUAL_INSTANTIATION
//#include "rtkOpenCLFFTRampImageFilter.hxx" //do it the cuda way in convolution
// filter only #endif

#endif
