/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef __itkOpenCLImageToImageFilter_h
#define __itkOpenCLImageToImageFilter_h

#include "itkImageToImageFilter.h"
//#include "itkOpenCLKernelManager.h"

namespace itk {

/** \class OpenCLImageToImageFilter
 *
 * \brief class to abstract the behaviour of the OpenCL filters.
 *
 * OpenCLImageToImageFilter is the OpenCL version of ImageToImageFilter.
 * This class can accept both CPU and GPU image as input and output,
 * and apply filter accordingly. If OpenCL is available for use, then
 * GPUGenerateData() is called. Otherwise, GenerateData() in the
 * parent class (i.e., ImageToImageFilter) will be called.
 *
 * \ingroup ITKOpenCLCommon
 */
template <class TInputImage, class TOutputImage,
          class TParentImageFilter =
              ImageToImageFilter<TInputImage, TOutputImage>>
class ITK_EXPORT OpenCLImageToImageFilter : public TParentImageFilter {
public:
  /** Standard class typedefs. */
  typedef OpenCLImageToImageFilter Self;
  typedef TParentImageFilter Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OpenCLImageToImageFilter, TParentImageFilter);

  /** Superclass typedefs. */
  // typedef typename Superclass::DataObjectIdentifierType
  // DataObjectIdentifierType;
  typedef unsigned int DataObjectIdentifierType;

  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType OutputImagePixelType;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename InputImageType::PixelType InputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  // macro to set if OpenCL is used
  itkSetMacro(GPUEnabled, bool);
  itkGetConstMacro(GPUEnabled, bool);
  itkBooleanMacro(GPUEnabled);

  void GenerateData();

protected:
  OpenCLImageToImageFilter();
  ~OpenCLImageToImageFilter();

  virtual void PrintSelf(std::ostream &os, Indent indent) const;

  virtual void GPUGenerateData() {}

  // OpenCL kernel manager
  // typename OpenCLKernelManager::Pointer m_OpenCLKernelManager;

private:
  OpenCLImageToImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);           // purposely not implemented

  bool m_GPUEnabled;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOpenCLImageToImageFilter.hxx"
#endif

#endif
