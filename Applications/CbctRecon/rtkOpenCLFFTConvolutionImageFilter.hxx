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

#ifndef rtkOpenCLFFTConvolutionImageFilter_hxx
#define rtkOpenCLFFTConvolutionImageFilter_hxx

#include "rtkOpenCLFFTConvolutionImageFilter.h"

// Use local RTK FFTW files taken from Gaëtan Lehmann's code for
// thread safety: http://hdl.handle.net/10380/3154
// #include <itkRealToHalfHermitianForwardFFTImageFilter.h> // to-be-replaced
// with clFFT #include <itkHalfHermitianToRealInverseFFTImageFilter.h>

#include "OpenCLFFTFilter.h"
#include <itkCropImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

namespace rtk {

template <class TParentImageFilter>
OpenCLFFTConvolutionImageFilter<
    TParentImageFilter>::OpenCLFFTConvolutionImageFilter() {
// We use FFTW for the kernel so we need to do the same thing as in the parent
#if defined(USE_FFTWF)
  this->SetGreatestPrimeFactor(13);
#endif
}

template <class TParentImageFilter>
typename OpenCLFFTConvolutionImageFilter<
    TParentImageFilter>::FFTInputImagePointer
OpenCLFFTConvolutionImageFilter<TParentImageFilter>::PadInputImageRegion(
    const RegionType &inputRegion) {
  FloatImageType::RegionType inBuffRegion =
      this->GetInput()->GetBufferedRegion();
  if (inBuffRegion != this->GetInput()->GetRequestedRegion()) {
    itkExceptionMacro(<< "OpenCLFFTConvolutionImageFilter assumes that input "
                         "requested and buffered regions are equal.");
  }

  TParentImageFilter::UpdateTruncationMirrorWeights();
  RegionType paddedRegion =
      TParentImageFilter::GetPaddedImageRegion(inputRegion);

  // Create padded image (spacing and origin do not matter)
  typename FFTInputImageType::Pointer paddedImage = FFTInputImageType::New();
  paddedImage->SetRegions(paddedRegion);
  paddedImage->Allocate();

  cl_uint4 sz, sz_i;
  cl_int4 idx;
  idx.x = inBuffRegion.GetIndex()[0] - paddedRegion.GetIndex()[0];
  idx.y = inBuffRegion.GetIndex()[1] - paddedRegion.GetIndex()[1];
  idx.z = inBuffRegion.GetIndex()[2] - paddedRegion.GetIndex()[2];
  sz.x = paddedRegion.GetSize()[0];
  sz.y = paddedRegion.GetSize()[1];
  sz.z = paddedRegion.GetSize()[2];
  sz_i.x = inBuffRegion.GetSize()[0];
  sz_i.y = inBuffRegion.GetSize()[1];
  sz_i.z = inBuffRegion.GetSize()[2];

  OpenCL_padding(idx, sz, sz_i, this->GetInput()->GetBufferPointer(),
                 paddedImage->GetBufferPointer(),
                 TParentImageFilter::m_TruncationMirrorWeights);

  return paddedImage.GetPointer();
}

template <class TParentImageFilter>
void OpenCLFFTConvolutionImageFilter<TParentImageFilter>::GPUGenerateData() {
  // Pad image region
  FFTInputImagePointer paddedImage =
      PadInputImageRegion(this->GetInput()->GetRequestedRegion());

  cl_int4 inputDimension;
  inputDimension.x = paddedImage->GetBufferedRegion().GetSize()[0];
  inputDimension.y = paddedImage->GetBufferedRegion().GetSize()[1];
  inputDimension.z = paddedImage->GetBufferedRegion().GetSize()[2];

  typename FFTInputImageType::RegionType inputreg;
  inputreg.SetSize(paddedImage->GetBufferedRegion().GetSize());
  inputreg.SetIndex(paddedImage->GetBufferedRegion().GetIndex());

  if (inputDimension.y == 1 &&
      inputDimension.z > 1) // Troubles cuda 3.2 and 4.0
    std::swap(inputDimension.y, inputDimension.z);

  typename OpenCLFFTOutputImageType::SizeType s =
      paddedImage->GetLargestPossibleRegion().GetSize();
  this->UpdateFFTConvolutionKernel(s);
  if (this->m_KernelFFTOpenCL.GetPointer() == ITK_NULLPTR ||
      this->m_KernelFFTOpenCL->GetTimeStamp() <
          this->m_KernelFFT->GetTimeStamp()) {

    // Create the Image holding the kernel
    typename OpenCLFFTOutputImageType::RegionType kreg =
        this->m_KernelFFT->GetLargestPossibleRegion();

    this->m_KernelFFTOpenCL = OpenCLFFTOutputImageType::New();
    this->m_KernelFFTOpenCL->SetRegions(kreg);
    this->m_KernelFFTOpenCL->Allocate();

    // clFFT scales by the number of element, correct for it in kernel.
    itk::ImageRegionIterator<typename TParentImageFilter::FFTOutputImageType>
        itKI(this->m_KernelFFT, kreg);
    itk::ImageRegionIterator<OpenCLFFTOutputImageType> itKO(
        this->m_KernelFFTOpenCL, kreg);
    FFTInputPixelType invNPixels;
    invNPixels =
        1.0 / static_cast<double>(
                  paddedImage->GetBufferedRegion().GetNumberOfPixels());
    while (!itKO.IsAtEnd()) {
      itKO.Set(itKI.Get() * invNPixels);
      ++itKI;
      ++itKO;
    }
  }

  // FloatImageType *PadImgP = paddedImage.GetPointer();

  cl_int2 kernelDimension;
  kernelDimension.x = this->m_KernelFFT->GetBufferedRegion().GetSize()[0];
  kernelDimension.y = this->m_KernelFFT->GetBufferedRegion().GetSize()[1];
  std::complex<float> *p_KernelFFTOpenCL =
      this->m_KernelFFTOpenCL
          ->GetBufferPointer(); // reinterpret_cast<std::complex<float>*>(

  // std::cout << "Before Conv.:: Some proj value: " <<
  // paddedImage->GetBufferPointer()[1000] << std::endl;  for (int idz = 0; idz
  // < inputDimension.z; idz++)
  //{
  itk::Index<3U> cur_idx = {0, 0, 0};
  OpenCL_fft_convolution(inputDimension, kernelDimension,
                         &paddedImage->GetPixel(cur_idx), p_KernelFFTOpenCL);
  //}
  std::cout << std::endl;
  // std::cout << "After  Conv.:: Some proj value: " <<
  // paddedImage->GetBufferPointer()[1000] << std::endl;

  // CUDA Cropping and Graft Output
  typedef itk::CropImageFilter<FFTInputImageType, itk::Image<float, 3U>>
      CropFilter;
  typename CropFilter::Pointer cf = CropFilter::New();
  typename Superclass::OutputImageType::SizeType upCropSize, lowCropSize;
  for (unsigned int i = 0; i < 3; i++) {
    lowCropSize[i] = this->GetOutput()->GetRequestedRegion().GetIndex()[i] -
                     paddedImage->GetLargestPossibleRegion().GetIndex()[i];
    upCropSize[i] = paddedImage->GetLargestPossibleRegion().GetSize()[i] -
                    this->GetOutput()->GetRequestedRegion().GetSize()[i] -
                    lowCropSize[i];
  }
  cf->SetUpperBoundaryCropSize(upCropSize);
  cf->SetLowerBoundaryCropSize(lowCropSize);
  cf->SetInput(paddedImage);
  cf->Update();

  // We only want to graft the data. To do so, we copy the rest before grafting.
  cf->GetOutput()->CopyInformation(this->GetOutput());
  cf->GetOutput()->SetBufferedRegion(this->GetOutput()->GetBufferedRegion());
  cf->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
  // Crop and paste result
  itk::ImageRegionConstIterator<itk::Image<float, 3U>> itS(
      cf->GetOutput(), this->GetOutput()->GetBufferedRegion());
  itk::ImageRegionIterator<itk::Image<float, 3U>> itD(
      this->GetOutput(), this->GetOutput()->GetBufferedRegion());
  itS.GoToBegin();
  itD.GoToBegin();
  while (!itS.IsAtEnd()) {
    itD.Set(itS.Get());
    ++itS;
    ++itD;
  }
}
} // end namespace rtk
#endif
