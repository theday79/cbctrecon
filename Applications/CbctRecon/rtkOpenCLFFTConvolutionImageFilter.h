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

#ifndef rtkOpenCLFFTConvolutionImageFilter_h
#define rtkOpenCLFFTConvolutionImageFilter_h

#include <itkImage.h>
#include <itkOpenCLImageToImageFilter.h>

namespace rtk
{

	/** \class OpenCLFFTConvolutionImageFilter
	* \brief Implements 1D or 2D FFT convolution.
	*
	* This filter implements a convolution using FFT of the input image. The
	* convolution kernel must be defined in the parent class, passed via the
	* template argument. The template argument must be a child of
	* rtk::FFTConvolutionImageFilter.
	*
	* \see rtk::FFTConvolutionImageFilter
	*
	* \test rtkrampfiltertest.cxx, rtkrampfiltertest2.cxx
	*
	* \author Simon Rit
	*
	* \ingroup CudaImageToImageFilter
	*/

	template< class TParentImageFilter >
	class ITK_EXPORT OpenCLFFTConvolutionImageFilter :
		public itk::OpenCLImageToImageFilter< itk::Image<float, 3U>,
		                                    itk::Image<float, 3U>,
		                                    TParentImageFilter >
		// Well this is just a way of inheritance, not that cuda specific, but convenient for GPU stuff.
	{
	public:
		/** Standard class typedefs. */
		typedef OpenCLFFTConvolutionImageFilter                    Self;
		typedef TParentImageFilter                                 Superclass;
		typedef itk::SmartPointer<Self>                            Pointer;
		typedef itk::SmartPointer<const Self>                      ConstPointer;

		/** Some convenient typedefs. */
		typedef typename TParentImageFilter::RegionType           RegionType;
		typedef typename TParentImageFilter::FFTInputImagePointer FFTInputImagePointer;
		typedef itk::Image<float, 3U>                             FloatImageType;
		
		typedef itk::Image<std::complex<float>, 3U>               OpenCLFFTOutputImageType;
		typedef typename OpenCLFFTOutputImageType::Pointer        OpenCLFFTOutputImagePointer;

		/** Runtime information support. */
		itkTypeMacro(OpenCLFFTConvolutionImageFilter, TParentImageFilter);

	protected:
		OpenCLFFTConvolutionImageFilter();
		~OpenCLFFTConvolutionImageFilter(){}

		virtual void GPUGenerateData();

		/** Pad the inputRegion region of the input image and returns a pointer to the new padded image.
		  * Padding includes a correction for truncation [Ohnesorge, Med Phys, 2000].
		  * centralRegion is the region of the returned image which corresponds to inputRegion.
		  */
		virtual FFTInputImagePointer PadInputImageRegion(const RegionType &inputRegion);
		

	private:
		OpenCLFFTConvolutionImageFilter(const Self&); //purposely not implemented
		void operator=(const Self&);            //purposely not implemented

		OpenCLFFTOutputImagePointer m_KernelFFTOpenCL;
	}; // end of class

} // end namespace rtk

#ifndef ITK_MANUAL_INSTANTIATION
#include "rtkOpenCLFFTConvolutionImageFilter.hxx"
#endif

#endif
