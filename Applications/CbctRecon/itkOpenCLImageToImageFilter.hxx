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
#ifndef __itkOpenCLImageToImageFilter_hxx
#define __itkOpenCLImageToImageFilter_hxx

#include "itkOpenCLImageToImageFilter.h"

namespace itk
{
	/*template <class T>
	class OpenCLTraits
	{
	public:
		typedef T Type;
	};

	template <class TPixelType, unsigned int NDimension>
	class OpenCLTraits< Image< TPixelType, NDimension > >
	{
	public:
		typedef itk::Image<TPixelType, NDimension> Type;
	};*/

	template< class TInputImage, class TOutputImage, class TParentImageFilter >
	OpenCLImageToImageFilter< TInputImage, TOutputImage, TParentImageFilter >::OpenCLImageToImageFilter() : m_GPUEnabled(true)
	{
	//	m_OpenCLKernelManager = OpenCLKernelManager::New();
	}

	template< class TInputImage, class TOutputImage, class TParentImageFilter >
	OpenCLImageToImageFilter< TInputImage, TOutputImage, TParentImageFilter >::~OpenCLImageToImageFilter()
	{
	}

	template< class TInputImage, class TOutputImage, class TParentImageFilter >
	void
		OpenCLImageToImageFilter< TInputImage, TOutputImage, TParentImageFilter >::PrintSelf(std::ostream & os,
			Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "GPU: " << (m_GPUEnabled ? "Enabled" : "Disabled") << std::endl;
	}

	template< class TInputImage, class TOutputImage, class TParentImageFilter >
	void
		OpenCLImageToImageFilter< TInputImage, TOutputImage, TParentImageFilter >::GenerateData()
	{
		if (!m_GPUEnabled) // call CPU update function
		{
			Superclass::GenerateData();
		}
		else // call OpenCL update function
		{
			// Call a method to allocate memory for the filter's outputs
			this->AllocateOutputs();

			GPUGenerateData();
		}
	}
	
	/*template< class TInputImage, class TOutputImage, class TParentImageFilter >
	void
		OpenCLImageToImageFilter< TInputImage, TOutputImage, TParentImageFilter >::GraftOutput(DataObject *output)
	{
		typedef typename itk::OpenCLTraits< TOutputImage >::Type OpenCLOutputImage;
		typename OpenCLOutputImage::Pointer otPtr = dynamic_cast<OpenCLOutputImage *>(this->GetOutput());
		//typename TOutputImage::Pointer otPtr = this->GetOutput();
		//otPtr = output;
		//otPtr->GraftOutput(output);
	}*/
	/*
	template< class TInputImage, class TOutputImage, class TParentImageFilter >
	void
		OpenCLImageToImageFilter< TInputImage, TOutputImage, TParentImageFilter >::GraftOutput(const DataObjectIdentifierType & key, DataObject *output)
	{
		// Until I choose to implement double, which will probably never happen.
		//typedef typename itk::OpenCLTraits< TOutputImage >::Type OpenCLOutputImage;
		//typename OpenCLOutputImage::Pointer otPtr = dynamic_cast<OpenCLOutputImage *>(this->ProcessObject::GetOutput(key));
		Superclass::Graft(output);
		//otPtr->Graft(output);
	}*/
	
} // end namespace itk

#endif
