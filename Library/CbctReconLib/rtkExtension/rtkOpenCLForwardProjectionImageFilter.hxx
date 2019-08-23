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

#ifndef rtkOpenCLForwardProjectionImageFilter_hxx
#define rtkOpenCLForwardProjectionImageFilter_hxx

#include "rtkConfiguration.h"
// Conditional definition of the class to pass ITKHeaderTest
#ifdef RTK_USE_OPENCL

#include "rtkOpenCLForwardProjectionImageFilter.h"

#include <itkCastImageFilter.h>
#include <itkMacro.h>

#include "OpenCL/ForwardProjectionImageFilter.hpp"

namespace rtk {

template <class TInputImage, class TOutputImage>
OpenCLForwardProjectionImageFilter<
    TInputImage, TOutputImage>::OpenCLForwardProjectionImageFilter()
    : m_StepSize(1) {}

template <class TInputImage, class TOutputImage>
void OpenCLForwardProjectionImageFilter<TInputImage,
                                        TOutputImage>::GenerateData() {
  if (!this->GetGeometry()->GetSourceToDetectorDistances().empty() &&
      this->GetGeometry()->GetSourceToDetectorDistances()[0] == 0) {
    itkGenericExceptionMacro(
        << "Parallel geometry is not handled by OpenCL forward projector.");
  }

  const typename Superclass::GeometryType *geometry = this->GetGeometry();
  const unsigned int Dimension = TInputImage::ImageDimension;
  const size_t iFirstProj =
      this->GetInput(0)->GetRequestedRegion().GetIndex(Dimension - 1);
  const size_t nProj =
      this->GetInput(0)->GetRequestedRegion().GetSize(Dimension - 1);

  OpenCL_forwardProject_options fwd_opts;

  fwd_opts.vectorLength =
      itk::PixelTraits<typename TInputImage::PixelType>::Dimension;

  fwd_opts.projSize.at(0) = this->GetOutput()->GetBufferedRegion().GetSize()[0];
  fwd_opts.projSize.at(1) = this->GetOutput()->GetBufferedRegion().GetSize()[1];

  const auto nPixelsPerProj =
      fwd_opts.projSize.at(0) * fwd_opts.projSize.at(1) * fwd_opts.vectorLength;

  auto largestReg = this->GetOutput()->GetLargestPossibleRegion();
  this->GetOutput()->SetRegions(largestReg);
  this->GetOutput()->Allocate();

  // Setting BoxMin and BoxMax
  // SR: we are using textures (read_imagef sampling) where the pixel definition
  // is not center but corner. Therefore, we set the box limits from index to
  // index+size instead of, for ITK, index-0.5 to index+size-0.5.
  auto &vol_buffer_reg = this->GetInput(1)->GetBufferedRegion();
  for (unsigned int i = 0; i < 3; i++) {
    fwd_opts.box_min.at(i) = vol_buffer_reg.GetIndex()[i] + 0.5;
    fwd_opts.box_max.at(i) =
        fwd_opts.box_min.at(i) + vol_buffer_reg.GetSize()[i] - 1.0;
  }

  // Getting Spacing
  for (unsigned int i = 0; i < 3; i++) {
    fwd_opts.spacing.at(i) = this->GetInput(1)->GetSpacing()[i];
  }
  fwd_opts.t_step = m_StepSize;

  // Getting dimensions
  fwd_opts.volSize.at(0) = vol_buffer_reg.GetSize()[0];
  fwd_opts.volSize.at(1) = vol_buffer_reg.GetSize()[1];
  fwd_opts.volSize.at(2) = vol_buffer_reg.GetSize()[2];

  auto cast_filter_0 =
      itk::CastImageFilter<TInputImage, itk::Image<float, 3>>::New();
  cast_filter_0->SetInput(this->GetInput(0));
  cast_filter_0->Update();
  auto *pin = cast_filter_0->GetOutput()->GetBufferPointer();

  auto cast_filter_1 =
      itk::CastImageFilter<TInputImage, itk::Image<float, 3>>::New();
  cast_filter_1->SetInput(this->GetInput(1));
  cast_filter_1->Update();
  auto *pvol = cast_filter_1->GetOutput()->GetBufferPointer();

  auto *pout = this->GetOutput()->GetBufferPointer();

  // Account for system rotations
  typename Superclass::GeometryType::ThreeDHomogeneousMatrixType volPPToIndex;
  volPPToIndex = GetPhysicalPointToIndexMatrix(this->GetInput(1));

  // Compute matrix to translate the pixel indices on the volume and the
  // detector if the Requested region has non-zero index
  typename Superclass::GeometryType::ThreeDHomogeneousMatrixType
      projIndexTranslation,
      volIndexTranslation;
  projIndexTranslation.SetIdentity();
  volIndexTranslation.SetIdentity();
  for (unsigned int i = 0; i < 3; i++) {
    projIndexTranslation[i][3] =
        this->GetOutput()->GetRequestedRegion().GetIndex(i);
    volIndexTranslation[i][3] = -vol_buffer_reg.GetIndex(i);

    // Adding 0.5 offset to change from the centered pixel convention (ITK)
    // to the corner pixel convention (OpenCL).
    volPPToIndex[i][3] += 0.5;
  }

  // Compute matrices to transform projection index to volume index, one per
  // projection
  auto translatedProjectionIndexTransformMatrices =
      std::vector<float>(nProj * 12);
  auto translatedVolumeTransformMatrices = std::vector<float>(nProj * 12);
  auto source_positions = std::vector<float>(nProj * 3);

  fwd_opts.radiusCylindricalDetector = geometry->GetRadiusCylindricalDetector();

  // Go over each projection
  for (auto iProj = iFirstProj; iProj < iFirstProj + nProj; iProj++) {
    const auto new_proj_index = iProj - iFirstProj;

    // The matrices required depend on the type of detector
    if (fabs(fwd_opts.radiusCylindricalDetector) < 0.001) {
      const auto translatedProjectionIndexTransformMatrix =
          volIndexTranslation.GetVnlMatrix() * volPPToIndex.GetVnlMatrix() *
          geometry->GetProjectionCoordinatesToFixedSystemMatrix(iProj)
              .GetVnlMatrix() *
          rtk::GetIndexToPhysicalPointMatrix(this->GetInput()).GetVnlMatrix() *
          projIndexTranslation.GetVnlMatrix();
      for (auto j = 0; j < 3; j++) { // Ignore the 4th row
        for (auto k = 0; k < 4; k++) {
          translatedProjectionIndexTransformMatrices.at(new_proj_index * 12 +
                                                        j * 4 + k) =
              static_cast<float>(
                  translatedProjectionIndexTransformMatrix[j][k]);
        }
      }
    } else {
      const auto translatedProjectionIndexTransformMatrix =
          geometry->GetProjectionCoordinatesToDetectorSystemMatrix(iProj)
              .GetVnlMatrix() *
          rtk::GetIndexToPhysicalPointMatrix(this->GetInput()).GetVnlMatrix() *
          projIndexTranslation.GetVnlMatrix();
      for (auto j = 0; j < 3; j++) { // Ignore the 4th row
        for (auto k = 0; k < 4; k++) {
          translatedProjectionIndexTransformMatrices.at(new_proj_index * 12 +
                                                        j * 4 + k) =
              static_cast<float>(
                  translatedProjectionIndexTransformMatrix[j][k]);
        }
      }

      const auto translatedVolumeTransformMatrix =
          volIndexTranslation.GetVnlMatrix() * volPPToIndex.GetVnlMatrix() *
          geometry->GetRotationMatrices()[iProj].GetInverse();
      for (auto j = 0; j < 3; j++) { // Ignore the 4th row
        for (auto k = 0; k < 4; k++) {
          translatedVolumeTransformMatrices.at(new_proj_index * 12 + j * 4 +
                                               k) =
              static_cast<float>(translatedVolumeTransformMatrix[j][k]);
        }
      }
    }

    // Compute source position in volume indices
    auto source_position = volPPToIndex * geometry->GetSourcePosition(iProj);

    // Copy it into a single large array
    for (unsigned int d = 0; d < 3; d++) {
      source_positions.at(new_proj_index * 3 + d) =
          source_position[d]; // Ignore the 4th component
    }
  }

  for (size_t i = 0; i < nProj; i += static_cast<size_t>(SLAB_SIZE)) {
    // If nProj is not a multiple of SLAB_SIZE, the last slab will contain less
    // than SLAB_SIZE projections
    const auto cur_slab = std::min(nProj - i, static_cast<size_t>(SLAB_SIZE));
    fwd_opts.projSize[2] = cur_slab;
    const auto projectionOffset =
        iFirstProj + i - this->GetOutput()->GetBufferedRegion().GetIndex(2);

    // Get the geometry corresonding to the relevant projections:
    std::copy(translatedProjectionIndexTransformMatrices.begin() + i * 12,
              translatedProjectionIndexTransformMatrices.begin() +
                  (i + cur_slab) * 12,
              std::back_inserter(
                  fwd_opts.translatedProjectionIndexTransformMatrices));

    std::copy(translatedVolumeTransformMatrices.begin() + i * 12,
              translatedVolumeTransformMatrices.begin() + (i + cur_slab) * 12,
              std::back_inserter(fwd_opts.translatedVolumeTransformMatrices));

    std::copy(source_positions.begin() + i * 3,
              source_positions.begin() + (i + cur_slab) * 3,
              std::back_inserter(fwd_opts.source_positions));

    // Run the forward projection with a slab of SLAB_SIZE or less projections
    OpenCL_forward_project(pin + (nPixelsPerProj * projectionOffset),
                           pout + (nPixelsPerProj * projectionOffset), pvol,
                           fwd_opts);

    fwd_opts.translatedProjectionIndexTransformMatrices.clear();
    fwd_opts.translatedVolumeTransformMatrices.clear();
    fwd_opts.source_positions.clear();
  }
}

} // end namespace rtk

#endif // end conditional definition of the class

#endif
