/*Utility functions for cbctrecon*/
#include "cbctrecon.h"

// std
#include <algorithm>                              // for max
#include <iostream>                               // for operator<<, endl

// Qt
#include <qstring.h>                              // for QString
#include <qfileinfo.h>
#include <qdir.h>

// ITK
#include "itkImage.h"                             // for Image<>::Pointer
#include "itkImageRegion.h"                       // for operator<<, Image...
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageConstIteratorWithIndex.h"       // for ImageConstIterato...
#include "itkSmartPointer.h"                      // for SmartPointer

// RTK
#include "rtkFieldOfViewImageFilter.h"
#include "rtkThreeDCircularProjectionGeometry.h"  // for ThreeDCircularProje...

template <typename ImageType>
double
GetFOVRadius(const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
             const typename ImageType::Pointer &ProjStack) {

  using FOVfilterType = rtk::FieldOfViewImageFilter<ImageType, ImageType>;
  typename FOVfilterType::Pointer FOVfilter = FOVfilterType::New();
  FOVfilter->SetGeometry(geometry);
  FOVfilter->SetProjectionsStack(ProjStack.GetPointer());
  double x, z;
  double r_inf = -1.0;
  bool hasOverlap = FOVfilter->ComputeFOVRadius(
      FOVfilterType::FOVRadiusType::RADIUSINF, x, z, r_inf);
  // halffan gives r(BOTH)~25, r(SUP)~25, r(INF)~232 -> RADIUSINF also seems to
  // work for fullfan, so we'll use that.

  if (hasOverlap) {
    std::cout << "FOV (inf) radius was found: r=" << r_inf << std::endl;
  }
  double r_sup = -1.0;
  hasOverlap = FOVfilter->ComputeFOVRadius(
      FOVfilterType::FOVRadiusType::RADIUSSUP, x, z, r_sup);
  if (hasOverlap) {
    std::cout << "FOV (sup) radius was found: r=" << r_sup << std::endl;
  }
  return std::max(r_inf, r_sup);
}

template <typename T, typename ImageType>
bool GetOutputResolutionFromFOV(
    typename T::SizeType &sizeOutput, typename T::SpacingType &spacing,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
    const typename ImageType::Pointer &ProjStack,
    const QString &outputFilePath) {

  QFileInfo outFileInfo(outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    double radius = GetFOVRadius<ImageType>(geometry, ProjStack);
    if (radius != -1.0) {
      sizeOutput[0] = 512; // AP
      sizeOutput[1] = 200; // SI
      sizeOutput[2] = 512; // LR
      spacing[0] = 2.0 * radius / sizeOutput[0];
      spacing[1] = 1.0;
      spacing[2] = 2.0 * radius / sizeOutput[2];
      return true;
    }
  }

  return false;
}

double
CbctRecon::GetValueFrom3DImageFloat(int reqX, int reqY, int reqZ,
                                    FloatImageType::Pointer &sp3DFloatImage) {
  if (sp3DFloatImage == nullptr) {
    return -1.0;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it(
      sp3DFloatImage, sp3DFloatImage->GetBufferedRegion());

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int idxX, idxY, idxZ;
  idxZ = 0;

  while (!it.IsAtEnd()) {
    if (idxZ == reqZ) {
      idxY = 0;
      while (!it.IsAtEndOfSlice()) {
        if (idxY == reqY) {
          idxX = 0;
          while (!it.IsAtEndOfLine()) {
            if (idxX == reqX) {
              double tmpVal = it.Get();
              return tmpVal;
            }
            ++it;
            idxX++;
          }
          break;
        }
        it.NextLine();
        idxY++;
      }
      break;
    }
    it.NextSlice();
    idxZ++;
  }

  return -2.0;
}

double CbctRecon::GetValueFrom3DImageUshort(
    int reqX, int reqY, int reqZ, UShortImageType::Pointer &sp3DUshortImage) {
  if (sp3DUshortImage == nullptr) {
    return 0;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      sp3DUshortImage, sp3DUshortImage->GetBufferedRegion());

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int idxX, idxY, idxZ;
  idxZ = 0;

  while (!it.IsAtEnd()) {
    if (idxZ == reqZ) {
      idxY = 0;
      while (!it.IsAtEndOfSlice()) {
        if (idxY == reqY) {

          idxX = 0;
          while (!it.IsAtEndOfLine()) {
            if (idxX == reqX) {
              unsigned short tmpVal = it.Get();
              return tmpVal;
            }
            ++it;
            idxX++;
          }
          break;
        }
        it.NextLine();
        idxY++;
      }
      break;
    }
    it.NextSlice();
    idxZ++;
  }
  return 65535;
}
