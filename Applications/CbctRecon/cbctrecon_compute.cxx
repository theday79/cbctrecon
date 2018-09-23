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

#include "cbctrecon_compute.h"

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
