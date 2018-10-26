#ifndef CBCTRECON_COMPUTE_H
#define CBCTRECON_COMPUTE_H

#include "cbctrecon_config.h"

#include <algorithm> // for std::max
#include <iostream>

#include <qdir.h>
#include <qfileinfo.h>
#include <qstring.h>

#include "itkImage.h"

#include "rtkFieldOfViewImageFilter.h"
#include "rtkThreeDCircularProjectionGeometry.h"

#include "cbctrecon_types.h"

class QXmlStreamReader;

CBCTRECON_API void ApplyBowtie(ProjReaderType::Pointer &reader,
                               FilterReaderType::Pointer &bowtie_reader);

CBCTRECON_API double GetMaxAndMinValueOfProjectionImage(
    double &fProjImgValueMax, double &fProjImgValueMin,
    const FloatImageType::Pointer &projImage); // , double theoreticalMin);

CBCTRECON_API void Get2DFrom3D(UShortImageType::Pointer &spSrcImg3D,
                               FloatImage2DType::Pointer &spTargetImg2D,
                               int idx, enPLANE iDirection);

CBCTRECON_API double
CalculateIntensityScaleFactorFromMeans(UShortImageType::Pointer &spProjRaw3D,
                                       UShortImageType::Pointer &spProjCT3D);
CBCTRECON_API double GetRawIntensityScaleFactor(QString &strRef_mAs,
                                                QString &strCur_mAs);
CBCTRECON_API void TransformationRTK2IEC(FloatImageType::Pointer &spSrcTarg);
CBCTRECON_API QString XML_GetSingleItemString(QXmlStreamReader &xml);

CBCTRECON_API bool GetXrayParamFromINI(QString &strPathINI, float &kVp,
                                       float &mA, float &ms);
CBCTRECON_API void AddConstHU(UShortImageType::Pointer &spImg, int HUval);
// Read long INIXVI text file and read couch shift values. apply cm -> mm
// conversion (multiply 10). NO sign changes.
CBCTRECON_API bool GetCouchShiftFromINIXVI(QString &strPathINIXVI,
                                           VEC3D *pTrans, VEC3D *pRot);
// This function came from the tracking project. trans values are all in mm,
// DICOM x, y, z
CBCTRECON_API void
ImageTransformUsingCouchCorrection(UShortImageType::Pointer &spUshortInput,
                                   UShortImageType::Pointer &spUshortOutput,
                                   VEC3D couch_trans, VEC3D couch_rot);

CBCTRECON_API void RotateImgBeforeFwd(UShortImageType::Pointer &spInputImgUS,
                                      UShortImageType::Pointer &spOutputImgUS);

CBCTRECON_API void
ConvertUshort2AttFloat(UShortImageType::Pointer &spImgUshort,
                       FloatImageType::Pointer &spAttImgFloat);

template <typename RefImageType, typename TargetImageType>
void AllocateByRef(typename RefImageType::Pointer &spRefImg3D,
                   typename TargetImageType::Pointer &spTarImg3D) {
  const auto sizeSrc = spRefImg3D->GetBufferedRegion().GetSize();
  const auto startSrc = spRefImg3D->GetBufferedRegion().GetIndex();

  const auto spacingSrc = spRefImg3D->GetSpacing();
  const auto originSrc = spRefImg3D->GetOrigin();

  typename RefImageType::RegionType region;
  region.SetSize(sizeSrc);
  region.SetIndex(startSrc);

  spTarImg3D = TargetImageType::New();

  spTarImg3D->SetRegions(region);
  spTarImg3D->SetSpacing(spacingSrc);
  spTarImg3D->SetOrigin(originSrc);

  spTarImg3D->Allocate();
  spTarImg3D->FillBuffer(0);
}

template <typename T>
T GetValueFrom3DImage(const unsigned int reqX, const unsigned int reqY,
                      const unsigned int reqZ,
                      typename itk::Image<T, 3>::Pointer &sp3DImage) {
  if (sp3DImage == nullptr) {
    return 0;
  }

  const auto idx = FloatImageType::IndexType{reqX, reqY, reqZ};
  if (!sp3DImage->GetBufferedRegion().IsInside(idx)) {
    return static_cast<T>(-1); // underflow on purpose
  }

  return sp3DImage->GetPixel(idx);
}

template <typename ImageType>
double
GetFOVRadius(const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
             const typename ImageType::Pointer &ProjStack) {

  using FOVfilterType = rtk::FieldOfViewImageFilter<ImageType, ImageType>;
  typename FOVfilterType::Pointer FOVfilter = FOVfilterType::New();
  FOVfilter->SetGeometry(geometry);
  FOVfilter->SetProjectionsStack(ProjStack.GetPointer());
  double x, z;
  auto r_inf = -1.0;
  bool hasOverlap = FOVfilter->ComputeFOVRadius(
      FOVfilterType::FOVRadiusType::RADIUSINF, x, z, r_inf);
  // halffan gives r(BOTH)~25, r(SUP)~25, r(INF)~232 -> RADIUSINF also seems to
  // work for fullfan, so we'll use that.

  if (hasOverlap) {
    std::cout << "FOV (inf) radius was found: r=" << r_inf << std::endl;
  }
  auto r_sup = -1.0;
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
  auto outFileDir = outFileInfo.absoluteDir();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    const double radius = GetFOVRadius<ImageType>(geometry, ProjStack);
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

#endif // CBCTRECON_COMPUTE_H
