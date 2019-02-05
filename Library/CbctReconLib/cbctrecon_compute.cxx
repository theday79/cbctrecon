// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

/*Utility functions for cbctrecon*/
#include "cbctrecon_compute.h"
#include "cbctrecon.h"

#ifdef USE_OPENMP
#include <omp.h> // Resharper lies, omp is used
#endif

// Qt
#include <qxmlstream.h>

// ITK
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkImage.h"       // for Image<>::Pointer
#include "itkImageRegion.h" // for operator<<, Image...
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSmartPointer.h" // for SmartPointer
#include "itkStatisticsImageFilter.h"

#include "OpenCLFFTFilter.h"

void ApplyBowtie(ProjReaderType::Pointer &reader,
                 FilterReaderType::Pointer &bowtie_reader) {

  OpenCL_subtract3Dfrom2DbySlice_InPlace(
      static_cast<cl_float *>(reader->GetOutput()->GetBufferPointer()),
      static_cast<cl_float *>(bowtie_reader->GetOutput()->GetBufferPointer()),
      reader->GetOutput()->GetLargestPossibleRegion().GetSize(),
      bowtie_reader->GetOutput()->GetLargestPossibleRegion().GetSize());
}

double GetMaxAndMinValueOfProjectionImage(
    double &fProjImgValueMax, double &fProjImgValueMin,
    const FloatImageType::Pointer &projImage) //, const double theoreticalMin)
{
  if (projImage == nullptr) {
    fProjImgValueMax = -1.0;
    fProjImgValueMin = -1.0;
    return 20.0;
  }

  using MinMaxCalcType = itk::MinimumMaximumImageCalculator<FloatImageType>;
  auto MinMaxFilter = MinMaxCalcType::New();
  MinMaxFilter->SetImage(projImage);
  MinMaxFilter->Compute();
  fProjImgValueMin = MinMaxFilter->GetMinimum();
  fProjImgValueMax = MinMaxFilter->GetMaximum();

  std::cout << "Min: " << fProjImgValueMin << " max: " << fProjImgValueMax
            << std::endl;

  return fProjImgValueMin;
}

/*
ui.lineEdit_CurmAs->setText(QString("%1,20").arg((64 * 40 / 20) / ScaleFactor));
ui.lineEdit_RefmAs->setText(QString("64,40"));
*/
double
CalculateIntensityScaleFactorFromMeans(UShortImageType::Pointer &spProjRaw3D,
                                       UShortImageType::Pointer &spProjCT3D) {
  using StatFilterType = itk::StatisticsImageFilter<UShortImageType>;
  auto raw_mean = 0.0;
  auto ctMean = 0.0;
#pragma omp parallel sections
  {
#pragma omp section
    {
      auto statFilter = StatFilterType::New();
      statFilter->SetInput(spProjRaw3D);
      statFilter->Update();

      raw_mean = statFilter->GetMean();
    }
#pragma omp section
    {
      auto statFilter = StatFilterType::New();
      statFilter->SetInput(spProjCT3D);
      statFilter->Update();

      ctMean = statFilter->GetMean();
    }
  }
  return ctMean / raw_mean;
}

// spSrcImg3D: usually projImage in USHORT type
void Get2DFrom3D(UShortImageType::Pointer &spSrcImg3D,
                 FloatImage2DType::Pointer &spTargetImg2D, const int idx,
                 const enPLANE iDirection) {
  if (spSrcImg3D == nullptr) {
    return;
  }

  auto idx_hor = 0, idxVer = 0, idxZ = 0;

  switch (iDirection) {
  case PLANE_AXIAL:
    idx_hor = 0;
    idxVer = 1;
    idxZ = 2;
    break;
  case PLANE_FRONTAL:
    idx_hor = 0;
    idxVer = 2;
    idxZ = 1;
    break;
  case PLANE_SAGITTAL:
    idx_hor = 1;
    idxVer = 2;
    idxZ = 0;
    break;
  }

  // Create 2D target image based on geometry of 3D
  auto imgDim = spSrcImg3D->GetBufferedRegion().GetSize();
  auto spacing = spSrcImg3D->GetSpacing();
  // UShortImageType::PointType origin = spSrcImg3D->GetOrigin();

  // int width = imgDim[idxHor];
  // int height = imgDim[idxVer];
  const int z_size = imgDim[idxZ];
  // std::cout << "Get2DFrom3D zSize = " << zSize << std::endl;

  if (idx < 0 || idx >= z_size) {
    std::cout << "Error! idx is out of the range" << std::endl;
    return;
  }

  FloatImage2DType::IndexType idxStart{};
  idxStart[0] = 0;
  idxStart[1] = 0;

  FloatImage2DType::SizeType size2D{};
  size2D[0] = imgDim[idx_hor];
  size2D[1] = imgDim[idxVer];

  FloatImage2DType::SpacingType spacing2D;
  spacing2D[0] = spacing[idx_hor];
  spacing2D[1] = spacing[idxVer];

  FloatImage2DType::PointType origin2D;
  //  origin2D[0] = origin[idxHor];
  //  origin2D[1] = origin[idxVer];
  origin2D[0] = size2D[0] * spacing2D[0] / -2.0;
  origin2D[1] = size2D[1] * spacing2D[1] / -2.0;

  FloatImage2DType::RegionType region;
  region.SetSize(size2D);
  region.SetIndex(idxStart);

  // spTargetImg2D is supposed to be empty.
  if (spTargetImg2D != nullptr) {
    std::cout
        << "something is here in target image. is it gonna be overwritten?"
        << std::endl;
  }

  spTargetImg2D = FloatImage2DType::New();
  spTargetImg2D->SetRegions(region);
  spTargetImg2D->SetSpacing(spacing2D);
  spTargetImg2D->SetOrigin(origin2D);

  spTargetImg2D->Allocate();
  spTargetImg2D->FillBuffer(0);

  // std::cout << "src size = " << spSrcImg3D->GetRequestedRegion().GetSize() <<
  // " " << std::endl;  std::cout << "target image size = " <<
  // spTargetImg2D->GetRequestedRegion().GetSize() << " " << std::endl;

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it_3D(
      spSrcImg3D, spSrcImg3D->GetRequestedRegion());
  // itk::ImageRegionIteratorWithIndex<FloatImageType2D> it_2D (spTargetImg2D,
  // spTargetImg2D->GetRequestedRegion());
  itk::ImageRegionIterator<FloatImage2DType> it_2D(
      spTargetImg2D, spTargetImg2D->GetRequestedRegion());

  it_3D.SetFirstDirection(idx_hor);
  it_3D.SetSecondDirection(idxVer);

  it_3D.GoToBegin();
  it_2D.GoToBegin();

  for (auto i = 0; i < z_size && !it_3D.IsAtEnd(); i++) {
    /*QFileInfo crntFileInfo(arrYKImage[i].m_strFilePath);
    QString crntFileName = crntFileInfo.fileName();
    QString crntPath = strSavingFolder + "/" + crntFileName;*/
    // Search matching slice using slice iterator for m_spProjCTImg
    // std::cout << "Get2DFrom3D: Slide= " << i  << " ";

    if (i == idx) {
      while (!it_3D.IsAtEndOfSlice()) // Error here why?
      {
        while (!it_3D.IsAtEndOfLine()) {
          auto tmpVal = static_cast<float>(
              it_3D.Get()); // in proj image case, this is intensity
          it_2D.Set(tmpVal);
          ++it_2D;
          ++it_3D;
        } // while2
        it_3D.NextLine();
      } // while1
      break;
    } // end if
    it_3D.NextSlice();
  } // end of for

  // std::cout << "cnt = " << cnt << " TotCnt " << cntTot << std::endl;
  /*YK16GrayImage tmpYK;
  tmpYK.UpdateFromItkImageFloat(spTargetImg2D);
  QString str = QString("D:\\testYK\\InsideFunc_%1.raw").arg(idx);
  tmpYK.SaveDataAsRaw(str.toLocal8Bit().constData());*/
}

double GetRawIntensityScaleFactor(QString &strRef_mAs, QString &strCur_mAs) {
  // GetRawIntensity Scale Factor
  auto rawIntensityScaleF = 1.0;

  auto fRef_mAs = 0.0;
  auto fCur_mAs = 0.0;
  // QString strRef_mAs = ui.lineEdit_RefmAs->text();
  auto listmAsRef = strRef_mAs.split(",");
  if (listmAsRef.length() == 2) {
    fRef_mAs = listmAsRef.at(0).toDouble() * listmAsRef.at(1).toDouble();
  }
  // QString strCur_mAs = ui.lineEdit_CurmAs->text();
  auto listmAsCur = strCur_mAs.split(",");
  if (listmAsCur.length() == 2) {
    fCur_mAs = listmAsCur.at(0).toDouble() * listmAsCur.at(1).toDouble();
  }
  if (fRef_mAs * fCur_mAs != 0) {
    rawIntensityScaleF = fRef_mAs / fCur_mAs;
  }

  return rawIntensityScaleF;
  // if 64 40 ref, 40 40 cur --> scaleF = 1.6
  // raw intensity X scaleF ==> raw intensity increased --> this avoids negative
  // scatter map
}

void TransformationRTK2IEC(FloatImageType::Pointer &spSrcTarg) {
  auto sizeOutput = spSrcTarg->GetBufferedRegion().GetSize();
  auto spacingOutput = spSrcTarg->GetSpacing();

  // Transformation is applied
  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;
  // Same image type from original image -3D & float
  FloatImageType::IndexType start_trans{};
  start_trans[0] = 0;
  start_trans[1] = 0;
  start_trans[2] = 0;

  FloatImageType::SizeType size_trans{};
  size_trans[0] = sizeOutput[0]; // X //410
  size_trans[1] = sizeOutput[2]; // Y  // 410
  size_trans[2] = sizeOutput[1]; // Z // 120?

  FloatImageType::SpacingType spacing_trans;
  spacing_trans[0] = spacingOutput[0];
  spacing_trans[1] = spacingOutput[2];
  spacing_trans[2] = spacingOutput[1];

  FloatImageType::PointType Origin_trans;
  Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
  Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
  Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

  FloatImageType::RegionType region_trans;
  region_trans.SetSize(size_trans);
  region_trans.SetIndex(start_trans);

  /* 2) Prepare Target image */
  const auto &targetImg = spSrcTarg;

  /* 3) Configure transform */
  using TransformType = itk::Euler3DTransform<double>;
  auto transform = TransformType::New();

  TransformType::ParametersType param;
  param.SetSize(6);
  // MAXIMUM PARAM NUMBER: 6!!!
  param.put(0, 0.0);                  // rot X // 0.5 = PI/2
  param.put(1, itk::Math::pi / 2.0);  // rot Y
  param.put(2, itk::Math::pi / -2.0); // rot Z
  param.put(3, 0.0);                  // Trans X mm
  param.put(4, 0.0);                  // Trans Y mm
  param.put(5, 0.0);                  // Trans Z mm

  TransformType::ParametersType fixedParam(3); // rotation center
  fixedParam.put(0, 0);
  fixedParam.put(1, 0);
  fixedParam.put(2, 0);

  transform->SetParameters(param);
  transform->SetFixedParameters(fixedParam); // Center of the Transform

  std::cout << "Transform matrix:"
            << "	" << std::endl;
  std::cout << transform->GetMatrix() << std::endl;

  using ResampleFilterType =
      itk::ResampleImageFilter<FloatImageType, FloatImageType>;
  auto resampler = ResampleFilterType::New();
  // FloatImageType::RegionType fixedImg_Region =
  // fixedImg->GetLargestPossibleRegion().GetSize();

  resampler->SetInput(targetImg);
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(targetImg->GetDirection()); // image normal?
  resampler->SetTransform(transform);

  // LR flip

  std::cout << "LR flip filter is being applied" << std::endl;

  using FilterType = itk::FlipImageFilter<FloatImageType>;

  auto flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;

  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;

  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  flipFilter->Update();

  spSrcTarg = flipFilter->GetOutput();
}

QString XML_GetSingleItemString(QXmlStreamReader &xml) {
  QString strResult = "";
  /* We need a start element, like <foo> */
  if (xml.tokenType() != QXmlStreamReader::StartElement) {
    return strResult;
  }

  /* Let's read the name... */
  // auto elementName = xml.name().toString();
  /* ...go to the next. */
  xml.readNext();
  /*
   * This elements needs to contain Characters so we know it's
   * actually data, if it's not we'll leave.
   */
  if (xml.tokenType() != QXmlStreamReader::Characters) {
    return strResult;
  }
  strResult = xml.text().toString();
  return strResult;
}

void AddConstHU(UShortImageType::Pointer &spImg, const int HUval) {

  using iteratorType = itk::ImageRegionIteratorWithIndex<UShortImageType>;
  iteratorType it(spImg, spImg->GetRequestedRegion());

  it.GoToBegin();

  while (!it.IsAtEnd()) {
    const auto crnt_val = static_cast<int>(it.Get());

    auto new_val = HUval + crnt_val;

    if (new_val <= 0) {
      new_val = 0;
    }

    if (new_val >= 4095) {
      new_val = 4095;
    }

    it.Set(static_cast<unsigned short>(new_val));
    ++it;
  }
}

// trans: mm, dicom order
// COuch shift values: directlry come from the INI.XVI file only multiplied
// by 10.0
void ImageTransformUsingCouchCorrection(
    UShortImageType::Pointer &spUshortInput,
    UShortImageType::Pointer &spUshortOutput, const VEC3D couch_trans,
    const VEC3D couch_rot) {
  // couch_trans, couch_rot--> as it is from the text file. only x 10.0 was
  // applied
  if (spUshortInput == nullptr) {
    return;
  }

  using FilterType = itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  auto filter = FilterType::New();

  using TransformType = itk::AffineTransform<double, 3>;
  auto transform = TransformType::New();
  filter->SetTransform(transform);
  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, double>;

  const auto interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);

  filter->SetDefaultPixelValue(0);

  //  const double spacing[3] = { 1.0, 1.0, 1.0 };
  const auto spacing = spUshortInput->GetSpacing();

  filter->SetOutputSpacing(spacing);

  const auto origin = spUshortInput->GetOrigin();

  filter->SetOutputOrigin(origin);

  UShortImageType::DirectionType direction;
  direction.SetIdentity();
  filter->SetOutputDirection(direction);

  const auto size = spUshortInput->GetLargestPossibleRegion().GetSize();
  filter->SetSize(size);
  filter->SetInput(spUshortInput);

  // NOTE: In couch shift reading
  // pTrans->x = couch_Lat_cm*10.0; //sign should be checked
  // pTrans->y = couch_Vert_cm*10.0; //sign should be checked // IEC-->DICOM is
  // already accounted for..but sign!  pTrans->z = couch_Long_cm*10.0; //sign
  // should be checked  pRot->x = couch_Pitch;  pRot->y = couch_Yaw;  pRot->z =
  // couch_Roll;

  TransformType::OutputVectorType translation;
  translation[0] = -couch_trans.x; // X translation in millimeters
  // translation[1] = +couch_trans.y; //so far so good// This is because when
  // IEC->DICOM, sign was not changed during reading the text file
  translation[1] = -couch_trans.y; // Consistent with Tracking software
  translation[2] = -couch_trans.z; // empirically found

  TransformType::OutputVectorType rotation;
  rotation[0] = -couch_rot.x; // X translation in millimeters
  rotation[1] = -couch_rot.y;
  rotation[2] = -couch_rot.z;

  transform->Translate(
      translation); // original position - (couch shift value in DICOM)
  // transform->Rotate3D(rotation);
  filter->Update();

  spUshortOutput = filter->GetOutput();
  // std::cout << "affine transform is successfully done" << std::endl;
}

void RotateImgBeforeFwd(UShortImageType::Pointer &spInputImgUS,
                        UShortImageType::Pointer &spOutputImgUS) {
  if (spInputImgUS == nullptr) {
    std::cout << "ERROR! No 3D image file" << std::endl;
    return;
  }
  // 1) Transform
  auto size_original = spInputImgUS->GetLargestPossibleRegion().GetSize();
  auto spacing_original = spInputImgUS->GetSpacing();

  // Same image type from original image -3D & float
  UShortImageType::IndexType start_trans{};
  start_trans[0] = 0;
  start_trans[1] = 0;
  start_trans[2] = 0;

  UShortImageType::SizeType size_trans{};
  size_trans[0] = size_original[1]; // X //512
  size_trans[1] = size_original[2]; // Y  //512
  size_trans[2] = size_original[0]; // Z //300

  UShortImageType::SpacingType spacing_trans;
  spacing_trans[0] = spacing_original[1];
  spacing_trans[1] = spacing_original[2];
  spacing_trans[2] = spacing_original[0];

  UShortImageType::PointType Origin_trans;
  Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
  Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
  Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

  UShortImageType::RegionType region_trans;
  region_trans.SetSize(size_trans);
  region_trans.SetIndex(start_trans);

  using FilterType = itk::FlipImageFilter<UShortImageType>;
  auto flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;

  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;

  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(spInputImgUS); // plan CT, USHORT image

  using TransformType = itk::Euler3DTransform<double>;
  auto transform = TransformType::New();

  TransformType::ParametersType param;
  param.SetSize(6);
  param.put(0, itk::Math::pi / -2.0); // rot X // 0.5 = PI/2
  param.put(1, 0);                    // rot Y
  param.put(2, itk::Math::pi / 2.0);  // rot Z
  param.put(3, 0.0);                  // Trans X mm
  param.put(4, 0.0);                  // Trans Y mm
  param.put(5, 0.0);                  // Trans Z mm

  TransformType::ParametersType fixedParam(3); // rotation center
  fixedParam.put(0, 0);
  fixedParam.put(1, 0);
  fixedParam.put(2, 0);

  transform->SetParameters(param);
  transform->SetFixedParameters(fixedParam); // Center of the Transform

  /*std::cout << "Transform matrix:" << "	" << std::endl;
  std::cout << transform->GetMatrix() << std::endl;*/

  using ResampleFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  auto resampler = ResampleFilterType::New();

  resampler->SetInput(flipFilter->GetOutput());
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(
      flipFilter->GetOutput()->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  resampler->Update();

  spOutputImgUS = resampler->GetOutput();
}

void ConvertUshort2AttFloat(UShortImageType::Pointer &spImgUshort,
                            FloatImageType::Pointer &spAttImgFloat) {
  using CastFilterType =
      itk::CastImageFilter<UShortImageType,
                           FloatImageType>; // Maybe not inplace filter
  auto castFilter = CastFilterType::New();
  castFilter->SetInput(spImgUshort);

  // Default value
  const auto calibF_A = 1.0;
  const auto calibF_B = 0.0;

  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  auto multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(castFilter->GetOutput());
  multiplyImageFilter->SetConstant(calibF_A / 65535.0);

  using AddImageFilterType =
      itk::AddImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  auto addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput1(multiplyImageFilter->GetOutput());
  const auto addingVal = calibF_B / 65535.0;
  addImageFilter->SetConstant2(addingVal);
  addImageFilter->Update(); // will generate map of real_mu (att.coeff)

  // FloatImageType::Pointer spCTImg_mu;
  spAttImgFloat = multiplyImageFilter->GetOutput();
}
