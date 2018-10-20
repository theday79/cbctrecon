#ifndef CBCTRECON_FDK_HXX
#define CBCTRECON_FDK_HXX

// Qt
#include <QDir>
#include <QFileInfo>
#include <QString> // for QString

// ITK
#include "itkAbsImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkMedianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkStreamingImageFilter.h"
#include "itkThresholdImageFilter.h"

// RTK
#include "rtkBackProjectionImageFilter.h" // for BackProje...
#include "rtkConstantImageSource.h"
#include "rtkDisplacedDetectorImageFilter.h"
#include "rtkFDKConeBeamReconstructionFilter.h" // for FDKConeBeamReconstr...
#include "rtkFFTProjectionsConvolutionImageFilter.h" // for FFTProjec...
#include "rtkFFTRampImageFilter.h"                   // for FFTRampImageFilter
#include "rtkFieldOfViewImageFilter.h"
#include "rtkMacro.h" // for TRY_AND_EXIT_ON_ITK_EXCEPTION
#include "rtkParkerShortScanImageFilter.h"
#include "rtkThreeDCircularProjectionGeometry.h" // for ThreeDCircularProje...

#if USE_OPENCL_RTK
#include "rtkOpenCLFDKConeBeamReconstructionFilter.h"
#endif

#if USE_CUDA
#include <rtkCudaDisplacedDetectorImageFilter.h>
#include <rtkCudaFDKConeBeamReconstructionFilter.h>
#include <rtkCudaParkerShortScanImageFilter.h>
#endif // USE_CUDA

#include "YK16GrayImage.h"
#include "cbctrecon_compute.h"


template <typename ImageType>
typename ImageType::Pointer RTKOpenCLFDK(
    const typename ImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    typename ImageType::SpacingType spacing,
    typename ImageType::SizeType sizeOutput, FDK_options fdk_options) {

  try {
    spCurImg->Update();
  } catch (const std::exception &err) {
    std::cerr << "Couldn't update spCurImg: " << err.what() << std::endl;
  }
  using CastFilterType = itk::CastImageFilter<ImageType, FloatImageType>;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(spCurImg);
  castFilter->Update();

  // Generate image sources for cone beam CT reconstruction
  using ConstantImageSourceType = rtk::ConstantImageSource<ImageType>;
  typename ConstantImageSourceType::PointType origin;

  origin[0] = -0.5 * sizeOutput[0] * spacing[0]; // Y in DCM?
  origin[1] = -0.5 * sizeOutput[1] * spacing[1]; // Z in DCM?
  origin[2] = -0.5 * sizeOutput[2] * spacing[2]; // X in DCM?

  typename ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();
  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);
  constantImageSource->SetSize(sizeOutput);
  constantImageSource->SetConstant(0.0); // initial value
  // constantImageSource->Update();

  // FDK reconstruction filtering
  using FDKOPENCLType = rtk::OpenCLFDKConeBeamReconstructionFilter;
  FDKOPENCLType::Pointer feldkampOCL = FDKOPENCLType::New();

  feldkampOCL->SetInput(0, constantImageSource->GetOutput());
  feldkampOCL->SetInput(1, castFilter->GetOutput());
  feldkampOCL->SetGeometry(m_spCustomGeometry);
  feldkampOCL->GetRampFilter()->SetTruncationCorrection(
      fdk_options.TruncCorFactor);
  feldkampOCL->GetRampFilter()->SetHannCutFrequency(fdk_options.HannCutX);
  feldkampOCL->GetRampFilter()->SetHannCutFrequencyY(fdk_options.HannCutY);
  feldkampOCL->GetRampFilter()->SetCosineCutFrequency(fdk_options.CosCut);
  feldkampOCL->GetRampFilter()->SetHammingFrequency(fdk_options.HammCut);

  // feldkampOCL->Update();
  // feldkampOCL->PrintTiming(std::cout); Deprecated in rtk 1.4

  using CastFilterType2 = itk::CastImageFilter<FloatImageType, ImageType>;
  typename CastFilterType2::Pointer castFilter2 = CastFilterType2::New();
  castFilter2->SetInput(feldkampOCL->GetOutput());
  castFilter2->Update();
  return castFilter2->GetOutput();
}

template <enDeviceType Tdev>
void CbctRecon::DoReconstructionFDK(enREGI_IMAGES target,
                                    FDK_options fdk_options) {
  if (Tdev == CUDA_DEVT) {
#if USE_CUDA
    using CUDAFloatImageType = itk::CudaImage<float, 3U>;
    using DDFType = rtk::CudaDisplacedDetectorImageFilter;
    using PSSFType = rtk::CudaParkerShortScanImageFilter;
    using FDKType = rtk::CudaFDKConeBeamReconstructionFilter;
    DoReconstructionFDK<Tdev, CUDAFloatImageType, DDFType, PSSFType, FDKType>(
        target, fdk_options);
#endif
  } else {
    using DDFType = rtk::DisplacedDetectorImageFilter<FloatImageType>;
    using PSSFType = rtk::ParkerShortScanImageFilter<FloatImageType>;
    using FDKType = rtk::FDKConeBeamReconstructionFilter<FloatImageType>;
    DoReconstructionFDK<Tdev, FloatImageType, DDFType, PSSFType, FDKType>(
        target, fdk_options);
  }
}

template <enDeviceType Tdev, typename ImageType, typename DDFType,
          typename PSSFType, typename FDKType>
void CbctRecon::DoReconstructionFDK(enREGI_IMAGES target,
                                    FDK_options fdk_options) {

  using castToImageType = itk::CastImageFilter<FloatImageType, ImageType>;
  typename castToImageType::Pointer castfilter = castToImageType::New();
  castfilter->SetInput(m_spProjImg3DFloat);
  castfilter->Update();
  auto cuda_spProjImg3DFloat = castfilter->GetOutput();

  if (cuda_spProjImg3DFloat == nullptr) {
    std::cout << "processed Projection image is not ready yet" << std::endl;
    return;
  }

  std::cout << "CUDA method will be used..." << std::endl;

  using DuplicatorType = itk::ImageDuplicator<ImageType>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(cuda_spProjImg3DFloat);
  duplicator->Update();

  typename ImageType::Pointer spCurImg =
      duplicator->GetOutput(); // already down sampled

  // Displaced detector weighting // set pipeline //inplace filter
  typename DDFType::Pointer ddf = DDFType::New();
  if (fdk_options.displacedDetectorFilter) {
    ddf->SetInput(spCurImg);
    ddf->SetGeometry(m_spCustomGeometry);
    std::cout << "DDF was set in pipeline" << std::endl;

    if (fdk_options.updateAfterDDF) {
      ddf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = ddf->GetOutput();
  }

  typename PSSFType::Pointer pssf = PSSFType::New();
  if (fdk_options.ParkerShortScan) {
    // Short scan image filter
    pssf->SetInput(spCurImg);
    // pssf->SetGeometry( geometryReader->GetOutputObject() );
    pssf->SetGeometry(m_spCustomGeometry);
    // pssf->InPlaceOff(); //YKComments: Do not overwrite input image buffer for
    // output

    pssf->GlobalWarningDisplayOff(); // I don't care about any of the potential
                                     // warnings

    if (fdk_options.updateAfterDDF) {
      pssf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = pssf->GetOutput();
  }

  // Just Before going to FDK recon,
  // Update Projection data and delete old data.

  // Let's duplicate this
  // Original m_spProjImg3D will be deleted after update

  if (fdk_options.updateAfterDDF) {
    using DuplicatorType = itk::ImageDuplicator<ImageType>;
    typename DuplicatorType::Pointer ImDuplicator = DuplicatorType::New();
    ImDuplicator->SetInputImage(spCurImg);
    ImDuplicator->Update();

    using CastFilterType = itk::CastImageFilter<ImageType, FloatImageType>;
    typename CastFilterType::Pointer CastFilter = CastFilterType::New();
    CastFilter->SetInput(ImDuplicator->GetOutput());
    CastFilter->Update();

    m_spProjImg3DFloat = CastFilter->GetOutput();

    SetMaxAndMinValueOfProjectionImage(); // scan m_spProjImg3D and update
                                          // m_fProjImgValueMin, max
  }

  // Generate image sources for cone beam CT reconstruction
  using ConstantImageSourceType = rtk::ConstantImageSource<ImageType>;

  typename ConstantImageSourceType::SizeType sizeOutput{};
  sizeOutput[0] = fdk_options.ct_size[0]; // pixel
  sizeOutput[1] =
      fdk_options
          .ct_size[1]; // Caution!: direction is different in NKI SCAN FIle
  sizeOutput[2] = fdk_options.ct_size[2];

  typename ConstantImageSourceType::SpacingType spacing;
  spacing[0] = fdk_options.ct_spacing[0];
  spacing[1] = fdk_options.ct_spacing[1];
  spacing[2] = fdk_options.ct_spacing[2];

  if (GetOutputResolutionFromFOV<ConstantImageSourceType, FloatImageType>(
          sizeOutput, spacing, m_spCustomGeometry, m_spProjImg3DFloat,
          fdk_options.outputFilePath)) {
    std::cout << "Reconstruction resolution and image size were set "
                 "automatically, as no outputfilepath was given."
              << std::endl;
  }

  if (fdk_options.TruncCorFactor > 0.0 && target == REGISTER_COR_CBCT) {
    std::cout << "Warning! Truncation factor is " << fdk_options.TruncCorFactor
              << ". Regardless of previous setting, this factor should not be "
                 "0 for scatter corrected CBCT. Now zero value is applied."
              << std::endl;
    fdk_options.TruncCorFactor = 0.0;
  }

  typename ImageType::Pointer targetImg;

  if (Tdev == OPENCL_DEVT) {

    std::cout << "Starting RTK fdk" << std::endl;
    targetImg = RTKOpenCLFDK<ImageType>(spCurImg, m_spCustomGeometry, spacing,
                                        sizeOutput, fdk_options);
  } else {

    typename ConstantImageSourceType::PointType origin;
    origin[0] = -0.5 * sizeOutput[0] * spacing[0]; // Y in DCM?
    origin[1] = -0.5 * sizeOutput[1] * spacing[1]; // Z in DCM?
    origin[2] = -0.5 * sizeOutput[2] * spacing[2]; // X in DCM?

    typename ConstantImageSourceType::Pointer constantImageSource =
        ConstantImageSourceType::New();
    constantImageSource->SetOrigin(origin);
    constantImageSource->SetSpacing(spacing);
    constantImageSource->SetSize(sizeOutput);
    constantImageSource->SetConstant(0.0); // initial value

    using StreamerType = itk::StreamingImageFilter<ImageType, ImageType>;
    // FDK reconstruction filtering

    typename FDKType::Pointer feldkamp = FDKType::New();

    std::cout << "CUDA will be used for FDK reconstruction" << std::endl;
    feldkamp->SetInput(0, constantImageSource->GetOutput());
    feldkamp->SetInput(1, spCurImg);
    feldkamp->SetGeometry(m_spCustomGeometry);
    feldkamp->GetRampFilter()->SetTruncationCorrection(
        fdk_options.TruncCorFactor);
    feldkamp->GetRampFilter()->SetHannCutFrequency(fdk_options.HannCutX);
    feldkamp->GetRampFilter()->SetCosineCutFrequency(fdk_options.CosCut);
    feldkamp->GetRampFilter()->SetHammingFrequency(fdk_options.HammCut);
    feldkamp->GetRampFilter()->SetHannCutFrequencyY(fdk_options.HannCutY);

    if (fdk_options.updateAfterDDF) {
      feldkamp->Update();
    } else {
      feldkamp->UpdateOutputInformation();
    }

    // Motion compensated CBCT settings
    // if(args_info.signal_given && args_info.dvf_given)
    //{
    // dvfReader->SetFileName(args_info.dvf_arg);
    // def->SetSignalFilename(args_info.signal_arg);
    // dynamic_cast<FDKCPUType*>(feldkamp.GetPointer())->SetBackProjectionFilter(
    // bp.GetPointer() );
    //}

    std::cout << "Cone beam reconstruction pipeline is ready" << std::endl;

    // Streaming depending on streaming capability of writer --> not affect the
    // calc. speed

    using StreamerType = itk::StreamingImageFilter<ImageType, ImageType>;
    typename StreamerType::Pointer streamerBP = StreamerType::New();

    streamerBP->SetInput(feldkamp->GetOutput());
    streamerBP->SetNumberOfStreamDivisions(4); // YK: 1 in example code from
                                               // "rtkfdk" //AG: stated in test:
                                               // 4 for ITK MAJOR >= 4
    targetImg = streamerBP->GetOutput();
  }
  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;

  using FOVfilterType = rtk::FieldOfViewImageFilter<ImageType, ImageType>;
  typename FOVfilterType::Pointer FOVfilter = FOVfilterType::New();
  FOVfilter->SetGeometry(m_spCustomGeometry);
  FOVfilter->SetInput(0, targetImg);
  FOVfilter->SetProjectionsStack(spCurImg.GetPointer());

  /* RTK-produced 3D Volume should be changed in coordination of itk */
  /* Coordination transformation using Euler 3D transformation */

  // 1) Prepare Canvas parameter
  // FloatImageType::Pointer fixedImg = FloatImageType::New();
  // start index: What is the index of Left Top Inferior corner in DICOM
  // coordinate?

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
  spacing_trans[0] = spacing[0];
  spacing_trans[1] = spacing[2];
  spacing_trans[2] = spacing[1];

  FloatImageType::PointType Origin_trans;
  Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
  Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
  Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

  FloatImageType::RegionType region_trans;
  region_trans.SetSize(size_trans);
  region_trans.SetIndex(start_trans);

  /* 2) Prepare Target image */
  // ImageType::Pointer
  targetImg = FOVfilter->GetOutput(); // streamerBP->GetOutput();

  /* 3) Configure transform */
  //                            | 0  1  0 |
  // Rz(-pi/2)*Rx(0)*Ry(pi/2) = | 0  0 -1 |
  //                            |-1  0  0 |
  itk::Matrix<double, 3, 3> CoordChangeMatrix;
  // 1st row
  CoordChangeMatrix[0][0] = 0.0;
  CoordChangeMatrix[0][1] = 1.0;
  CoordChangeMatrix[0][2] = 0.0;
  // 2nd row
  CoordChangeMatrix[1][0] = 0.0;
  CoordChangeMatrix[1][1] = 0.0;
  CoordChangeMatrix[1][2] = -1.0;
  // 3rd row
  CoordChangeMatrix[2][0] = -1.0;
  CoordChangeMatrix[2][1] = 0.0;
  CoordChangeMatrix[2][2] = 0.0;

  const itk::Vector<double, 3U> offset(0.0);

  using TransformType = itk::MatrixOffsetTransformBase<double, 3U, 3U>;
  TransformType::Pointer transform = TransformType::New();
  transform->SetMatrix(CoordChangeMatrix);
  transform->SetOffset(offset);

  using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetInput(targetImg);
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(targetImg->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  // resampler->Update();//yktemp Error 2

  // LR flip

  std::cout << "Flip filter is being applied" << std::endl;
  using FlipFilterType = itk::FlipImageFilter<ImageType>;
  typename FlipFilterType::Pointer flipFilter = FlipFilterType::New();
  using FlipAxesArrayType = typename FlipFilterType::FlipAxesArrayType;
  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;
  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  using AbsImageFilterType = itk::AbsImageFilter<ImageType, ImageType>;
  typename AbsImageFilterType::Pointer absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(
      flipFilter
          ->GetOutput()); // 20140206 modified it was a buug
                          // absImgFilter->SetInput(resampler->GetOutput());

  using MultiplyImageFilterType = itk::MultiplyImageFilter<ImageType>;
  typename MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  using CastFilterType = itk::CastImageFilter<ImageType, UShortImageType>;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());
  try {
    castFilter->Update(); // YK20150109
  } catch (const std::exception &e) {
    std::cerr << "std::exception thrown: " << e.what() << std::endl;
    return;
  }

  UShortImageType::Pointer tmpReconImg;
  // if all 0 0 0 don't do the median filtering

  std::cout << "Reconstructing the image.. please wait..." << std::endl;
  UShortImageType::SizeType indexRadius{};
  indexRadius[0] = fdk_options.medianRadius[0]; // radius along x
  indexRadius[1] = fdk_options.medianRadius[1]; // radius along y
  indexRadius[2] = fdk_options.medianRadius[2]; // radius along z
  if (fdk_options.medianFilter &&
      (indexRadius[0] != 0 || indexRadius[1] != 0 || indexRadius[2] != 0)) {
    using MedianFilterType =
        itk::MedianImageFilter<UShortImageType,
                               UShortImageType>; // TODO(AGA): CUDA THIS!!

    MedianFilterType::Pointer medFilter = MedianFilterType::New();
    // this is radius. 1 --> median window 3
    std::cout << "Post median(3D) filtering is in the pipeline..Size(radius X "
                 "Y Z) is = "
              << indexRadius << std::endl;

    medFilter->SetRadius(indexRadius);
    medFilter->SetInput(castFilter->GetOutput());
    medFilter->Update(); // Error here!g

    tmpReconImg = medFilter->GetOutput();
    std::cout << "median filtering has been done" << std::endl;
  } else {
    std::cout << "No post median filtering is used" << std::endl;
    castFilter->Update();
    tmpReconImg = castFilter->GetOutput();
  }

  using ThresholdImageFilterType = itk::ThresholdImageFilter<UShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilterAbove =
      ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(tmpReconImg);
  thresholdFilterAbove->ThresholdAbove(4095);
  thresholdFilterAbove->SetOutsideValue(4095);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
      ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(0);
  thresholdFilterBelow->SetOutsideValue(0);
  thresholdFilterBelow->Update();
  tmpReconImg = thresholdFilterBelow->GetOutput();

  std::cout << "After Filtering" << std::endl;

  // By default CanRunInPlace checks whether the input and output image type
  // match.
  switch (target) {
  case REGISTER_RAW_CBCT:
    m_spRawReconImg = tmpReconImg; // Checked.. successfully alive.
    m_spCrntReconImg = m_spRawReconImg;
    break;
  case REGISTER_COR_CBCT:
    m_spScatCorrReconImg = tmpReconImg; // Checked.. successfully alive.
    m_spCrntReconImg = m_spScatCorrReconImg;
    break;
  default:
    std::cerr << "You are using a non-valid target!" << std::endl;
    return;
  }

  QFileInfo outFileInfo(fdk_options.outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  if (fdk_options.outputFilePath.length() < 2 || !outFileDir.exists()) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(fdk_options.outputFilePath.toLocal8Bit().constData());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(m_spCrntReconImg);

    std::cout << "Writing the image to: "
              << fdk_options.outputFilePath.toLocal8Bit().constData()
              << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());
    std::cout << std::endl;
    std::cout << "Output generation was succeeded" << std::endl;
  }

  m_dspYKReconImage->CreateImage(size_trans[0], size_trans[1], 0);

  using CastBackType = itk::CastImageFilter<ImageType, FloatImageType>;
  typename CastBackType::Pointer castBack = CastBackType::New();
  castBack->SetInput(cuda_spProjImg3DFloat);
  castBack->Update();
  m_spProjImg3DFloat = castBack->GetOutput();

  std::cout << "FINISHED!: FDK CBCT reconstruction" << std::endl;
}

#endif
