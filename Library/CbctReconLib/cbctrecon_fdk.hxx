#ifndef CBCTRECON_FDK_HXX
#define CBCTRECON_FDK_HXX

// Qt
#include <QString> // for QString

// ITK
#include "itkAbsImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkMedianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkStreamingImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkTimeProbe.h"

// RTK
#include "rtkConstantImageSource.h"
#include "rtkDisplacedDetectorImageFilter.h"
#include "rtkFDKConeBeamReconstructionFilter.h" // for FDKConeBeamReconstr...
#include "rtkFFTProjectionsConvolutionImageFilter.h" // for FFTProjec...
#include "rtkFFTRampImageFilter.h"                   // for FFTRampImageFilter
#include "rtkFieldOfViewImageFilter.h"
#include "rtkMacro.h" // for TRY_AND_EXIT_ON_ITK_EXCEPTION
#include "rtkParkerShortScanImageFilter.h"
#include "rtkThreeDCircularProjectionGeometry.h" // for ThreeDCircularProje...

#if RTK_USE_OPENCL
#include "rtkOpenCLFDKConeBeamReconstructionFilter.h"
#include "rtkOpenCLForwardProjectionImageFilter.h"
#else
#include "rtkJosephForwardProjectionImageFilter.h"
#endif

#if USE_CUDA
#include "rtkCudaForwardProjectionImageFilter.h"
#include <rtkCudaDisplacedDetectorImageFilter.h>
#include <rtkCudaFDKConeBeamReconstructionFilter.h>
#include <rtkCudaParkerShortScanImageFilter.h>
#endif // USE_CUDA

#include "YK16GrayImage.h"
#include "cbctrecon_compute.h"
#include "cbctrecon_io.h"

using namespace std::literals;

#ifdef RTK_USE_OPENCL
template <typename ImageType>
typename ImageType::Pointer RTKOpenCLFDK(
    const typename ImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    typename ImageType::SpacingType spacing,
    typename ImageType::SizeType sizeOutput, FDK_options const &fdk_options) {

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
  auto feldkampOCL = FDKOPENCLType::New();

  feldkampOCL->SetInput(0, constantImageSource->GetOutput());
  feldkampOCL->SetInput(1, castFilter->GetOutput());
  feldkampOCL->SetGeometry(m_spCustomGeometry);
  auto p_rampfilter = feldkampOCL->GetRampFilter().GetPointer();
  p_rampfilter->SetTruncationCorrection(fdk_options.TruncCorFactor);
  p_rampfilter->SetHannCutFrequency(fdk_options.HannCutX);
  p_rampfilter->SetHannCutFrequencyY(fdk_options.HannCutY);
  p_rampfilter->SetCosineCutFrequency(fdk_options.CosCut);
  p_rampfilter->SetHammingFrequency(fdk_options.HammCut);

  // feldkampOCL->Update();
  // feldkampOCL->PrintTiming(std::cout); Deprecated in rtk 1.4

  using CastFilterType2 = itk::CastImageFilter<FloatImageType, ImageType>;
  typename CastFilterType2::Pointer castFilter2 = CastFilterType2::New();
  castFilter2->SetInput(feldkampOCL->GetOutput());
  castFilter2->Update();
  return castFilter2->GetOutput();
}
#endif

template <enDeviceType Tdev>
void CbctRecon::DoReconstructionFDK(const enREGI_IMAGES target,
                                    const FDK_options &fdk_options) {
  if (Tdev == enDeviceType::CUDA_DEVT) {
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
void CbctRecon::DoReconstructionFDK(const enREGI_IMAGES target,
                                    const FDK_options &fdk_options) {

  using castToImageType = itk::CastImageFilter<FloatImageType, ImageType>;
  typename castToImageType::Pointer castfilter = castToImageType::New();
  castfilter->SetInput(m_spProjImg3DFloat);
  castfilter->Update();
  auto cuda_spProjImg3DFloat = castfilter->GetOutput();

  if (cuda_spProjImg3DFloat == nullptr) {
    std::cout << "processed Projection image is not ready yet" << std::endl;
    return;
  }

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

  if (crl::GetOutputResolutionFromFOV<ConstantImageSourceType, FloatImageType>(
          sizeOutput, spacing, m_spCustomGeometry, m_spProjImg3DFloat,
          fdk_options.outputFilePath)) {
    std::cout << "Reconstruction resolution and image size were set "
                 "automatically, as no outputfilepath was given."
              << std::endl;
  }

  auto trunc_factor = fdk_options.TruncCorFactor;
  if (fdk_options.TruncCorFactor > 0.0 &&
      target == enREGI_IMAGES::REGISTER_COR_CBCT) {
    std::cout << "Warning! Truncation factor is " << fdk_options.TruncCorFactor
              << ". Regardless of previous setting, this factor should not be "
                 "0 for scatter corrected CBCT. Now zero value is applied."
              << std::endl;
    trunc_factor = 0.0;
  }

  typename ImageType::Pointer targetImg;

  if (Tdev == enDeviceType::OPENCL_DEVT) {
#ifdef RTK_USE_OPENCL
    std::cout << "Starting RTK fdk" << std::endl;
    targetImg = RTKOpenCLFDK<ImageType>(spCurImg, m_spCustomGeometry, spacing,
                                        sizeOutput, fdk_options);
#else
    std::cerr << "You did not compile with RTK_USE_OPENCL=ON!\n";
    return;
#endif
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
    auto p_rampfilter = feldkamp->GetRampFilter().GetPointer();
    p_rampfilter->SetTruncationCorrection(trunc_factor);
    p_rampfilter->SetHannCutFrequency(fdk_options.HannCutX);
    p_rampfilter->SetHannCutFrequencyY(fdk_options.HannCutY);
    p_rampfilter->SetCosineCutFrequency(fdk_options.CosCut);
    p_rampfilter->SetHammingFrequency(fdk_options.HammCut);

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
    streamerBP->Update();
    targetImg = streamerBP->GetOutput();
  }
  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;

  using FOVfilterType = rtk::FieldOfViewImageFilter<ImageType, FloatImageType>;
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
  // FloatImageType::Pointer
  // targetImg = FOVfilter->GetOutput(); // streamerBP->GetOutput();

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
  auto transform = TransformType::New();
  transform->SetMatrix(CoordChangeMatrix);
  transform->SetOffset(offset);

  using ResampleFilterType =
      itk::ResampleImageFilter<FloatImageType, FloatImageType>;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetInput(FOVfilter->GetOutput());
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(
      FOVfilter->GetOutput()->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  // resampler->Update();//yktemp Error 2

  // LR flip

  std::cout << "Flip filter is being applied" << std::endl;
  using FlipFilterType = itk::FlipImageFilter<FloatImageType>;
  typename FlipFilterType::Pointer flipFilter = FlipFilterType::New();
  using FlipAxesArrayType = typename FlipFilterType::FlipAxesArrayType;
  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;
  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  using AbsImageFilterType =
      itk::AbsImageFilter<FloatImageType, FloatImageType>;
  typename AbsImageFilterType::Pointer absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(
      flipFilter
          ->GetOutput()); // 20140206 modified it was a buug
                          // absImgFilter->SetInput(resampler->GetOutput());

  using MultiplyImageFilterType = itk::MultiplyImageFilter<FloatImageType>;
  typename MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
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
                               UShortImageType>; // TODO(AGA): OpenCL THIS!!

    auto medFilter = MedianFilterType::New();
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
  auto thresholdFilterAbove = ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(tmpReconImg);
  thresholdFilterAbove->ThresholdAbove(4095);
  thresholdFilterAbove->SetOutsideValue(4095);

  auto thresholdFilterBelow = ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(0);
  thresholdFilterBelow->SetOutsideValue(0);
  thresholdFilterBelow->Update();
  tmpReconImg = thresholdFilterBelow->GetOutput();

  std::cout << "After Filtering" << std::endl;

  // By default CanRunInPlace checks whether the input and output image type
  // match.
  switch (target) {
  case enREGI_IMAGES::REGISTER_RAW_CBCT:
    m_spRawReconImg = std::move(tmpReconImg); // Checked.. successfully alive.
    m_spCrntReconImg = m_spRawReconImg;
    break;
  case enREGI_IMAGES::REGISTER_COR_CBCT:
    m_spScatCorrReconImg =
        std::move(tmpReconImg); // Checked.. successfully alive.
    m_spCrntReconImg = m_spScatCorrReconImg;
    break;
  default:
    std::cerr << "You are using a non-valid target!" << std::endl;
    return;
  }

  fs::path outFileInfo = fdk_options.outputFilePath;
  auto outFileDir = fs::absolute(outFileInfo);

  if (fdk_options.outputFilePath.empty() || !fs::exists(outFileDir)) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    auto writer = WriterType::New();
    writer->SetFileName(fdk_options.outputFilePath.string());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(m_spCrntReconImg);

    std::cout << "Writing the image to: " << fdk_options.outputFilePath
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

// output spProjCT3D => intensity value, not line integral
template <typename CTImageType>
FloatImageType::Pointer
CbctRecon::ForwardProjection_master(typename CTImageType::Pointer &spVolImg3D,
                                    GeometryType::Pointer &spGeometry,
                                    const bool bSave, const bool use_cuda) {
  if (spVolImg3D == nullptr) {
    std::cout << "ERROR! No 3D-CT file. Load 3D CT file first" << std::endl;
    return nullptr;
  }

  if (this->m_iCntSelectedProj < 1 && bSave) {
    std::cout << "Error! No projection image is loaded" << std::endl;
    return nullptr;
  }

  if (spGeometry->GetGantryAngles().empty()) {
    std::cout << "No geometry!" << std::endl;
    return nullptr;
  }

  FloatImageType::Pointer spProj3D;
#if USE_CUDA
  if (use_cuda) {
    spProj3D = this->ForwardProjection<CUDAFloatImageType>(
        spVolImg3D, spGeometry); // final moving image
  } else
#else
  if (use_cuda) {
    std::cerr << "USE_CUDA not defined at compiletime, using CPU or OpenCL "
                 "implementation!\n";
  }

#endif
  {
    spProj3D = this->ForwardProjection<FloatImageType>(
        spVolImg3D, spGeometry); // final moving image
  }
  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << " Patient DIR Path: " << this->m_strPathPatientDir
              << std::endl;

    auto manuallySelectedDir = true; // <- just to make sure I don't break
                                     // usecases of the older version.

    if (this->m_strPathPatientDir.empty()) {
      std::cerr << "No patient DIR name, using plm_tmp\n";
      manuallySelectedDir = false;
    }

    if (!manuallySelectedDir) {
      this->m_strPathPatientDir = "plm_tmp";
    }

    // Get current folder
    const auto subdir_images = "IMAGES"s;
    const auto strCrntDir =
        this->m_strPathPatientDir / subdir_images; // current Proj folder

    // Make a sub directory
    auto crntDir = strCrntDir;

    if (!fs::exists(crntDir) && m_projFormat == enProjFormat::HIS_FORMAT) {
      if (manuallySelectedDir) {
        auto current_dir = this->m_strPathPatientDir;
        const auto success = fs::create_directory(current_dir / subdir_images);
        if (!success) {
          std::cerr << "Could not create subfolder IMAGES in given directory"
                    << std::endl;
          return spProj3D;
        }
      } else if (m_projFormat == enProjFormat::HIS_FORMAT) {
        std::cout << "File save error: The specified folder does not exist."
                  << std::endl;
        return spProj3D;
      }
    } else {
      // Odds are that just the IMAGES subdirectory didn't exists
      crntDir = this->m_strPathPatientDir;
      if (!fs::exists(crntDir)) {
        std::cerr << "Non-existent path provided, not saving fwd projs\n";
        return spProj3D;
      }
    }

    const auto fwdDirName = "fwd_" + this->m_strDCMUID;

    const auto tmpResult = fs::create_directory(
        crntDir / fwdDirName); // what if the directory exists?

    if (!tmpResult) {
      std::cout << "FwdProj directory seems to exist already. Files will be "
                   "overwritten."
                << std::endl;
    }

    auto strSavingFolder = fs::absolute(crntDir) / fwdDirName;
    if (m_projFormat == enProjFormat::HIS_FORMAT) {
      this->SaveProjImageAsHIS(spProj3D, this->m_arrYKBufProj, strSavingFolder,
                               this->m_fResampleF);
    } else {
      auto fn_fwd_prj = "fwd_proj.mha";
      if (m_spRawReconImg->GetBufferPointer() ==
          spVolImg3D->GetBufferPointer()) {
        fn_fwd_prj = "fwd_proj_rawrec.mha";
      }
      crl::saveImageAsMHA<FloatImageType>(
          spProj3D, fs::absolute(strSavingFolder / fn_fwd_prj).string());
    }
  }
  return spProj3D;
}

template <typename ImageType> struct forward_projector {
#ifdef RTK_USE_OPENCL
  using type = rtk::OpenCLForwardProjectionImageFilter<ImageType, ImageType>;
#else
  using type = rtk::JosephForwardProjectionImageFilter<ImageType, ImageType>;
#endif
  // forwardProjection =
  // rtk::RayCastInterpolatorForwardProjectionImageFilter<FloatImageType,
  // FloatImageType>::New();
};
#if USE_CUDA
template <> struct forward_projector<CUDAFloatImageType> {
  using type = rtk::CudaForwardProjectionImageFilter<CUDAFloatImageType,
                                                     CUDAFloatImageType>;
};
#endif

template <typename DevFloatImageType>
FloatImageType::Pointer
CbctRecon::ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                             GeometryType::Pointer &spGeometry) const {

  // m_spProjCTImg --> spProjCT3D

  FloatImageType::Pointer spResultProjImageFloat;
  // Euler Transformation for RTK's weird orientation

  // int iNumOfProjections = 0;

  {
    // 0) CT image Transformation
    auto size_original = spVolImg3D->GetLargestPossibleRegion().GetSize();
    auto spacing_original = spVolImg3D->GetSpacing();

    // Same image type from original image -3D & float
    UShortImageType::IndexType start_trans{};
    start_trans[0] = 0;
    start_trans[1] = 0;
    start_trans[2] = 0;

    UShortImageType::SizeType size_trans{};
    size_trans[0] = size_original[1]; // X //512
    size_trans[1] = size_original[2]; // Y  //512
    size_trans[2] = size_original[0]; // Z //300

    std::cout << " size_trans" << size_trans << std::endl;

    UShortImageType::SpacingType spacing_trans;
    spacing_trans[0] = spacing_original[1];
    spacing_trans[1] = spacing_original[2];
    spacing_trans[2] = spacing_original[0];

    std::cout << " spacing_trans" << spacing_trans << std::endl;

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
    flipFilter->SetInput(spVolImg3D); // plan CT, USHORT image

    //                            | 0  0 -1 |
    // Rz(pi/2)*Rx(-pi/2)*Ry(0) = | 1  0  0 |
    //                            | 0 -1  0 |
    itk::Matrix<double, 3, 3> CoordChangeMatrix;
    // 1st row
    CoordChangeMatrix[0][0] = 0.0;
    CoordChangeMatrix[0][1] = 0.0;
    CoordChangeMatrix[0][2] = -1.0;
    // 2nd row
    CoordChangeMatrix[1][0] = 1.0;
    CoordChangeMatrix[1][1] = 0.0;
    CoordChangeMatrix[1][2] = 0.0;
    // 3rd row
    CoordChangeMatrix[2][0] = 0.0;
    CoordChangeMatrix[2][1] = -1.0;
    CoordChangeMatrix[2][2] = 0.0;

    const itk::Vector<double, 3U> offset(0.0);

    using TransformType = itk::MatrixOffsetTransformBase<double, 3U, 3U>;
    auto transform = TransformType::New();
    transform->SetMatrix(CoordChangeMatrix);
    transform->SetOffset(offset);

    std::cout << "Transform matrix:"
              << "\n"
              << transform->GetMatrix() << std::endl;

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

    using CastFilterType =
        itk::CastImageFilter<UShortImageType,
                             FloatImageType>; // Maybe not inplace filter
    auto castFilter = CastFilterType::New();
    castFilter->SetInput(resampler->GetOutput());

    // Default value
    const auto calibF_A = 1.0 / std::numeric_limits<unsigned short>::max();

    using MultiplyImageFilterType =
        itk::MultiplyImageFilter<FloatImageType, FloatImageType,
                                 DevFloatImageType>;
    auto multiplyImageFilter = MultiplyImageFilterType::New();
    multiplyImageFilter->SetInput(castFilter->GetOutput());
    multiplyImageFilter->SetConstant(calibF_A);
    multiplyImageFilter->Update(); // will generate map of real_mu (att.coeff)

    const auto spCTImg_mu = multiplyImageFilter->GetOutput();

    // 2) Prepare empty projection images //Should be corresonponding to raw
    // projection images

    // Create a stack of empty projection images
    using ConstantImageSourceType =
        rtk::ConstantImageSource<DevFloatImageType>; // Output: FLoat image =
                                                     // may be mu_t =
                                                     // log(I_0/I)
    auto constantImageSource = ConstantImageSourceType::New();

    typename ConstantImageSourceType::SizeType size{};
    typename ConstantImageSourceType::SpacingType spacing;
    typename ConstantImageSourceType::PointType origin;

    std::cout << "Setting-up vacant projection image data" << std::endl;

    // a) size
    // std::cout << "chk1" << std::endl;
    size[0] = m_spProjImg3DFloat->GetBufferedRegion()
                  .GetSize()[0]; // crl::ce_round((double)DEFAULT_W*m_fResampleF);
    size[1] = m_spProjImg3DFloat->GetBufferedRegion()
                  .GetSize()[1]; // crl::ce_round((double)DEFAULT_H*m_fResampleF);
    size[2] = spGeometry->GetGantryAngles().size();
    // iNumOfProjections = size[2];

    // b) spacing
    spacing[0] =
        m_spProjImg3DFloat->GetSpacing()[0]; // m_fProjSpacingX / m_fResampleF;
                                             // // typical HIS file
    spacing[1] =
        m_spProjImg3DFloat->GetSpacing()[1]; // m_fProjSpacingY / m_fResampleF;
    spacing[2] = 1.0;

    // c) Origin: can center be the image center? or should be related to the CT
    // image???
    origin[0] = spacing[0] * (size[0] - 1) * -0.5;
    origin[1] = spacing[1] * (size[1] - 1) * -0.5;
    origin[2] = 0.0;

    constantImageSource->SetOrigin(origin);
    constantImageSource->SetSpacing(spacing);

    FloatImageType::DirectionType imageDirection;
    imageDirection.SetIdentity(); // no effect
    constantImageSource->SetDirection(imageDirection);
    constantImageSource->SetSize(size);
    constantImageSource->SetConstant(1.0);
    constantImageSource->UpdateOutputInformation();
    std::cout << "Canvas for projection image is ready to write" << std::endl;

    // 4) Prepare CT image to be projected
    auto ForwardProjection = forward_projector<DevFloatImageType>::type::New();

    itk::TimeProbe projProbe;
    std::cout << "Forward projection is now ongoing" << std::endl;

    ForwardProjection->SetInput(
        constantImageSource
            ->GetOutput()); // Canvas. projection image will be saved here.
    ForwardProjection->SetInput(1, spCTImg_mu); // reference plan CT image
    ForwardProjection->SetGeometry(spGeometry);
    projProbe.Start();
    TRY_AND_EXIT_ON_ITK_EXCEPTION(ForwardProjection->Update());
    projProbe.Stop();
    spResultProjImageFloat = ForwardProjection->GetOutput();
    std::cout << "Forward projection done in:	" << projProbe.GetMean() << ' '
              << projProbe.GetUnit() << '.' << std::endl;
  } // release all the memory

  auto stat_filter = itk::StatisticsImageFilter<FloatImageType>::New();
  stat_filter->SetInput(spResultProjImageFloat);
  stat_filter->Update();
  std::cerr << "Fwd Proj, mean: " << stat_filter->GetMean()
            << " min: " << stat_filter->GetMinimum()
            << " max: " << stat_filter->GetMaximum() << "\n";

  return spResultProjImageFloat;
  // From Float to USHORT and line integral to intensity
  /*
    spProjCT3D = UShortImageType::New(); // later
    const auto projCT_size = spResultProjImageFloat->GetLargestPossibleRegion()
                                 .GetSize(); // 1024 1024 350
    const auto projCT_idxStart =
        spResultProjImageFloat->GetLargestPossibleRegion().GetIndex(); // 0 0 0
    const auto projCT_spacing =
        spResultProjImageFloat->GetSpacing(); // 0.4 0.4 1.0
    const auto projCT_origin =
        spResultProjImageFloat->GetOrigin(); //-204.6 -204.6 -174.5

    // Copy informations from spResultProjImageFloat
    FloatImageType::RegionType projCT_region;
    projCT_region.SetSize(projCT_size);
    projCT_region.SetIndex(projCT_idxStart);

    spProjCT3D->SetRegions(projCT_region);
    spProjCT3D->SetSpacing(projCT_spacing);
    spProjCT3D->SetOrigin(projCT_origin);

    spProjCT3D->Allocate();
    spProjCT3D->FillBuffer(0);

    // Calculation process
    itk::ImageRegionConstIterator<FloatImageType> itSrc(
        spResultProjImageFloat, spResultProjImageFloat->GetRequestedRegion());
    itk::ImageRegionIterator<UShortImageType> itTarg(
        spProjCT3D, spProjCT3D->GetRequestedRegion()); // writing

    itSrc.GoToBegin();
    itTarg.GoToBegin();

    // Convert line integral to intensity value (I0/I = exp(mu_t)) --> I =
    // I0/exp(mu_t)
    const auto ushort_max = std::numeric_limits<unsigned short>::max();
    while (!itSrc.IsAtEnd() && !itTarg.IsAtEnd()) {
      const auto fProjVal = itSrc.Get();                  // mu_t //63.5
    --> 6.35 const auto tmpConvVal = ushort_max / exp(fProjVal); // physically
    true

      if (tmpConvVal <= 0.0) {
        itTarg.Set(0);
      } else if (tmpConvVal > ushort_max) {
        itTarg.Set(ushort_max);
      } else {
        itTarg.Set(static_cast<unsigned short>(tmpConvVal));
      }

      ++itSrc;
      ++itTarg;
    }

    // spProjCT3D: USHORT IMAGE of intensity. Not inverted (physical intensity)
    */
}

template <typename LikeImageType, typename OutImageType>
typename OutImageType::Pointer
create_empty_projections(typename LikeImageType::Pointer &spProjImg3D) {

  // Create a stack of empty projection images
  using ConstantImageSourceType =
      rtk::ConstantImageSource<OutImageType>; // Output: FLoat image = may be
                                              // mu_t = log(I_0/I)

  typename ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();

  typename ConstantImageSourceType::SizeType size{};
  typename ConstantImageSourceType::SpacingType spacing;
  typename ConstantImageSourceType::PointType origin;

  // std::cout << "Setting-up vacant projection image data" << std::endl;

  // a) size
  // std::cout << "chk1" << std::endl;
  size[0] = spProjImg3D->GetBufferedRegion().GetSize()[0];
  size[1] = spProjImg3D->GetBufferedRegion().GetSize()[1];
  size[2] = 1;

  // b) spacing
  spacing[0] = spProjImg3D->GetSpacing()[0];
  spacing[1] = spProjImg3D->GetSpacing()[1];
  spacing[2] = 1.0;

  // c) Origin: can center be the image center? or should be related to the CT
  // image???
  /*origin[0] = spacing[0] * (size[0] - 1) * -0.5;
  origin[1] = spacing[1] * (size[1] - 1) * -0.5;
  origin[2] = 0.0;*/

  origin[0] = spProjImg3D->GetOrigin()[0];
  origin[1] = spProjImg3D->GetOrigin()[1];
  origin[2] = 0.0;

  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);

  typename OutImageType::DirectionType imageDirection;
  imageDirection.SetIdentity(); // no effect
  constantImageSource->SetDirection(imageDirection);
  constantImageSource->SetSize(size);
  constantImageSource->SetConstant(1.0);
  constantImageSource->UpdateOutputInformation();

  return constantImageSource->GetOutput();
}

template <typename ImageType> auto fwd_method() {
  return std::string("CPU Joseph Fwd projection");
}
#if USE_CUDA
template <> inline auto fwd_method<CUDAFloatImageType>() {
  return std::string("CUDA Fwd projection");
}
#endif
// refer to YKPRoc later
// void YKPROC::ForwardProjection(FloatImageType::Pointer& spVolImgFloat, float
// fMVGanAngle, float panelOffsetX, float panelOffsetY, ,
// UShortImageType::Pointer& spProj3D)  iSliceIdx == Proj index 0 - 364
template <typename DevFloatImageType>
void CbctRecon::SingleForwardProjection(FloatImageType::Pointer &spVolImgFloat,
                                        const float fMVGanAngle,
                                        const float panelOffsetX,
                                        const float panelOffsetY,
                                        FloatImageType::Pointer &spProjImg3D,
                                        const int iSliceIdx) const {
  if (spVolImgFloat == nullptr) {
    return;
  }
  if (spProjImg3D == nullptr) {
    return;
  }

  // 2) Prepare empty projection images //Should be corresonponding to raw
  // projection images
  const int totalProjSize = spProjImg3D->GetBufferedRegion().GetSize()[2];
  if (iSliceIdx >= totalProjSize) {
    std::cout << "Error! totalProjSize= " << totalProjSize
              << " iSliceIdx= " << iSliceIdx << std::endl;
  }

  // std::cout << "Canvas for projection image is ready to write" << std::endl;

  // 4) Prepare CT image to be projected
  //    std::cout << "projection algorithm (0:Joseph, 1: CUDA, 2:RayCast ): " <<
  //    fwdMethod << std::endl;

  // Create forward projection image filter
  auto forward_projection = forward_projector<DevFloatImageType>::type::New();

  auto spGeometry = GeometryType::New();

  // 9 parameters are required
  const auto curSAD = 1000.0; // SourceToIsocenterDistances
  const auto curSDD = 1536.0;
  const double curGantryAngle = fMVGanAngle; // MV

  const double curProjOffsetX = panelOffsetX;
  const double curProjOffsetY = panelOffsetY;

  const auto curOutOfPlaneAngles = 0.0;
  const auto curInPlaneAngles = 0.0;

  const auto curSrcOffsetX = 0.0;
  const auto curSrcOffsetY = 0.0;

  spGeometry->AddProjection(
      curSAD, curSDD, curGantryAngle, curProjOffsetX, curProjOffsetY, // Flexmap
      curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
      curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0

  itk::TimeProbe projProbe;
  auto proj_husk = // Don't use const to allow for in-place filter
      create_empty_projections<FloatImageType, DevFloatImageType>(spProjImg3D);
  forward_projection->SetInput(
      proj_husk); // Canvas. projection image will be saved here.
  using Caster = itk::CastImageFilter<FloatImageType, DevFloatImageType>;
  auto caster = Caster::New();
  caster->SetInput(spVolImgFloat);
  caster->Update();
  forward_projection->SetInput(1,
                               caster->GetOutput()); // reference plan CT image
  forward_projection->SetGeometry(spGeometry);

  projProbe.Start();
  TRY_AND_EXIT_ON_ITK_EXCEPTION(forward_projection->Update())
  projProbe.Stop();

  const FloatImageType::Pointer resultFwdImg = forward_projection->GetOutput();

  std::cout << "Forward projection done by " << fwd_method<DevFloatImageType>()
            << " in: " << projProbe.GetMean() << ' ' << projProbe.GetUnit()
            << '.' << std::endl;

  auto stat_filter = itk::StatisticsImageFilter<FloatImageType>::New();
  stat_filter->SetInput(resultFwdImg);
  stat_filter->Update();
  std::cerr << "Fwd Proj, mean: " << stat_filter->GetMean()
            << " min: " << stat_filter->GetMinimum()
            << " max: " << stat_filter->GetMaximum() << "\n";

  // normalization or shift

  // typedef itk::MinimumMaximumImageCalculator<FloatImageType>
  // MinMaxCalculatorType;  MinMaxCalculatorType::Pointer spCalculator =
  // MinMaxCalculatorType::New();  spCalculator->SetImage(resultFwdImg);
  // spCalculator->Compute();

  // float minValAtt = spCalculator->GetMinimum();
  // float maxValAtt = spCalculator->GetMaximum();

  // float maxValProj = (65535.0 / exp(minValAtt));
  // float minValProj = (65535.0 / exp(maxValAtt)); //physically true

  // float valOffset = maxValProj - 65535.0; //not possible! always <=65535
  // if (valOffset < 0)
  //    valOffset = 0.0;

  // std::cout << "MaxValProj=" << maxValProj << " MInval=" << minValProj << "
  // ValOffset = " << valOffset << std::endl;

  itk::ImageRegionConstIterator<FloatImageType> itSrc(
      resultFwdImg, resultFwdImg->GetBufferedRegion()); // 2D

  // Convert line integral to intensity value (I0/I = exp(mu_t)) --> I =
  // I0/exp(mu_t)

  /*if (resultFwdImg->GetBufferedRegion().GetSize()[0] != pYKImage2D->m_iWidth)
      return;*/

  itSrc.GoToBegin();

  itk::ImageSliceIteratorWithIndex<FloatImageType> it_FwdProj3D(
      spProjImg3D, spProjImg3D->GetBufferedRegion());

  it_FwdProj3D.SetFirstDirection(0);
  it_FwdProj3D.SetSecondDirection(1);
  it_FwdProj3D.GoToBegin();

  auto curSliceIdx = 0;

  while (!it_FwdProj3D.IsAtEnd()) {
    if (curSliceIdx == iSliceIdx) {
      // Search matching slice using slice iterator for m_spProjCTImg
      while (!it_FwdProj3D.IsAtEndOfSlice() && !itSrc.IsAtEnd()) {
        while (!it_FwdProj3D.IsAtEndOfLine() && !itSrc.IsAtEnd()) {
          const auto fProjVal = itSrc.Get();
          // mu_t, the lower means the higher attn.
          // mu_t //63.5 --> 6.35
          const auto tmpConvVal = 65535.0 / exp(fProjVal); // intensity value

          /*
          unsigned short val = 0;
          if (tmpConvVal <= 0.0) {
            val = 0;
          } else if (tmpConvVal > 65535.0) {
            val = 65535;
          } else {
            val = static_cast<unsigned short>(tmpConvVal);
          }

          // unsigned short tmpVal = (unsigned short)(it_FwdProj3D.Get());
          // tmpVal = 65535 - tmpVal; //inverse is done here

          it_FwdProj3D.Set(val);
          */
          it_FwdProj3D.Set(tmpConvVal);

          ++it_FwdProj3D;
          ++itSrc;
        }
        it_FwdProj3D.NextLine();
      }
      it_FwdProj3D.NextSlice();
    }
    it_FwdProj3D.NextSlice();
    curSliceIdx++;
  }
  // Save this file

} // release all the memory

#endif
