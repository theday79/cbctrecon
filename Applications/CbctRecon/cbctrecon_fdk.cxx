/* Do reconstruction algorithms */
#include "cbctrecon.h"

// std
#include <math.h>                                           // for asin, cos
#include <stdlib.h>
#include <iostream>              // for operator<<, endl, basic_ostream, cout
#include <memory>                // for unique_ptr
#include <algorithm>                                        // for copy, min
#include <complex>                                          // for operator*
#include <string>                                           // for operator<<
#include <utility>                                          // for swap

// Qt
#include <QDir>
#include <QCheckBox>           // for QCheckBox
#include <QtGlobal>             // for qRound
#include <QLineEdit>           // for QLineEdit
#include <QSpinBox>            // for QSpinBox
#include <QString>             // for QString
#include <QByteArray>                                     // for QByteArray
#include <QFileInfo>

// ITK
#include "itkArray.h"                                     // for Array::Ar...
#include "itkArray2D.h"                                   // for Array2D::...
#include "itkBinaryFunctorImageFilter.h"                  // for BinaryFun...
#include "itkBoxImageFilter.h"                            // for BoxImageF...
#include "itkConstNeighborhoodIterator.h"                 // for ConstNeig...
#include "itkDataObjectDecorator.h"                       // for DataObjec...
#include "itkDiffusionTensor3D.h"                         // for Diffusion...
#include "itkEuler3DTransform.h"                          // for Euler3DTr...
#include "itkExtractImageFilter.h"                        // for ExtractIm...
#include "itkFixedArray.h"                                // for FixedArra...
#include "itkHalfHermitianToRealInverseFFTImageFilter.h"  // for HalfHermi...
#include "itkImage.h"                                     // for Image::Al...
#include "itkImageBase.h"                                 // for ImageBase...
#include "itkImageConstIteratorWithIndex.h"               // for ImageCons...
#include "itkImageDuplicator.h"                           // for ImageDupl...
#include "itkImageFileWriter.h"                           // for ImageFile...
#include "itkImageFunction.h"                             // for ImageFunc...
#include "itkImageRegion.h"                               // for ImageRegi...
#include "itkImageRegionConstIteratorWithIndex.h"         // for ImageRegi...
#include "itkImageRegionIterator.h"                       // for ImageRegi...
#include "itkImageRegionIteratorWithIndex.h"              // for ImageRegi...
#include "itkImageScanlineIterator.h"                     // for ImageScan...
#include "itkImageSource.h"                               // for ImageSour...
#include "itkImageToImageFilter.h"                        // for ImageToIm...
#include "itkImportImageContainer.h"                      // for ImportIma...
#include "itkInPlaceImageFilter.h"                        // for InPlaceIm...
#include "itkMatrix.h"                                    // for Matrix::o...
#include "itkMatrixOffsetTransformBase.h"                 // for MatrixOff...
#include "itkNeighborhood.h"                              // for Neighborh...
#include "itkNeighborhoodAlgorithm.h"                     // for ImageBoun...
#include "itkOptimizerParameters.h"                       // for Optimizer...
#include "itkPoint.h"                                     // for operator<<
#include "itkRealToHalfHermitianForwardFFTImageFilter.h"  // for RealToHal...
#include "itkResourceProbe.h"                             // for ResourceP...
#include "itkRigid3DTransform.h"                          // for Rigid3DTr...
#include "itkSimpleDataObjectDecorator.h"                 // for SimpleDat...
#include "itkSymmetricEigenAnalysis.h"                    // for Symmetric...
#include "itkSymmetricSecondRankTensor.h"                 // for Symmetric...
#include "itkThresholdImageFilter.h"                      // for Threshold...
#include "itkTransform.h"                                 // for Transform...
#include "itkUnaryFunctorImageFilter.h"                   // for UnaryFunc...
#include "itkVariableLengthVector.h"                      // for VariableL...
#include "itkVector.h"                                    // for operator<<
#include "itkAbsImageFilter.h"
#include "itkCastImageFilter.h"                             // for CastImage...
#include "itkEuler3DTransform.h"                          // for Euler3DTr...
#include "itkExtrapolateImageFunction.h"          // for ExtrapolateImageFun...
#include "itkFixedArray.h"                                // for FixedArra...
#include "itkFlipImageFilter.h"
#include "itkIdentityTransform.h"                           // for IdentityT...
#include "itkImage.h"                                     // for Image::Al...
#include "itkImageDuplicator.h"                           // for ImageDupl...
#include "itkImageFileWriter.h"                           // for ImageFile...
#include "itkImageIORegion.h"                               // for operator<<
#include "itkImageToImageFilter.h"                        // for ImageToIm...
#include "itkIndex.h"                                       // for operator<<
#include "itkInterpolateImageFunction.h"          // for InterpolateImageFun...
#include "itkLightObject.h"                 // for LightObject::P...
#include "itkLinearInterpolateImageFunction.h"    // for LinearInterpolateIm...
#include "itkMath.h"                        // for pi
#include "itkMatrix.h"                                    // for Matrix::o...
#include "itkMatrixOffsetTransformBase.h"                 // for MatrixOff...
#include "itkMedianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkOffset.h"                                      // for operator<<
#include "itkProcessObject.h"               // for ProcessObject:...
#include "itkRegion.h"                      // for Region::Region...
#include "itkResampleImageFilter.h"
#include "itkSize.h"                                        // for operator<<
#include "itkSmartPointer.h"                // for SmartPointer
#include "itkStreamingImageFilter.h"
#include "itkThresholdImageFilter.h"                      // for Threshold...
#include "itkTimeProbe.h"
#include "itkTransform.h"                                 // for Transform...
#include "itkZeroFluxNeumannBoundaryCondition.h"            // for ZeroFluxN...
#include "vnl_matrix_fixed.h"                             // for vnl_matri...
#include "vnl_vector.h"                                     // for operator*

// RTK
#include "rtkMacro.h"            // for TRY_AND_EXIT_ON_ITK_EXCEPTION
#include "rtkConstantImageSource.h"
#include "rtkFDKConeBeamReconstructionFilter.h"   // for FDKConeBeamReconstr...
#include "rtkFFTRampImageFilter.h"                // for FFTRampImageFilter
#include "rtkThreeDCircularProjectionGeometry.h"  // for ThreeDCircularProje...
#include "rtkParkerShortScanImageFilter.h"
#include "rtkDisplacedDetectorImageFilter.h"
#include "rtkBackProjectionImageFilter.h"                 // for BackProje...
#include "rtkFFTProjectionsConvolutionImageFilter.h"      // for FFTProjec...
#include "rtkFieldOfViewImageFilter.h"

#if USE_CUDA
#include <rtkCudaDisplacedDetectorImageFilter.h>
#include <rtkCudaFDKConeBeamReconstructionFilter.h>
#include <rtkCudaParkerShortScanImageFilter.h>
#endif // USE_CUDA

// PLM
#if USE_OPENCL_PLM
#include <autotune_opencl.h>
#include <fdk.h>
#include <fdk_opencl.h>
#include <opencl_util.h>
#include <plm_image.h>
#include <plmreconstruct_config.h>
#include <proj_image.h>
#include <proj_image_filter.h>
#include <proj_matrix.h>
#endif // USE_OPENCL_PLM

#if USE_OPENCL_RTK
#include "rtkOpenCLFDKConeBeamReconstructionFilter.h"
#endif

// Local
#include "cbctrecon_compute.h"
#include "YK16GrayImage.h"       // for YK16GrayImage
#include "ui_cbctrecon.h"        // for CbctReconClass

#if USE_CUDA
using CUDAFloatImageType = itk::CudaImage<float, 3U>;

void CbctRecon::CudaDoReconstructionFDK(enREGI_IMAGES target) {
  using castToCudaType =
      itk::CastImageFilter<FloatImageType, CUDAFloatImageType>;
  castToCudaType::Pointer castToCuda = castToCudaType::New();
  castToCuda->SetInput(m_spProjImg3DFloat);
  castToCuda->Update();
  CUDAFloatImageType::Pointer cuda_spProjImg3DFloat = castToCuda->GetOutput();

  if (cuda_spProjImg3DFloat == nullptr) {
    std::cout << "processed Projection image is not ready yet" << std::endl;
    return;
  }

  std::cout << "CUDA method will be used..." << std::endl;

  using DuplicatorType = itk::ImageDuplicator<CUDAFloatImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(cuda_spProjImg3DFloat);
  duplicator->Update();

  CUDAFloatImageType::Pointer spCurImg =
      duplicator->GetOutput(); // already down sampled

  // Displaced detector weighting // set pipeline //inplace filter
  using DDFType = rtk::CudaDisplacedDetectorImageFilter;
  DDFType::Pointer ddf = DDFType::New();
  if (ui.checkBox_UseDDF->isChecked()) {
    ddf->SetInput(spCurImg);
    ddf->SetGeometry(m_spCustomGeometry);
    std::cout << "DDF was set in pipeline" << std::endl;

    if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
      ddf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = ddf->GetOutput();
  }

  using PSSFType = rtk::CudaParkerShortScanImageFilter;
  PSSFType::Pointer pssf = PSSFType::New();
  if (ui.checkBox_UsePSSF->isChecked()) {
    // Short scan image filter
    pssf->SetInput(spCurImg);
    // pssf->SetGeometry( geometryReader->GetOutputObject() );
    pssf->SetGeometry(m_spCustomGeometry);
    // pssf->InPlaceOff(); //YKComments: Do not overwrite input image buffer for
    // output

    pssf->GlobalWarningDisplayOff(); // I don't care about any of the potential
                                     // warnings

    if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
      pssf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = pssf->GetOutput();
  }

  // Just Before going to FDK recon,
  // Update Projection data and delete old data.

  // Let's duplicate this
  // Original m_spProjImg3D will be deleted after update

  if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
    using DuplicatorType = itk::ImageDuplicator<CUDAFloatImageType>;
    DuplicatorType::Pointer ImDuplicator = DuplicatorType::New();
    ImDuplicator->SetInputImage(spCurImg);
    ImDuplicator->Update();

    using CastFilterType =
        itk::CastImageFilter<CUDAFloatImageType, FloatImageType>;
    CastFilterType::Pointer CastFilter = CastFilterType::New();
    CastFilter->SetInput(ImDuplicator->GetOutput());
    CastFilter->Update();

    m_spProjImg3DFloat = CastFilter->GetOutput();

    SetMaxAndMinValueOfProjectionImage(); // scan m_spProjImg3D and update
                                          // m_fProjImgValueMin, max
    SLT_DrawProjImages();
  }

  // Generate image sources for cone beam CT reconstruction
  using ConstantImageSourceType = rtk::ConstantImageSource<CUDAFloatImageType>;

  ConstantImageSourceType::SizeType sizeOutput{};
  sizeOutput[0] = ui.lineEdit_outImgDim_AP->text().toInt(); // pixel
  sizeOutput[1] =
      ui.lineEdit_outImgDim_SI->text()
          .toInt(); // Caution!: direction is different in NKI SCAN FIle
  sizeOutput[2] = ui.lineEdit_outImgDim_LR->text().toInt();

  ConstantImageSourceType::SpacingType spacing;
  spacing[0] = ui.lineEdit_outImgSp_AP->text().toDouble();
  spacing[1] = ui.lineEdit_outImgSp_SI->text().toDouble();
  spacing[2] = ui.lineEdit_outImgSp_LR->text().toDouble();
  if (GetOutputResolutionFromFOV<ConstantImageSourceType, FloatImageType>(
          sizeOutput, spacing, m_spCustomGeometry, m_spProjImg3DFloat,
          ui.lineEdit_OutputFilePath->text())) {
    std::cout << "Reconstruction resolution and image size were set "
                 "automatically, as no outputfilepath was given."
              << std::endl;
  }

  ConstantImageSourceType::PointType origin;
  origin[0] = -0.5 * sizeOutput[0] * spacing[0]; // Y in DCM?
  origin[1] = -0.5 * sizeOutput[1] * spacing[1]; // Z in DCM?
  origin[2] = -0.5 * sizeOutput[2] * spacing[2]; // X in DCM?

  ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();
  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);
  constantImageSource->SetSize(sizeOutput);
  constantImageSource->SetConstant(0.0); // initial value

  double fTruncCorFactor =
      ui.lineEdit_Ramp_TruncationCorrection->text().toDouble();
  const double fHannCut = ui.lineEdit_Ramp_HannCut->text().toDouble();
  const double fCosineCut = ui.lineEdit_Ramp_CosineCut->text().toDouble();
  const double fHamming = ui.lineEdit_Ramp_Hamming->text().toDouble();
  const double fHannCutY = ui.lineEdit_Ramp_HannCutY->text().toDouble();

  if (fTruncCorFactor > 0.0 && target == REGISTER_COR_CBCT) {
    std::cout << "Warning! Truncation factor is " << fTruncCorFactor
              << ". Regardless of previous setting, this factor should not be "
                 "0 for scatter corrected CBCT. Now zero value is applied."
              << std::endl;
    fTruncCorFactor = 0.0;
  }

  using StreamerType =
      itk::StreamingImageFilter<CUDAFloatImageType, CUDAFloatImageType>;
  // FDK reconstruction filtering
  using FDKCUDAType = rtk::CudaFDKConeBeamReconstructionFilter;
  FDKCUDAType::Pointer CUDAfeldkamp = FDKCUDAType::New();

  std::cout << "CUDA will be used for FDK reconstruction" << std::endl;
  CUDAfeldkamp->SetInput(0, constantImageSource->GetOutput());
  CUDAfeldkamp->SetInput(1, spCurImg);
  CUDAfeldkamp->SetGeometry(m_spCustomGeometry);
  CUDAfeldkamp->GetRampFilter()->SetTruncationCorrection(fTruncCorFactor);
  CUDAfeldkamp->GetRampFilter()->SetHannCutFrequency(fHannCut);
  CUDAfeldkamp->GetRampFilter()->SetCosineCutFrequency(fCosineCut);
  CUDAfeldkamp->GetRampFilter()->SetHammingFrequency(fHamming);
  CUDAfeldkamp->GetRampFilter()->SetHannCutFrequencyY(fHannCutY);

  if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
    CUDAfeldkamp->Update();
  } else {
    CUDAfeldkamp->UpdateOutputInformation();
  }

  std::cout << "Cone beam reconstruction pipeline is ready" << std::endl;

  // Streaming depending on streaming capability of writer --> not affect the
  // calc. speed

  using StreamerType =
      itk::StreamingImageFilter<CUDAFloatImageType, CUDAFloatImageType>;
  StreamerType::Pointer streamerBP = StreamerType::New();

  streamerBP->SetInput(CUDAfeldkamp->GetOutput());
  streamerBP->SetNumberOfStreamDivisions(4); // YK: 1 in example code from
                                             // "rtkfdk" //AG: stated in test: 4
                                             // for ITK MAJOR >= 4

  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;

  using FOVfilterType =
      rtk::FieldOfViewImageFilter<CUDAFloatImageType, CUDAFloatImageType>;
  typename FOVfilterType::Pointer FOVfilter = FOVfilterType::New();
  FOVfilter->SetGeometry(m_spCustomGeometry);
  FOVfilter->SetInput(0, streamerBP->GetOutput());
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
  CUDAFloatImageType::Pointer targetImg =
      FOVfilter->GetOutput(); // streamerBP->GetOutput();

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

  itk::Vector<double, 3U> offset(0.0);

  using TransformType = itk::MatrixOffsetTransformBase<double, 3U, 3U>;
  TransformType::Pointer transform = TransformType::New();
  transform->SetMatrix(CoordChangeMatrix);
  transform->SetOffset(offset);

  using ResampleFilterType =
      itk::ResampleImageFilter<CUDAFloatImageType, CUDAFloatImageType>;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetInput(targetImg);
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(targetImg->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  // resampler->Update();//yktemp Error 2

  // LR flip

  std::cout << "Flip filter is being applied" << std::endl;
  using FilterType = itk::FlipImageFilter<CUDAFloatImageType>;
  FilterType::Pointer flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;
  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;
  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  using AbsImageFilterType =
      itk::AbsImageFilter<CUDAFloatImageType, CUDAFloatImageType>;
  AbsImageFilterType::Pointer absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(
      flipFilter
          ->GetOutput()); // 20140206 modified it was a bug
                          // absImgFilter->SetInput(resampler->GetOutput());

  using MultiplyImageFilterType = itk::MultiplyImageFilter<CUDAFloatImageType>;
  MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  using CastFilterType =
      itk::CastImageFilter<CUDAFloatImageType, UShortImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());
  try {
    castFilter->Update(); // YK20150109
  } catch (const std::exception &e) {
    std::cerr << "std::exception thrown: " << e.what() << std::endl;
    return;
  }

  UShortImageType::Pointer tmpReconImg;
  // if all 0 0 0 don't do the median filtering

  itk::TimeProbe reconTimeProbe;
  reconTimeProbe.Start();

  std::cout << "Reconstructing the image.. please wait..." << std::endl;
  UShortImageType::SizeType indexRadius{};
  indexRadius[0] = ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
  indexRadius[1] = ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
  indexRadius[2] = ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along y
  if (ui.checkBox_PostMedianOn->isChecked() &&
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

  reconTimeProbe.Stop();
  std::cout << "It took " << reconTimeProbe.GetMean() << ' '
            << reconTimeProbe.GetUnit() << std::endl;
  ui.lineEdit_ReconstructionTime->setText(
      QString("%1").arg(reconTimeProbe.GetMean()));

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
  QString outputFilePath = ui.lineEdit_OutputFilePath->text();

  QFileInfo outFileInfo(outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
    ui.lineEdit_Cur3DFileName->setText("FDK-reconstructed volume");
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFilePath.toLocal8Bit().constData());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(m_spCrntReconImg);

    std::cout << "Writing the image to: "
              << outputFilePath.toLocal8Bit().constData() << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

    ui.lineEdit_Cur3DFileName->setText(outputFilePath);
    std::cout << std::endl;
    std::cout << "Output generation was succeeded" << std::endl;
  }

  m_dspYKReconImage->CreateImage(size_trans[0], size_trans[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(size_trans[2] - 1);
  ui.spinBoxReconImgSliceNo->setValue(
      qRound(size_trans[2] / 2.0)); // DrawReconImage is called automatically


  using CastBackType = itk::CastImageFilter<CUDAFloatImageType, FloatImageType>;
  CastBackType::Pointer castBack = CastBackType::New();
  castBack->SetInput(cuda_spProjImg3DFloat);
  castBack->Update();
  m_spProjImg3DFloat = castBack->GetOutput();

  std::cout << "FINISHED!: FDK CBCT reconstruction" << std::endl;
}
#else
void CbctRecon::CudaDoReconstructionFDK(enREGI_IMAGES target) {
  std::cout << "This program were not compiled with the CUDA option, please "
               "select OpenCL or CPU instead!"
            << std::endl;
  return;
}
#endif

#if USE_OPENCL_PLM
FloatImageType::Pointer PlastimatchOpenCLFDK(
    const FloatImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    std::array<float, 3> spacing, std::array<plm_long, 3> sizeOutput) {
  spCurImg->Update();

  //////////////////////// start PLM fdk /////////////////////////////
  // Copy of opencl_reconstruct_conebeam from fdk_opencl.cxx in plastimatch
  // (because we already have loaded projections into memory)
  Opencl_device ocl_dev{};
  Opencl_buf *ocl_buf_vol;
  Opencl_buf *ocl_buf_img;
  Opencl_buf *ocl_buf_matrix;
  cl_int4 ocl_vol_dim{};
  Opencl_buf *ocl_buf_vol_origin;
  Opencl_buf *ocl_buf_vol_spacing;
  cl_int2 ocl_proj_dim{};
  Opencl_buf *ocl_buf_nrm;
  Opencl_buf *ocl_buf_ic;
  // Proj_image *proj;
  // int image_num;
  // float scale;
  std::cout << "Loading plmopencl..." << std::endl;
  LOAD_LIBRARY_SAFE(libplmopencl);
  LOAD_SYMBOL(opencl_open_device, libplmopencl);
  LOAD_SYMBOL(opencl_load_programs, libplmopencl);
  LOAD_SYMBOL(opencl_kernel_create, libplmopencl);
  LOAD_SYMBOL(opencl_buf_create, libplmopencl);
  LOAD_SYMBOL(opencl_buf_write, libplmopencl);
  LOAD_SYMBOL(opencl_set_kernel_args, libplmopencl);
  LOAD_SYMBOL(opencl_kernel_enqueue, libplmopencl);
  LOAD_SYMBOL(opencl_buf_read, libplmopencl);

  std::cout << "Start opencl environment and enqueue kernels..." << std::endl;
  /* Set up devices and kernels */
  opencl_open_device(&ocl_dev);
  opencl_load_programs(&ocl_dev, "fdk_opencl.cl");
  opencl_kernel_create(&ocl_dev, "fdk_kernel_nn");

  /* Retrieve 2D image to get dimensions */
  const itk::Size<3U> proj_dim = spCurImg->GetLargestPossibleRegion().GetSize();
  // proj = proj_dir->load_image(0);

  // Generate image sources for cone beam CT reconstruction
  std::array<float, 3> offset = {
      {static_cast<float>(-0.5 * sizeOutput[0] * spacing[0]),
       static_cast<float>(-0.5 * sizeOutput[1] * spacing[1]),
       static_cast<float>(-0.5 * sizeOutput[2] * spacing[2])}};

  auto *vol =
      new Volume(&sizeOutput[0], &offset[0], &spacing[0], nullptr, PT_FLOAT, 1);
  std::cout << "Buffer vol: " << vol->pix_size * vol->npix << std::endl;
  std::cout << "Buffer proj: " << proj_dim[1] * proj_dim[0] * sizeof(float)
            << std::endl;
  /* Set up device memory */
  ocl_buf_vol =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                        vol->pix_size * vol->npix, vol->img);

  ocl_buf_img =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        proj_dim[1] * proj_dim[0] * sizeof(float), // is 0
                        nullptr);

  ocl_buf_matrix =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        12 * sizeof(float), nullptr);

  /* Copy volume dim (convert from size_t to int) */
  ocl_vol_dim.x = vol->dim[0];
  ocl_vol_dim.y = vol->dim[1];
  ocl_vol_dim.z = vol->dim[2];

  ocl_buf_vol_origin =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                        3 * sizeof(float), &vol->origin[0]);

  ocl_buf_vol_spacing =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                        3 * sizeof(float), &vol->spacing[0]);

  /* Copy projection image dim (convert from size_t to int) */
  ocl_proj_dim.x = proj_dim[0];
  ocl_proj_dim.y = proj_dim[1];

  ocl_buf_nrm =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        3 * sizeof(float), nullptr);

  ocl_buf_ic =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        2 * sizeof(float), nullptr);

  /* Calculate the scale */
  const float scale =
      static_cast<float>(sqrt(3.0) / static_cast<double>(1 + proj_dim[2])) *
      10.0f;

  std::vector<double>::const_iterator itShiftX, itShiftY;
  itShiftX = m_spCustomGeometry->GetProjectionOffsetsX().begin();
  itShiftY = m_spCustomGeometry->GetProjectionOffsetsY().begin();
  std::cout << "Proj. # " << std::endl;
  /* Project each image into the volume one at a time */
  for (int image_num = 0; image_num < static_cast<int>(proj_dim[2]);
       image_num++) {
    /* Translate image and properties to PLM */
    auto *proj = new Proj_image;
    proj->init();
    proj->dim[0] = proj_dim[0];
    proj->dim[1] = proj_dim[1];
    proj->xy_offset[0] = *itShiftX;
    proj->xy_offset[1] = *itShiftY;
    itk::Index<3U> cur_idx = {{0, 0, image_num}};
    proj->img = &spCurImg->GetPixel(cur_idx);
    if (proj->img == ITK_NULLPTR) {
      std::cout << "No image!" << std::endl;
    }

    proj->pmat = new Proj_matrix;
    proj->pmat->ic[0] = 0.5 * proj_dim[0] - 0.5 + proj->xy_offset[0];
    proj->pmat->ic[1] = 0.5 * proj_dim[1] - 0.5 + proj->xy_offset[1];
    proj->pmat->sad =
        m_spCustomGeometry->GetSourceToIsocenterDistances()[image_num];

    const double sid =
        m_spCustomGeometry->GetSourceToDetectorDistances()[image_num];
    const double angle = m_spCustomGeometry->GetGantryAngles()[image_num];
    const double tgt[] = {0.0, 0.0, 0.0};
    const double cam[] = {tgt[0] + proj->pmat->sad * cos(angle),
                          tgt[1] - proj->pmat->sad * sin(angle), tgt[2]};
    const double vup[] = {0, 0, 1}; // R U sure about that?
    const double ps[] = {spCurImg->GetSpacing()[0], spCurImg->GetSpacing()[1]};

    proj->pmat->set(&cam[0], &tgt[0], &vup[0], sid, &proj->pmat->ic[0],
                    &ps[0]); // , proj->dim);

    if (image_num % (proj_dim[2] / 10) == 0) {
      std::cout << image_num << "::" << angle * 180 * itk::Math::one_over_pi
                << std::endl;
      proj->stats();
      std::cout << "Cam::" << proj->pmat->cam[0] << ", " << proj->pmat->cam[1]
                << ", " << proj->pmat->cam[2] << std::endl;
    }

    /* Apply ramp filter */
    proj_image_filter(proj);

    /* Copy image bytes to device */
    opencl_buf_write(&ocl_dev, ocl_buf_img,
                     proj_dim[1] * proj_dim[0] * sizeof(float), proj->img);

    /* Copy matrix to device (convert from double to float) */
    std::array<float, 12> matrix{};
    std::memcpy(matrix.data(), &proj->pmat->matrix[0], 12);

    opencl_buf_write(&ocl_dev, ocl_buf_matrix, 12 * sizeof(float), &matrix[0]);

    /* Copy ic to device (convert from double to float) */
    float ic[] = {static_cast<float>(proj->pmat->ic[0]),
                  static_cast<float>(proj->pmat->ic[1])};
    opencl_buf_write(&ocl_dev, ocl_buf_ic, 2 * sizeof(float), &ic[0]);

    /* Copy nrm to device (convert from double to float) */
    float nrm[] = {static_cast<float>(proj->pmat->nrm[0]),
                   static_cast<float>(proj->pmat->nrm[1]),
                   static_cast<float>(proj->pmat->nrm[2])};
    opencl_buf_write(&ocl_dev, ocl_buf_nrm, 3 * sizeof(float), &nrm[0]);

    /* Convert sad from double to float */
    const cl_float sad = proj->pmat->sad;

    /* Set fdk kernel arguments */
    opencl_set_kernel_args( // Marked as "dangerous" c-style function by
                            // cppcoreguidelines-pro-type-vararg
        &ocl_dev, sizeof(cl_mem), &ocl_buf_vol[0], sizeof(cl_mem),
        &ocl_buf_img[0], sizeof(cl_mem), &ocl_buf_matrix[0], sizeof(cl_int4),
        &ocl_vol_dim, sizeof(cl_mem), &ocl_buf_vol_origin[0], sizeof(cl_mem),
        &ocl_buf_vol_spacing[0], sizeof(cl_int2), &ocl_proj_dim, sizeof(cl_mem),
        &ocl_buf_nrm[0], sizeof(cl_mem), &ocl_buf_ic[0], sizeof(cl_float), &sad,
        sizeof(cl_float), &scale, static_cast<size_t>(0));

    /* Compute workgroup size */
    /* (Max local_work_size for (Greg's?) ATI RV710 is 128) */
    size_t local_work_size = 128;
    size_t global_work_size = static_cast<float>(vol->npix);

    /* Invoke kernel */
    opencl_kernel_enqueue(&ocl_dev, global_work_size, local_work_size);
    itShiftX++;
    itShiftY++;

    delete proj->pmat;
    delete proj;
  }
  std::cout << std::endl;
  std::cout << "Reading results..." << std::endl;
  /* Read back results */
  opencl_buf_read(&ocl_dev, ocl_buf_vol, vol->pix_size * vol->npix, vol->img);

  UNLOAD_LIBRARY(libplmopencl);
  return Plm_image(vol).itk_float();
  // End of opencl_reconstruct_conebeam
  //////////////////////// end PLM fdk /////////////////////////////
}
#endif

#if USE_OPENCL_RTK
FloatImageType::Pointer RTKOpenCLFDK(
    const FloatImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    FloatImageType::SpacingType spacing, FloatImageType::SizeType sizeOutput,
    std::array<const double, 5> fdk_options) {

  try {
    spCurImg->Update();
  } catch (const std::exception &err) {
    std::cerr << "Couldn't update spCurImg: " << err.what() << std::endl;
  }
  // Generate image sources for cone beam CT reconstruction
  using ConstantImageSourceType = rtk::ConstantImageSource<FloatImageType>;
  ConstantImageSourceType::PointType origin;

  origin[0] = -0.5 * sizeOutput[0] * spacing[0]; // Y in DCM?
  origin[1] = -0.5 * sizeOutput[1] * spacing[1]; // Z in DCM?
  origin[2] = -0.5 * sizeOutput[2] * spacing[2]; // X in DCM?

  ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();
  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);
  constantImageSource->SetSize(sizeOutput);
  constantImageSource->SetConstant(0.0); // initial value
  // constantImageSource->Update();

  // just for expressiveness
  const double fTruncCorFactor = fdk_options.at(0);
  const double fHannCut = fdk_options.at(1);
  const double fCosineCut = fdk_options.at(2);
  const double fHamming = fdk_options.at(3);
  const double fHannCutY = fdk_options.at(4);

  // FDK reconstruction filtering
  using FDKOPENCLType = rtk::OpenCLFDKConeBeamReconstructionFilter;
  FDKOPENCLType::Pointer feldkampOCL;
  feldkampOCL = FDKOPENCLType::New();

  feldkampOCL->SetInput(0, constantImageSource->GetOutput());
  feldkampOCL->SetInput(1, spCurImg);
  feldkampOCL->SetGeometry(m_spCustomGeometry);
  feldkampOCL->GetRampFilter()->SetTruncationCorrection(fTruncCorFactor);
  feldkampOCL->GetRampFilter()->SetHannCutFrequency(fHannCut);
  feldkampOCL->GetRampFilter()->SetHannCutFrequencyY(fHannCutY);
  feldkampOCL->GetRampFilter()->SetCosineCutFrequency(fCosineCut);
  feldkampOCL->GetRampFilter()->SetHammingFrequency(fHamming);

  feldkampOCL->Update();
  // feldkampOCL->PrintTiming(std::cout); Deprecated in rtk 1.4

  return feldkampOCL->GetOutput();
}
#endif // USE_OPENCL_RTK

#if (USE_OPENCL_PLM || USE_OPENCL_RTK)
// Using plastimatch and not RTK
void CbctRecon::OpenCLDoReconstructionFDK(enREGI_IMAGES target) {
  if (m_spProjImg3DFloat == nullptr) {
    std::cout << "processed Projection image is not ready yet" << std::endl;
    return;
  }

  using DuplicatorType = itk::ImageDuplicator<FloatImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(m_spProjImg3DFloat);
  duplicator->Update();

  FloatImageType::Pointer spCurImg =
      duplicator->GetOutput(); // already down sampled
  using DDFType = rtk::DisplacedDetectorImageFilter<FloatImageType>;
  DDFType::Pointer ddf = DDFType::New();

  if (ui.checkBox_UseDDF->isChecked()) {
    ddf->SetInput(spCurImg);
    ddf->SetGeometry(m_spCustomGeometry);
    std::cout << "DDF was set in pipeline" << std::endl;

    if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
      ddf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = ddf->GetOutput();
  }

  using PSSFType = rtk::ParkerShortScanImageFilter<FloatImageType>;
  PSSFType::Pointer pssf = PSSFType::New();

  if (ui.checkBox_UsePSSF->isChecked()) {
    // Short scan image filter
    pssf->SetInput(spCurImg);
    // pssf->SetGeometry( geometryReader->GetOutputObject() );
    pssf->SetGeometry(m_spCustomGeometry);
    // pssf->InPlaceOff(); //YKComments: Do not overwrite input image buffer for
    // output
    std::cout << "short scan image filter success" << std::endl;

    if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
      pssf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = pssf->GetOutput();
  }
  // Just Before going to FDK recon,
  // Update Projection data and delete old data.

  // Let's duplicate this
  // Original m_spProjImg3D will be deleted after update

  if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
    using DuplicatorType = itk::ImageDuplicator<FloatImageType>;
    DuplicatorType::Pointer ImDuplicator = DuplicatorType::New();
    ImDuplicator->SetInputImage(spCurImg);
    ImDuplicator->Update();
    m_spProjImg3DFloat = ImDuplicator->GetOutput();

    SetMaxAndMinValueOfProjectionImage(); // scan m_spProjImg3D and update
                                          // m_fProjImgValueMin, max
    SLT_DrawProjImages();
  }

#if DISABLE_USE_OPENCL_PLM
  float spacing[3];
  plm_long sizeOutput[3];

  sizeOutput[0] = ui.lineEdit_outImgDim_AP->text().toInt(); // pixel
  sizeOutput[1] =
      ui.lineEdit_outImgDim_SI->text()
          .toInt(); // Caution!: direction is different in NKI SCAN FIle
  sizeOutput[2] = ui.lineEdit_outImgDim_LR->text().toInt();

  spacing[0] = ui.lineEdit_outImgSp_AP->text().toDouble();
  spacing[1] = ui.lineEdit_outImgSp_SI->text().toDouble();
  spacing[2] = ui.lineEdit_outImgSp_LR->text().toDouble();
  FloatImageType::Pointer targetImg =
      PlastimatchOpenCLFDK(spCurImg, m_spCustomGeometry, spacing,
                           sizeOutput); // streamerBP->GetOutput();
#endif

#if USE_OPENCL_RTK
  FloatImageType::SpacingType spacing;
  FloatImageType::SizeType sizeOutput{};
  sizeOutput[0] = ui.lineEdit_outImgDim_AP->text().toInt(); // pixel
  sizeOutput[1] =
      ui.lineEdit_outImgDim_SI->text()
          .toInt(); // Caution!: direction is different in NKI SCAN FIle
  sizeOutput[2] = ui.lineEdit_outImgDim_LR->text().toInt();

  spacing[0] = ui.lineEdit_outImgSp_AP->text().toDouble();
  spacing[1] = ui.lineEdit_outImgSp_SI->text().toDouble();
  spacing[2] = ui.lineEdit_outImgSp_LR->text().toDouble();

  if (GetOutputResolutionFromFOV<FloatImageType, FloatImageType>(
          sizeOutput, spacing, m_spCustomGeometry, m_spProjImg3DFloat,
          ui.lineEdit_OutputFilePath->text())) {
    std::cout << "Reconstruction resolution and image size were set "
                 "automatically, as no outputfilepath was given."
              << std::endl;
  }

  double fTruncCorFactor =
      ui.lineEdit_Ramp_TruncationCorrection->text().toDouble();
  if (fTruncCorFactor > 0.0 && target == REGISTER_COR_CBCT) {
    std::cout << "Warning! Truncation factor is " << fTruncCorFactor
              << ". Regardless of previous setting, this factor should not be "
                 "0 for scatter corrected CBCT. Now zero value is applied."
              << std::endl;
    fTruncCorFactor = 0.0;
  }

  std::array<const double, 5> fdk_options = {
      {fTruncCorFactor, ui.lineEdit_Ramp_HannCut->text().toDouble(),
       ui.lineEdit_Ramp_CosineCut->text().toDouble(),
       ui.lineEdit_Ramp_Hamming->text().toDouble(),
       ui.lineEdit_Ramp_HannCutY->text().toDouble()}};

  std::cout << "Starting RTK fdk" << std::endl;
  FloatImageType::Pointer targetImg = RTKOpenCLFDK(
      spCurImg, m_spCustomGeometry, spacing, sizeOutput, fdk_options);
#endif

  using FOVfilterType =
      rtk::FieldOfViewImageFilter<FloatImageType, FloatImageType>;
  typename FOVfilterType::Pointer FOVfilter = FOVfilterType::New();
  FOVfilter->SetGeometry(m_spCustomGeometry);
  FOVfilter->SetInput(0, targetImg);
  FOVfilter->SetProjectionsStack(spCurImg.GetPointer());
  targetImg = FOVfilter->GetOutput();

  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;

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

  itk::Vector<double, 3U> offset(0.0);

  using TransformType = itk::MatrixOffsetTransformBase<double, 3U, 3U>;
  TransformType::Pointer transform = TransformType::New();
  transform->SetMatrix(CoordChangeMatrix);
  transform->SetOffset(offset);

  using ResampleFilterType =
      itk::ResampleImageFilter<FloatImageType, FloatImageType>;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  // FloatImageType::RegionType fixedImg_Region =
  // fixedImg->GetLargestPossibleRegion().GetSize();

  resampler->SetInput(targetImg);
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(targetImg->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  // resampler->Update();//yktemp Error 2

  // LR flip

  std::cout << "Flip filter is being applied" << std::endl;
  using FilterType = itk::FlipImageFilter<FloatImageType>;
  FilterType::Pointer flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;
  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;
  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  /*FloatImageType::Pointer floatImg = flipFilter->GetOutput();*/
  // const unsigned int Dimension = 3;
  // FinalImageType::Pointer finalImg ;

  using AbsImageFilterType =
      itk::AbsImageFilter<FloatImageType, FloatImageType>;
  AbsImageFilterType::Pointer absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(
      flipFilter->GetOutput()); // 20140206 modified it was a bug

  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  /*Total Variation Filter*/
  castFilter->SetInput(multiplyImageFilter->GetOutput());
  castFilter->Update(); // YK20150109

  UShortImageType::Pointer tmpReconImg;
  // if all 0 0 0 don't do the median filtering

  itk::TimeProbe reconTimeProbe;
  reconTimeProbe.Start();

  std::cout << "Reconstructing the image.. please wait..." << std::endl;
  UShortImageType::SizeType indexRadius{};
  indexRadius[0] = ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
  indexRadius[1] = ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
  indexRadius[2] = ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along y
  if (ui.checkBox_PostMedianOn->isChecked() &&
      (indexRadius[0] != 0 || indexRadius[1] != 0 || indexRadius[2] != 0)) {
    using MedianFilterType =
        itk::MedianImageFilter<UShortImageType,
                               UShortImageType>; // TODO(AGA): CUDA THIS!!

    MedianFilterType::Pointer medFilter = MedianFilterType::New();

    // YKTEMP20141218 S
    // typedef itk::MeanImageFilter<UShortImageType, UShortImageType >
    // FilterType; FilterType::Pointer medFilter = FilterType::New();
    // YKTEMP20141218 E

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

  reconTimeProbe.Stop();
  std::cout << "It took " << reconTimeProbe.GetMean() << ' '
            << reconTimeProbe.GetUnit() << std::endl;
  ui.lineEdit_ReconstructionTime->setText(
      QString("%1").arg(reconTimeProbe.GetMean()));

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
  QString outputFilePath = ui.lineEdit_OutputFilePath->text();

  QFileInfo outFileInfo(outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
    ui.lineEdit_Cur3DFileName->setText("FDK-reconstructed volume");
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFilePath.toLocal8Bit().constData());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(m_spCrntReconImg);

    std::cout << "Writing the image to: "
              << outputFilePath.toLocal8Bit().constData() << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

    ui.lineEdit_Cur3DFileName->setText(outputFilePath);
    std::cout << std::endl;
    std::cout << "Output generation was succeeded" << std::endl;
  }

  m_dspYKReconImage->CreateImage(size_trans[0], size_trans[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(size_trans[2] - 1);
  ui.spinBoxReconImgSliceNo->setValue(
      qRound(size_trans[2] / 2.0)); // DrawReconImage is called automatically

  // SLT_ViewRegistration();

  std::cout << "FINISHED!: FDK CBCT reconstruction" << std::endl;
  // if not found, just skip

  // SLT_DrawGraph();
  // 2) Load Geometry file.
  // 3) Prepare all parameters from GUI components
}
#else
void CbctRecon::OpenCLDoReconstructionFDK(enREGI_IMAGES target) {
  std::cout << "This program were not compiled with the CUDA option, please "
               "select OpenCL or CPU instead!"
            << std::endl;
  return;
}
#endif

void CbctRecon::DoReconstructionFDK(enREGI_IMAGES target) {
  if (m_spProjImg3DFloat == nullptr) {
    std::cout << "processed Projection image is not ready yet" << std::endl;
    return;
  }
  // Resampling first --> to save the recon time. 1024 --> 512

  using DuplicatorType = itk::ImageDuplicator<FloatImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(m_spProjImg3DFloat);
  duplicator->Update();

  FloatImageType::Pointer spCurImg =
      duplicator->GetOutput(); // already down sampled

  // Displaced detector weighting // set pipeline //inplace filter
  using DDFType = rtk::DisplacedDetectorImageFilter<FloatImageType>;
  DDFType::Pointer ddf = DDFType::New();

  if (ui.checkBox_UseDDF->isChecked()) {
    ddf->SetInput(spCurImg);
    // ddf->SetGeometry( geometryReader->GetOutputObject() );
    ddf->SetGeometry(m_spCustomGeometry);
    std::cout << "DDF was set in pipeline" << std::endl;

    if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
      ddf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = ddf->GetOutput();
  }

  using PSSFType = rtk::ParkerShortScanImageFilter<FloatImageType>;
  PSSFType::Pointer pssf = PSSFType::New();

  if (ui.checkBox_UsePSSF->isChecked()) {
    // Short scan image filter
    pssf->SetInput(spCurImg);
    // pssf->SetGeometry( geometryReader->GetOutputObject() );
    pssf->SetGeometry(m_spCustomGeometry);
    // pssf->InPlaceOff(); //YKComments: Do not overwrite input image buffer for
    // output
    std::cout << "short scan image filter success" << std::endl;

    if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
      pssf->Update(); // no mememory increas: InPlace Filter
    }

    spCurImg = pssf->GetOutput();
  }

  // Just Before going to FDK recon,
  // Update Projection data and delete old data.

  // Let's duplicate this
  // Original m_spProjImg3D will be deleted after update

  if (ui.checkBox_UpdateAfterFiltering->isChecked()) {
    using DuplicatorType = itk::ImageDuplicator<FloatImageType>;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(spCurImg);
    duplicator->Update();
    m_spProjImg3DFloat = duplicator->GetOutput();

    SetMaxAndMinValueOfProjectionImage(); // scan m_spProjImg3D and update
                                          // m_fProjImgValueMin, max
    SLT_DrawProjImages();
  }

  // Generate image sources for cone beam CT reconstruction
  using ConstantImageSourceType = rtk::ConstantImageSource<FloatImageType>;

  ConstantImageSourceType::PointType origin;
  ConstantImageSourceType::SpacingType spacing;
  ConstantImageSourceType::SizeType sizeOutput{};

  sizeOutput[0] = ui.lineEdit_outImgDim_AP->text().toInt(); // pixel
  sizeOutput[1] =
      ui.lineEdit_outImgDim_SI->text()
          .toInt(); // Caution!: direction is different in NKI SCAN FIle
  sizeOutput[2] = ui.lineEdit_outImgDim_LR->text().toInt();

  spacing[0] = ui.lineEdit_outImgSp_AP->text().toDouble();
  spacing[1] = ui.lineEdit_outImgSp_SI->text().toDouble();
  spacing[2] = ui.lineEdit_outImgSp_LR->text().toDouble();

  origin[0] = -0.5 * sizeOutput[0] * spacing[0]; // Y in DCM?
  origin[1] = -0.5 * sizeOutput[1] * spacing[1]; // Z in DCM?
  origin[2] = -0.5 * sizeOutput[2] * spacing[2]; // X in DCM?

  ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();
  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);
  constantImageSource->SetSize(sizeOutput);
  constantImageSource->SetConstant(0.0); // initial value

  double fTruncCorFactor =
      ui.lineEdit_Ramp_TruncationCorrection->text().toDouble();
  double fHannCut = ui.lineEdit_Ramp_HannCut->text().toDouble();
  double fCosineCut = ui.lineEdit_Ramp_CosineCut->text().toDouble();
  double fHamming = ui.lineEdit_Ramp_Hamming->text().toDouble();
  double fHannCutY = ui.lineEdit_Ramp_HannCutY->text().toDouble();

  if (fTruncCorFactor > 0.0 && target == REGISTER_COR_CBCT) {
    std::cout << "Warning! Truncation factor is " << fTruncCorFactor
              << ". Regardless of previous setting, this factor should not be "
                 "0 for scatter corrected CBCT. Now zero value is applied."
              << std::endl;
    fTruncCorFactor = 0.0;
  }

  // YKTEMP
  std::cout << "fTruncCorFactor =" << fTruncCorFactor << std::endl;
  // This macro sets options for fdk filter which I can not see how to do better
  // because TFFTPrecision is not the same, e.g. for CPU and CUDA (SR)
#define SET_FELDKAMP_OPTIONS(f)                                                \
  f->SetInput(0, constantImageSource->GetOutput());                            \
  f->SetInput(1, spCurImg);                                                    \
  f->SetGeometry(m_spCustomGeometry);                                          \
  f->GetRampFilter()->SetTruncationCorrection(fTruncCorFactor);                \
  f->GetRampFilter()->SetHannCutFrequency(fHannCut);                           \
  f->GetRampFilter()->SetCosineCutFrequency(fCosineCut);                       \
  f->GetRampFilter()->SetHammingFrequency(fHamming);                           \
  f->GetRampFilter()->SetHannCutFrequencyY(fHannCutY);

  // FDK reconstruction filtering
  itk::ImageToImageFilter<FloatImageType, FloatImageType>::Pointer feldkamp;
  using FDKCPUType = rtk::FDKConeBeamReconstructionFilter<FloatImageType>;

  feldkamp = FDKCPUType::New();
  SET_FELDKAMP_OPTIONS(dynamic_cast<FDKCPUType *>(feldkamp.GetPointer()));

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
  using StreamerType =
      itk::StreamingImageFilter<FloatImageType, FloatImageType>;
  StreamerType::Pointer streamerBP = StreamerType::New();
  streamerBP->SetInput(feldkamp->GetOutput());
  streamerBP->SetNumberOfStreamDivisions(
      1); // YK: 1 in example code from "rtkfdk"

  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;

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
  FloatImageType::Pointer targetImg = streamerBP->GetOutput();

  /* 3) Configure transform */
  using TransformType = itk::Euler3DTransform<double>;
  TransformType::Pointer transform = TransformType::New();

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
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  // FloatImageType::RegionType fixedImg_Region =
  // fixedImg->GetLargestPossibleRegion().GetSize();

  resampler->SetInput(targetImg);
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(targetImg->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  // resampler->Update();//yktemp Error 2

  // LR flip

  std::cout << "LR flip filter is being applied" << std::endl;

  using FilterType = itk::FlipImageFilter<FloatImageType>;

  FilterType::Pointer flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;

  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;

  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  /*FloatImageType::Pointer floatImg = flipFilter->GetOutput();*/
  // const unsigned int Dimension = 3;
  // FinalImageType::Pointer finalImg ;

  using AbsImageFilterType =
      itk::AbsImageFilter<FloatImageType, FloatImageType>;
  AbsImageFilterType::Pointer absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(
      flipFilter->GetOutput()); // 20140206 modified it was a bug
  // absImgFilter->SetInput(resampler->GetOutput());

  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  // typedef unsigned short FinalPixelType;
  // typedef itk::Image< FinalPixelType, 3 > FinalImageType;

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());
  // castFilter->Update(); //YK20150109

  UShortImageType::SizeType indexRadius{};
  indexRadius[0] = ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
  indexRadius[1] = ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
  indexRadius[2] = ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along y

  UShortImageType::Pointer tmpReconImg;
  // if all 0 0 0 don't do the median filtering

  itk::TimeProbe reconTimeProbe;
  reconTimeProbe.Start();

  std::cout << "Reconstructing the image.. please wait..." << std::endl;

  if (ui.checkBox_PostMedianOn->isChecked() &&
      (indexRadius[0] != 0 || indexRadius[1] != 0 || indexRadius[2] != 0)) {
    using MedianImageFilterType =
        itk::MedianImageFilter<UShortImageType, UShortImageType>;
    MedianImageFilterType::Pointer medFilter = MedianImageFilterType::New();

    // YKTEMP20141218 S
    // typedef itk::MeanImageFilter<UShortImageType, UShortImageType >
    // FilterType; FilterType::Pointer medFilter = FilterType::New();
    // YKTEMP20141218 E

    // this is radius. 1 --> median window 3
    std::cout << "Post median(3D) filtering is in the pipeline..Size(radius X "
                 "Y Z) is = "
              << indexRadius << std::endl;

    medFilter->SetRadius(indexRadius);
    medFilter->SetInput(castFilter->GetOutput());
    medFilter->Update(); // Error here!

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

  reconTimeProbe.Stop();
  std::cout << "It took " << reconTimeProbe.GetMean() << ' '
            << reconTimeProbe.GetUnit() << std::endl;
  ui.lineEdit_ReconstructionTime->setText(
      QString("%1").arg(reconTimeProbe.GetMean()));

  // std::cout << "TestYK" << std::endl;
  // std::cout << tmpReconImg->GetRequestedRegion().GetSize() << std::endl;
  // std::cout << "Before InPlace Off: " <<  castFilter->GetInPlace() <<
  // std::endl;	//0  castFilter->InPlaceOff();  std::cout << "After InPlace
  // Off: "
  // <<  castFilter->GetInPlace() << std::endl; //0  Because Input Output format
  // are different, this filter cannot do InPlace function.

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
  QString outputFilePath = ui.lineEdit_OutputFilePath->text();

  QFileInfo outFileInfo(outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
    ui.lineEdit_Cur3DFileName->setText("FDK-reconstructed volume");
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFilePath.toLocal8Bit().constData());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(m_spCrntReconImg);

    std::cout << "Writing the image to: "
              << outputFilePath.toLocal8Bit().constData() << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

    ui.lineEdit_Cur3DFileName->setText(outputFilePath);
    std::cout << std::endl;
    std::cout << "Output generation was succeeded" << std::endl;
  }

  // feldkamp->PrintTiming(std::cout); Deprecated in rtk 1.4
  m_dspYKReconImage->CreateImage(size_trans[0], size_trans[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(size_trans[2] - 1);
  ui.spinBoxReconImgSliceNo->setValue(
      qRound(size_trans[2] / 2.0)); // DrawReconImage is called automatically

  // For 2D image display
  // SLT_DrawReconImage();
  // m_spReconImg = writer->GetOutput();

  /*FloatImageType::SizeType AfterReconSize =
  m_spProjImg3D->GetBufferedRegion().GetSize(); std::cout << "AfterReconSize  "
  << AfterReconSize[0] << "	"
  << AfterReconSize[1] << "	"
  << AfterReconSize[2] << std::endl;*/

  /*ui.radioButton_graph_recon->setChecked(true);
  SLT_InitializeGraphLim();
  SLT_DrawReconImage();	*/

  // SLT_ViewRegistration();

  std::cout << "FINISHED!: FDK CBCT reconstruction" << std::endl;
  // if not found, just skip

  // SLT_DrawGraph();
  // 2) Load Geometry file.
  // 3) Prepare all parameters from GUI components
}
