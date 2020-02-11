// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#if USE_TINYREFL
#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#include <tinyrefl/api.hpp>
#else
#include "cbctrecon.h"
#endif

#define USE_AVX false
#if USE_AVX
#include <immintrin.h>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#endif

#include <charconv>
#include <cstdio>
#include <filesystem>
#include <valarray>

// ITK
#include <itkAbsImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryFillholeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkDOMNodeXMLReader.h>
#include <itkImageDuplicator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkMaskImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiplyImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkResampleImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkTimeProbe.h>

// RTK includes
#include <rtkConstantImageSource.h>
#include <rtkFieldOfViewImageFilter.h>
#include <rtkProjectionsReader.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include <rtkThreeDCircularProjectionGeometryXMLFileWriter.h>
#include <rtkVarianObiGeometryReader.h>
#include <rtkVarianProBeamGeometryReader.h>

#if USE_CUDA
#include <itkCudaImage.h>
using CUDAFloatImageType = itk::CudaImage<float, 3U>;
#endif // USE_CUDA

// Local
#include "AG17RGBAImage.h"
#include "OpenCL/ImageFilters.hpp"
#include "StructureSet.h"
#include "WEPL.h"
#include "YK16GrayImage.h"
#include "cbctrecon_compute.h"
#include "cbctrecon_io.h"
#include "free_functions.h"

namespace fs = std::filesystem;
using namespace std::literals;

CbctRecon::CbctRecon() {

  m_dspYKReconImage = std::make_unique<YK16GrayImage>();
  m_dspYKImgProj = std::make_unique<YK16GrayImage>();

  // Badpixmap;
  m_pImgOffset = std::make_unique<YK16GrayImage>(DEFAULT_ELEKTA_PROJ_WIDTH,
                                                 DEFAULT_ELEKTA_PROJ_HEIGHT);
  m_pImgGain = std::make_unique<YK16GrayImage>(DEFAULT_ELEKTA_PROJ_WIDTH,
                                               DEFAULT_ELEKTA_PROJ_HEIGHT);
  // Prepare Raw image

  m_iTmpIdx = 60;

  m_fProjImgValueMax = 0.0; // value of float image
  m_fProjImgValueMin = 0.0;

  m_arrYKBufProj.clear();
  m_iCntSelectedProj = 0;

  m_structures = std::make_unique<StructureSet>();

  m_iFixedOffset_ScatterMap = 10000; // fixed! allows negative value of scatter
  // m_iFixedOffset_ScatterMap = 0;//fixed! allows negative value of scatter
  m_fResampleF = 1.0;
  // m_fProjSpacingX = 0.4; // DEFAULT, will be updated during Load Proj
  // selected m_fProjSpacingY = 0.4;

  m_strPathDirDefault = fs::current_path();
  std::cout << "Current Default Dir: " << m_strPathDirDefault << std::endl;

  // shell test
  // std::string strCurFolder = "H:\\lib\\rtk\\NightlyBUILD64\\bin\\Release";
  // std::string strCommand = std::string("explorer %1").arg(strCurFolder);
  // std::cout << strCommand.toLocal8Bit().constData() << std::endl;
  //::system(strCommand.toLocal8Bit().constData());
  //::system("rtkfdk");

  //	if
  //(QProcess::execute(std::string("H:\\lib\\rtk\\NightlyBUILD64\\bin\\Release\\rtkfdk"))
  //< 0) 	qDebug() << "Failed to run";
  m_bMacroContinue = true;
}

void CbctRecon::ReleaseMemory() {
  if (!m_arrYKImage.empty()) {
    // m_iImgCnt = 0;
    m_arrYKImage.clear();
  }

  if (!m_arrYKBufProj.empty()) {
    m_arrYKBufProj.clear();
    m_iCntSelectedProj = 0;
  }
}

// Hexa name ->decimal name

void RenameFromHexToDecimal(const std::vector<fs::path> &filenameList) {
  const auto size = filenameList.size();

  for (auto i = 0; i < size; i++) {
    const auto &crntFilePath = filenameList.at(i);
    auto dir = fs::absolute(crntFilePath);
    auto fileBase = crntFilePath.stem();
    auto newBaseName = crl::HexStr2IntStr(fileBase.string());
    auto extStr = crntFilePath.extension();

    auto newFileName = newBaseName.append(".").append(extStr.string());
    auto newPath = fs::absolute(dir) / newFileName;

    // extract former part
    fs::rename(crntFilePath, newPath);
  }
  // Extract
}

bool CbctRecon::FillProjForDisplay(const int slice_number) {
  // Using slice iterator,
  // 1) Find the slice requested
  // 2) get dimension to create 2DYK16Image
  // 3) copy slice region to YK16 Image --> Cating: float to USHORT

  if (m_spProjImg3DFloat == nullptr) {
    return false;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it(
      m_spProjImg3DFloat, m_spProjImg3DFloat->GetBufferedRegion());

  auto imgSize =
      m_spProjImg3DFloat->GetBufferedRegion().GetSize(); // 1016x1016 x z

  const auto width = imgSize[0];
  const auto height = imgSize[1];
  m_dspYKImgProj->CreateImage(width, height, 0);

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  const auto realValGap = m_fProjImgValueMax - m_fProjImgValueMin;
  m_multiplyFactor = 0.0;

  if (realValGap > 0.0) {
    m_multiplyFactor = 65535.0 / realValGap;
  }

  size_t i_num_slice = 0;
  while (!it.IsAtEnd()) {
    if (i_num_slice == static_cast<size_t>(slice_number)) {
      size_t i_num_height = 0;
      while (!it.IsAtEndOfSlice()) {
        size_t i_num_width = 0;
        while (!it.IsAtEndOfLine()) {
          // double tmpVal = it.Get()*multiplyFactor;
          const double tmp_val = it.Get();

          m_dspYKImgProj->m_pData[i_num_width + width * i_num_height] =
              static_cast<unsigned short>((tmp_val - m_fProjImgValueMin) *
                                          m_multiplyFactor);
          // it.Set() doesn't exist in the Const Iterator
          ++it;
          ++i_num_width;
        }
        it.NextLine();
        ++i_num_height;
      }
      // break;
    }
    it.NextSlice();
    ++i_num_slice;
  }

  return true;
}

void CbctRecon::LoadCalibData(std::string &filepath,
                              const enCalibType calib_type) {
  switch (calib_type) {
  case enCalibType::GAIN_CALIB:
    m_pImgGain->LoadRawImage(filepath.c_str(), DEFAULT_ELEKTA_PROJ_WIDTH,
                             DEFAULT_ELEKTA_PROJ_HEIGHT);
    break;
  case enCalibType::OFFSET_CALIB:
    m_pImgOffset->LoadRawImage(filepath.c_str(), DEFAULT_ELEKTA_PROJ_WIDTH,
                               DEFAULT_ELEKTA_PROJ_HEIGHT);
    break;
  case enCalibType::BADPIXEL_CALIB:
    crl::LoadBadPixelMap(m_vPixelReplMap, filepath);
    break;
  default:
    break;
  }
}

void CbctRecon::SetProjDir(std::string &strProjPath) {
  m_strPathGeomXML.clear();
  m_strPathDirDefault = strProjPath;

  const UShortImageType::Pointer spNull;

  m_spCrntReconImg = spNull;     // fixed image // ID: RawCBCT
  m_spRawReconImg = spNull;      // just added --> when file is loaded
  m_spScatCorrReconImg = spNull; // just added --> after scatter correction

  FindAllRelevantPaths(strProjPath);
}

bool CbctRecon::LoadGeometry(fs::path &geomFileInfo,
                             std::vector<std::string> &names) {

  if (!fs::exists(geomFileInfo)) {
    std::cout << "Critical Error! geometry file is not existing. Please retry."
              << std::endl;
    return false;
  }

  const fs::path tmp_rtk_geom_file{"RTKgeometry.xml"};

  if (geomFileInfo.filename() == "_Frames.xml") // this is XVI XML.
  {
    std::cout << "XVI Geometry File was found. This will be temporarily used:"
              << geomFileInfo.filename() << std::endl;
    const auto success =
        LoadXVIGeometryFile(fs::absolute(geomFileInfo)
                                ); // will generate m_spFullGeometry
    if (!success) {
      return false;
    }
  } else if (geomFileInfo.filename() ==
             "ProjectionInfo.xml") // this is OBI XML.
  {
    std::cout
        << "Varian XML Geometry File was found. This will be temporarily used:"
        << geomFileInfo.filename() << std::endl;
    auto reader = rtk::VarianObiGeometryReader::New();
    reader->SetXMLFileName(geomFileInfo.filename().string());
    reader->SetProjectionsFileNames(names);
    reader->UpdateOutputData();
    // Write
    auto xmlWriter = rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
    xmlWriter->SetFilename(tmp_rtk_geom_file.string());
    xmlWriter->SetObject(reader->GetGeometry());
    TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
    std::cout << "RTK standard Geometry XML File was created:"
              << tmp_rtk_geom_file << std::endl;
    m_spFullGeometry = crl::LoadRTKGeometryFile(tmp_rtk_geom_file);
    // ::::::::::::::::::::::::::::LoadXMLGeometryFile(geomPath.toLocal8Bit().constData());
    // //will generate m_spFullGeometry
  } else if (geomFileInfo.filename() == "Scan.xml") // this is XIM XML.
  {
    std::cout << "Varian Xim XML Geometry File was found. This will be "
                 "temporarily used:"
              << geomFileInfo.filename() << std::endl;
    auto reader = rtk::VarianProBeamGeometryReader::New();
    reader->SetXMLFileName(
        fs::absolute(geomFileInfo).string());
    reader->SetProjectionsFileNames(names);
    reader->UpdateOutputData();
    // Write
    auto xmlWriter = rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
    xmlWriter->SetFilename(tmp_rtk_geom_file.string());
    xmlWriter->SetObject(reader->GetGeometry());
    xmlWriter->WriteFile();
    std::cout << "RTK standard Geometry XML File was created:"
              << tmp_rtk_geom_file << std::endl;
    m_spFullGeometry = crl::LoadRTKGeometryFile(tmp_rtk_geom_file);
    std::cout << "Done!";
    // ::::::::::::::::::::::::::::LoadXMLGeometryFile(geomPath.toLocal8Bit().constData());
    // //will generate m_spFullGeometry
  } else {
    std::cout << "RTK standard Geometry XML File was found:"
              << fs::absolute(geomFileInfo)
              << std::endl;
    m_spFullGeometry = crl::LoadRTKGeometryFile(
        fs::absolute(geomFileInfo)
            ); // will generate m_spFullGeometry
  }
  return true;
}

void CbctRecon::LoadSelectedProj(const std::vector<size_t> &exclude_ids,
                                 const std::vector<std::string> &names) {
  // 3) Seletively load projection file

  // Regenerate geometry object
  m_spCustomGeometry = GeometryType::New();

  // for (itIdx =vSelectedIdx.begin() ; itIdx != vSelectedIdx.end() ; itIdx++ )
  // #pragma omp parallel for private(itIdx) shared(m_spCustomGeometry)
  // schedule(static)
  for (auto &it_idx : exclude_ids) {
    // 9 parameters are required
    const auto curSID =
        m_spFullGeometry->GetSourceToIsocenterDistances().at(it_idx);
    const auto curSDD =
        m_spFullGeometry->GetSourceToDetectorDistances().at(it_idx);
    auto curGantryAngle = m_spFullGeometry->GetGantryAngles().at(it_idx);
    const auto kVAng =
        curGantryAngle * 180.0 *
        itk::Math::one_over_pi; // 360 / 2 = 180 radians to degrees
    auto MVAng =
        kVAng - (m_projFormat == enProjFormat::HIS_FORMAT ? 0.0 : 90.0);
    if (MVAng < 0.0) {
      MVAng = MVAng + 360.0;
    }
    curGantryAngle = MVAng;

    const auto curProjOffsetX =
        m_spFullGeometry->GetProjectionOffsetsX().at(it_idx);
    const auto curProjOffsetY =
        m_spFullGeometry->GetProjectionOffsetsY().at(it_idx);

    const auto curOutOfPlaneAngles =
        m_spFullGeometry->GetOutOfPlaneAngles().at(it_idx);
    const auto curInPlaneAngles =
        m_spFullGeometry->GetInPlaneAngles().at(it_idx);

    const auto curSrcOffsetX = m_spFullGeometry->GetSourceOffsetsX().at(it_idx);
    const auto curSrcOffsetY = m_spFullGeometry->GetSourceOffsetsY().at(it_idx);

    m_spCustomGeometry->AddProjection(
        curSID, curSDD, curGantryAngle, curProjOffsetX,
        curProjOffsetY, // Flexmap, For Xim, these are 0
        curOutOfPlaneAngles,
        curInPlaneAngles,              // In elekta and varian, these are 0
        curSrcOffsetX, curSrcOffsetY); // In elekta and varian, these are 0
  }

  std::cout << "Excluded proj count: " << m_vExcludeProjIdx.size() << std::endl;
  std::cout << "Final proj count: " << exclude_ids.size() << std::endl;

  // Regenerate fileNames and geometry object based on the selected indices.

  if (!m_vSelectedFileNames.empty()) {
    m_vSelectedFileNames.clear();
  }

  std::ofstream fout;
  fout.open("DebugFileNames.txt");

  for (auto &it_idx : exclude_ids) {
    const auto &cur_str = names.at(it_idx);
    m_vSelectedFileNames.push_back(cur_str);
    fout << cur_str.c_str() << std::endl;
  }

  fout.close();

  m_iCntSelectedProj =
      m_vSelectedFileNames.size(); // Used to check in ScatterCorr_PrioriCT
}

void CbctRecon::saveHisHeader() {
  if (m_projFormat == enProjFormat::HIS_FORMAT) {
    std::cout << "Copying the HIS info to buffer." << std::endl;
    m_arrYKBufProj.resize(m_vSelectedFileNames.size());
    auto it_selected = m_vSelectedFileNames.begin();
    for (auto it = m_arrYKBufProj.begin();
         it != m_arrYKBufProj.end() &&
         it_selected != m_vSelectedFileNames.end();
         ++it, ++it_selected) {
      it->m_strFilePath = it_selected->c_str();
      it->CopyHisHeader(it_selected->c_str());
    }
  }
}

void CbctRecon::NormalizeProjections(
    const FloatImageType::Pointer &reader_output) {

  auto originalMax = -1.0;
  auto originalMin = -1.0;
  //                   µ/rho (cm^2/g)  * rho (g/cm^3) * path (cm)
  const auto theoreticalMin =
      0.1541 * 1.225e-3 *
      sqrt(pow(m_spCustomGeometry->GetSourceToDetectorDistances()[0], 2) +
           pow(m_spCustomGeometry->GetSourceToDetectorDistances()[1], 2) +
           pow(m_spCustomGeometry->GetSourceToDetectorDistances()[2], 2)) *
      0.1; // mm -> cm

  const auto correctionValue = crl::GetMaxAndMinValueOfProjectionImage(
      originalMax, originalMin, reader_output); // , theoreticalMin);
  std::cout << "Reader Max, Min=" << originalMax << "	" << originalMin
            << std::endl;

  if (correctionValue > 1000.0) {
    auto add_filter = itk::AddImageFilter<FloatImageType, FloatImageType,
                                          FloatImageType>::New();
    add_filter->SetInput(reader_output);
    add_filter->SetConstant(-correctionValue);
    if (originalMax - originalMin > log(65535.0f) - theoreticalMin) {
      auto mul_filter = itk::MultiplyImageFilter<FloatImageType, FloatImageType,
                                                 FloatImageType>::New();
      mul_filter->SetInput(add_filter->GetOutput());
      mul_filter->SetConstant((log(65535.0f) - theoreticalMin) /
                              (originalMax - originalMin));
      mul_filter->Update();
      m_spProjImg3DFloat =
          mul_filter->GetOutput(); // 1024 1024, line integ image
    } else {
      add_filter->Update();
      m_spProjImg3DFloat =
          add_filter->GetOutput(); // 1024 1024, line integ image
    }

    /* OpenCL is slower than ITK for the realistic image sizes
    if (correctionValue > 1000.0) {
      if (originalMax - originalMin > log(65535.0f) - theoreticalMin) {
        OpenCL_AddConst_MulConst_InPlace(
            static_cast<cl_float *>(reader->GetOutput()->GetBufferPointer()),
            reader->GetOutput()->GetLargestPossibleRegion().GetSize(),
            static_cast<cl_float>(-correctionValue),
            static_cast<cl_float>((log(65535.0f) - theoreticalMin) /
                                  (originalMax - originalMin)));
      } else {
        OpenCL_AddConst_InPlace(
            static_cast<cl_float *>(reader->GetOutput()->GetBufferPointer()),
            reader->GetOutput()->GetLargestPossibleRegion().GetSize(),
            static_cast<cl_float>(-correctionValue));
      }
      */
    // Reset min max:
    originalMax = -1.0;
    originalMin = -1.0;
    if (crl::GetMaxAndMinValueOfProjectionImage(originalMax, originalMin,
                                                m_spProjImg3DFloat) > -1000.0) {
      std::cout << "Reader Max, Min=" << originalMax << "	" << originalMin
                << std::endl;
    }
  }
}

// True if projections were resampled
bool CbctRecon::ResampleProjections(double &resample_factor) {

  std::cout << "ProjectionReader Get Spacing : "
            << m_spProjImg3DFloat->GetSpacing() << std::endl;

  // m_fProjSpacingX = m_spProjImg3DFloat->GetSpacing()[0];
  // m_fProjSpacingY = m_spProjImg3DFloat->GetSpacing()[1];

  if (resample_factor > 1 || resample_factor <= 0) {
    std::cout << "wrong resample factor. reset to 1.0" << std::endl;
    resample_factor = 1.0;
    return false;
  }

  if (fabs(resample_factor - 1.0) > 0.001) {
    ResampleItkImage(
        m_spProjImg3DFloat, m_spProjImg3DFloat,
        resample_factor); // was! BROKEN AF for .his where input size
                          // != 1024 (tested with 1016) -> outputs
                          // offset -inputoffset/refactor^2 and 4
                          // pixels too few in x and y
  }
  return true;
}

void CbctRecon::GetExcludeIndexByNames(
    const fs::path &outlierListPath,
    std::vector<std::string> &vProjFileFullPath,
    std::vector<int> &vExcludeIdx) const {
  std::ifstream fin;
  fin.open(outlierListPath, std::ios::in);
  if (static_cast<int>(fin.fail()) == 1) {
    return;
  }

  char str[MAX_LINE_LENGTH];

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    if (strlen(&str[0]) < 1) {
      continue;
    }

    std::cout << "Outlier file: " << &str[0] << std::endl;

    auto curIdx = 0;
    for (auto &it : vProjFileFullPath) {
      if (strstr(it.c_str(), static_cast<const char *>(str)) != nullptr) {
        std::cout << "Detected in the list. Index = " << curIdx << std::endl;
        vExcludeIdx.push_back(curIdx);
        break;
      }
      curIdx++;
    }
  }
  fin.close();
}

void CbctRecon::SetMaxAndMinValueOfProjectionImage() // should be called
                                                     // whenever proj image is
                                                     // changed
{
  if (m_iImgCnt > 0) {
    m_fProjImgValueMax = 65535;
    m_fProjImgValueMin = 0;
    return;
  }

  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it(
      m_spProjImg3DFloat, m_spProjImg3DFloat->GetBufferedRegion());

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  m_fProjImgValueMin = 65535.0;
  m_fProjImgValueMax = -9999.0;

  while (!it.IsAtEnd()) {
    while (!it.IsAtEndOfSlice()) {
      while (!it.IsAtEndOfLine()) {
        const double tmpVal = it.Get();

        if (m_fProjImgValueMax < tmpVal) {
          m_fProjImgValueMax = tmpVal;
        }

        if (m_fProjImgValueMin > tmpVal) {
          m_fProjImgValueMin = tmpVal;
        }

        ++it;
      }
      it.NextLine();
    }
    it.NextSlice();
  }
}

void CbctRecon::PostApplyFOVDispParam(const float physPosX,
                                      const float physPosY,
                                      const float physRadius,
                                      const float physTablePosY) const {
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  auto origin = m_spCrntReconImg->GetOrigin();
  auto spacing = m_spCrntReconImg->GetSpacing();
  // UShortImageType::SizeType size =
  // m_spCrntReconImg->GetBufferedRegion().GetSize();

  const auto pixPosX = qRound((physPosX - static_cast<double>(origin[0])) /
                              static_cast<double>(spacing[0]));
  const auto pixPosY = qRound((physPosY - static_cast<double>(origin[1])) /
                              static_cast<double>(spacing[1]));

  const auto pixRadius = qRound(physRadius / static_cast<double>(spacing[0]));

  // int pixWidth = qRound((qreal) size[0]);
  // int pixHeight = qRound((qreal) size[1]);

  const auto pixTableY =
      qRound((physTablePosY - static_cast<double>(origin[1])) /
             static_cast<double>(spacing[1]));

  if (pixPosX >= 0 && pixPosY < m_dspYKReconImage->m_iWidth && pixPosY >= 0 &&
      pixPosY < m_dspYKReconImage->m_iHeight && pixRadius > 0 &&
      pixRadius < m_dspYKReconImage->m_iWidth && pixTableY >= 0 &&
      pixTableY < m_dspYKReconImage->m_iHeight) {
    m_dspYKReconImage->m_ptFOVCenter.setX(pixPosX); // data pos
    m_dspYKReconImage->m_ptFOVCenter.setY(pixPosY);

    m_dspYKReconImage->m_iFOVRadius = pixRadius;
    m_dspYKReconImage->m_iTableTopPos = pixTableY;
  }
}

void CbctRecon::CropSupInf(UShortImageType::Pointer &sp_Img,
                           float physPosInfCut, float physPosSupCut) {
  if (sp_Img == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  auto origin = sp_Img->GetOrigin();
  auto spacing = sp_Img->GetSpacing();
  auto size = sp_Img->GetBufferedRegion().GetSize();

  std::cout << "Old Origin" << origin << std::endl;
  std::cout << "Old spacing" << spacing << std::endl;
  std::cout << "Old size" << size << std::endl;

  UShortImageType::SizeType sizeLower{}, sizeUpper{}; // not index. this is
                                                      // width of pixels that
                                                      // will be taken away.
  sizeLower[0] = 0;
  sizeLower[1] = 0;
  sizeLower[2] = 0;
  /*indexUpper[0] = size[0] - 1;
  indexUpper[1] = size[1] - 1;
  indexUpper[2] = size[2] - 1;*/
  sizeUpper[0] = 0;
  sizeUpper[1] = 0;
  sizeUpper[2] = 0;

  const auto minPosSI = origin[2];
  const auto maxPosSI = origin[2] + (size[2] - 1) * spacing[2];

  if (minPosSI >= physPosInfCut) {
    physPosInfCut = minPosSI;
  }
  if (maxPosSI <= physPosSupCut) {
    physPosSupCut = maxPosSI;
  }

  if (physPosSupCut <= physPosInfCut) {
    return;
  }

  ////calc index
  sizeLower[2] = qRound((physPosInfCut - minPosSI) / spacing[2]);
  sizeUpper[2] = qRound((maxPosSI - physPosSupCut) / spacing[2]);
  //
  using CropImageFilterType =
      itk::CropImageFilter<UShortImageType, UShortImageType>;
  auto CropFilter = CropImageFilterType::New();

  CropFilter->SetInput(sp_Img);
  CropFilter->SetLowerBoundaryCropSize(sizeLower);
  CropFilter->SetUpperBoundaryCropSize(sizeUpper);

  CropFilter->Update();

  if (sp_Img == m_spRawReconImg) {
    sp_Img = CropFilter->GetOutput();
    m_spRawReconImg = sp_Img;
  }

  if (sp_Img == m_spRefCTImg) {
    sp_Img = CropFilter->GetOutput();
    m_spRefCTImg = sp_Img;
    m_spManualRigidCT = sp_Img;
  }

  const auto origin_new = sp_Img->GetOrigin();
  const auto spacing_new = sp_Img->GetSpacing();
  const auto size_new = sp_Img->GetBufferedRegion().GetSize();

  // origin_new[2] = physPosInfCut;
  // sp_Img->SetOrigin(origin_new);

  std::cout << "New Origin" << origin_new << std::endl;
  std::cout << "New spacing" << spacing_new << std::endl;
  std::cout << "New size" << size_new << std::endl;

  std::cout << "LowPos[mm, index] = " << physPosInfCut << ", " << sizeLower[2]
            << std::endl;
  std::cout << "UpperPos[mm, index] = " << physPosSupCut << ", " << sizeUpper[2]
            << std::endl;
  std::cout << "Cropping SI has been successfully done." << std::endl;

  // Result: same image after cropping
  /*
      sizeDiff[0] = CropFilter->GetOutput()->GetBufferedRegion().GetSize()[0] -
     pYKImageROI->m_iWidth; sizeDiff[1] =
     CropFilter->GetOutput()->GetBufferedRegion().GetSize()[1] -
     pYKImageROI->m_iHeight;

      if (sizeDiff[0] != 0 || sizeDiff[1] != 0)
      {
      std::cout << "Cross-correlation error! template size is not matching ROI
     image even after cropping" << std::endl; return;
      }

      rescaleFilter->SetInput(CropFilter->GetOutput());*/
}

// mm

void CbctRecon::DoBeamHardeningCorrection() const {
  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  // FloatImageType m_spProjImg3D: float image

  // double crntVal = 0.0;
  // double corrF = 0.0;

  // REMEMBER to change in the model functions in crl, here is only for
  // debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // HND_FORMAT:
  auto poly3_a = 6.0e-08;
  auto poly3_b = -1.0e-08;
  auto poly3_c = -5.0e-07;
  auto poly3_d = 8.0e-01;
  auto poly3_e = 1.47;

  switch (m_projFormat) {
  case enProjFormat::HIS_FORMAT:
    poly3_a = 9.321e-05;
    poly3_b = -2.609e-03;
    poly3_c = 3.374e-02;
    poly3_d = 9.691e-01;
    poly3_e = 0.0;
    break;
  case enProjFormat::HND_FORMAT:
    break; // used to initialize
  case enProjFormat::XIM_FORMAT:
    poly3_a = 6.0e-8;
    poly3_b = 9.0e-5;
    poly3_c = 1.0e-2;
    poly3_d = 0.8;
    poly3_e = -1.47;
    break;
  }
  // If Hnd:

  // double corrVal = 0.0;

  // Shortening factor 0.9 is applied
  /*double poly3_a = 11.24e-05;
  double poly3_b = -29.28e-04;
  double poly3_c = 35.48e-03;
  double poly3_d = 9.701e-01;*/

  // Shortening factor 0.7 is applied
  /*double poly3_a = 1.660e-04;
  double poly3_b = -3.699e-03;
  double poly3_c = 3.923e-02;
  double poly3_d = 9.727e-01;*/

  std::cout << "Beam hardening corrF poly curve:" << poly3_a << "	"
            << poly3_b << "	" << poly3_c << "	" << poly3_d << "  "
            << poly3_e << std::endl;

  const auto pImgBuffer = m_spProjImg3DFloat->GetBufferPointer();
  const auto pImgSize =
      m_spProjImg3DFloat->GetLargestPossibleRegion().GetSize();
  const unsigned int nPix = pImgSize[0] * pImgSize[1] * pImgSize[2];

  switch (m_projFormat) {
  case enProjFormat::HIS_FORMAT:
    crl::BeamHardening<enProjFormat::HIS_FORMAT>(pImgBuffer, nPix);
    break;
  case enProjFormat::HND_FORMAT:
    crl::BeamHardening<enProjFormat::HND_FORMAT>(pImgBuffer, nPix);
    break;
  case enProjFormat::XIM_FORMAT:
    crl::BeamHardening<enProjFormat::XIM_FORMAT>(pImgBuffer, nPix);
    break;
  }
}

void CbctRecon::BowtieByFit(const bool fullfan,
                            const std::vector<std::string> &params) const {
  if (params.size() != 4 && !fullfan) {
    std::cout << "Wrong number of arguments!" << std::endl
              << "Must be a;b;c;d -> d. / (1 + exp(-b.*(x - a))) + c"
              << std::endl;
    return;
  }
  if (params.size() != 5 && fullfan) {
    std::cout << "Wrong number of arguments!" << std::endl
              << "Must be a;b;c;d;e -> c - sqrt(abs(a^2-((x+/-e)*d-b)^2)) * "
                 "heaviside((x+/-e)*d-b+a) * heaviside(-((x+/-e)*d-b-a))"
              << std::endl;
    return;
  }
  const auto poly3_a = crl::from_string<double>(params.at(0))
                           .value_or(0.0); // 264.6; //comboBox_fBTcor
  const auto poly3_b =
      crl::from_string<double>(params.at(1)).value_or(0.0); // 0.06258;
  const auto poly3_c =
      crl::from_string<double>(params.at(2)).value_or(0.0); // 2.502;
  const auto poly3_d =
      crl::from_string<double>(params.at(3)).value_or(0.0); // 1.455;
  auto poly3_e = 0.0;
  if (fullfan) {
    poly3_e = crl::from_string<double>(params.at(4)).value_or(0.0);
  }

  auto imgSize = m_spProjImg3DFloat->GetLargestPossibleRegion().GetSize();
  if (fullfan) {
    std::cout << "Bow-tie correction curve:" << poly3_d << " / ( 1 + exp(-"
              << poly3_b << " * (x - " << poly3_a << "))) + " << poly3_c
              << " + " << poly3_d << " / ( 1 + exp(" << poly3_b << " * (x - "
              << poly3_e << ")))" << std::endl;
  } else {
    std::cout << "Bow-tie correction curve:" << poly3_d << " / ( 1 + exp(-"
              << poly3_b << " * (x - " << poly3_a << "))) + " << poly3_c
              << std::endl;
  }

  using iteratorType = itk::ImageRegionIteratorWithIndex<FloatImageType>;
  iteratorType it(m_spProjImg3DFloat, m_spProjImg3DFloat->GetBufferedRegion());

  it.GoToBegin();
  if (fullfan) {
    while (!it.IsAtEnd()) {
      const auto crntVal =
          static_cast<double>(it.Get()); // (65535 / exp(it.Get())); //raw
                                         // mu_t = ln(65535/I) <-> I =
                                         // 65535 / exp(mu_t)
      const auto x_idx =
          static_cast<double>(it.GetIndex()[0]) *
          (512.0 / imgSize[0]); // 512 from current fit -> conversion
                                // to be consistent with downResFactor
      // if (crntVal > (poly3_c - poly3_a)), negative values are fine, don't
      // worry
      const auto corrF = crl::fullFan_Function(poly3_a, poly3_b, poly3_c,
                                               poly3_d, poly3_e, x_idx);
      it.Set(static_cast<float>(crntVal -
                                corrF)); // (log(65535 / (crntVal - corrF))));
      ++it;
    }
  } else {
    while (!it.IsAtEnd()) {
      const auto crntVal =
          static_cast<double>(it.Get()); // (65535 / exp(it.Get()));
                                         // //raw mu_t = ln(65535/I)
                                         // <-> I = 65535 / exp(mu_t)
      const auto x_idx =
          static_cast<double>(it.GetIndex()[0]) *
          (512.0 / imgSize[0]); // 512 from current fit -> conversion
                                // to be consistent with downResFactor
      // if (crntVal > poly3_c), negative values are fine, don't worry
      const auto corrF =
          poly3_d / (1.0 + exp(-poly3_b * (x_idx - poly3_a))) + poly3_c;
      it.Set(static_cast<float>(crntVal -
                                corrF)); // (log(65535 / (crntVal - corrF))));
      ++it;
    }
  }
}

void CbctRecon::Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                                   UShortImageType::Pointer &spMovingImg,
                                   const enPLANE enPlane, const double pos,
                                   YK16GrayImage &YKFixed,
                                   YK16GrayImage &YKMoving) const {
  if (spFixedImg == nullptr || spMovingImg == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      spFixedImg, spFixedImg->GetBufferedRegion());

  auto imgSize = spFixedImg->GetBufferedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::SizeType imgSizeBuf =
  // spFixedImg->GetBufferedRegion().GetSize(); //1016x1016 x z
  // UShortImageType::SizeType imgSizeLargest =
  // spFixedImg->GetLargestPossibleRegion().GetSize(); //1016x1016 x z

  auto imgOrigin = spFixedImg->GetOrigin();
  auto imgSpacing = spFixedImg->GetSpacing();

  auto width = imgSize[0];
  auto height = imgSize[1];
  auto i_req_slice =
      static_cast<size_t>(qRound((pos - imgOrigin[2]) / imgSpacing[2]));
  auto i_cnt_slice = imgSize[2];

  // For moving image
  using ResampleFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  auto filter = ResampleFilterType::New();

  filter->SetInput(spMovingImg);

  using TransformType = itk::AffineTransform<double, 3>;
  const auto transform = TransformType::New();
  filter->SetTransform(transform);

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, double>;
  const auto interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);
  filter->SetDefaultPixelValue(0);

  UShortImageType::DirectionType direction;
  direction.SetIdentity();
  filter->SetOutputDirection(direction);

  auto movingSpacing = imgSpacing;
  auto movingOrigin = imgOrigin;
  auto movingSize = imgSize;

  switch (enPlane) {
  case enPLANE::PLANE_AXIAL:
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?

    movingSpacing[2] = 1.0;
    movingOrigin[2] = pos;
    movingSize[2] = 1;
    // Resample Here! make corresponding 2D image for Moving image
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);

    break;
  case enPLANE::PLANE_FRONTAL:
    width = imgSize[0];
    height = imgSize[2];
    i_cnt_slice = imgSize[1];
    i_req_slice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[1] = 1.0;
    movingOrigin[1] = pos;
    movingSize[1] = 1;

    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[2]);
    break;
  case enPLANE::PLANE_SAGITTAL:
    width = imgSize[1];
    height = imgSize[2];
    i_cnt_slice = imgSize[0];
    i_req_slice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
    it.SetFirstDirection(1);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[0] = 1.0;
    movingOrigin[0] = pos;
    movingSize[0] = 1;

    YKFixed.SetSpacing(imgSpacing[1], imgSpacing[2]);
    break;

  default:
    std::cout << "default should not passed by" << std::endl;
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  }

  filter->SetOutputSpacing(movingSpacing);
  filter->SetOutputOrigin(movingOrigin);
  filter->SetSize(movingSize);
  filter->Update();

  YKFixed.CreateImage(width, height, 0);
  // std::cout << "Before MovingImg Creation " << std::endl;

  YKMoving.CreateImage(width, height, 0); // exactly same dimension

  itk::ImageRegionConstIterator<UShortImageType> itMoving(
      filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
  size_t cnt = 0;

  // this simple code will cause flip of the image in frontal and sagittal image
  for (itMoving.GoToBegin(); !itMoving.IsAtEnd(); ++itMoving) {
    YKMoving.m_pData[cnt] = itMoving.Get();
    ++cnt;
  }
  if (enPlane != enPLANE::PLANE_AXIAL) {
    YKMoving.EditImage_Flip();
  }

  // std::cout << "tot pixel no: " << cnt << std::endl;

  // YK16GrayImage::CopyItkImage2YKImage(filter->GetOutput(), &YKMoving);

  // std::cout << "After MovingImg Creation " << std::endl;

  it.GoToBegin();

  size_t iNumSlice = 0;

  if (i_req_slice >= i_cnt_slice) {
    return;
  }

  while (!it.IsAtEnd()) {
    if (iNumSlice == i_req_slice) {
      size_t iNumHeight = 0;

      while (!it.IsAtEndOfSlice()) {
        size_t iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          const auto fixedImgVal = it.Get();

          if (enPlane == enPLANE::PLANE_AXIAL) {
            YKFixed.m_pData[iNumWidth + width * iNumHeight] = fixedImgVal;
          } else {
            YKFixed.m_pData[iNumWidth + width * (height - iNumHeight - 1)] =
                fixedImgVal;
          }

          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      break;
    }
    it.NextSlice();
    iNumSlice++;
  }

  //    std::cout << "YK Images were filled" << std::endl;
}

// Actually just an overload of the above function
void CbctRecon::Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                                   UShortImageType::Pointer &spMovingImg,
                                   const enPLANE enPlane, const double pos,
                                   AG17RGBAImage &YKFixed,
                                   AG17RGBAImage &YKMoving) const {
  if (spFixedImg == nullptr || spMovingImg == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      spFixedImg, spFixedImg->GetBufferedRegion());

  auto imgSize = spFixedImg->GetBufferedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::SizeType imgSizeBuf =
  // spFixedImg->GetBufferedRegion().GetSize(); //1016x1016 x z
  // UShortImageType::SizeType imgSizeLargest =
  // spFixedImg->GetLargestPossibleRegion().GetSize(); //1016x1016 x z

  auto imgOrigin = spFixedImg->GetOrigin();
  auto imgSpacing = spFixedImg->GetSpacing();

  auto width = imgSize[0];
  auto height = imgSize[1];
  auto iCntSlice = imgSize[2];
  auto iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);

  // For moving image
  using ResampleFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  auto filter = ResampleFilterType::New();

  filter->SetInput(spMovingImg);

  using TransformType = itk::AffineTransform<double, 3>;
  const auto transform = TransformType::New();
  filter->SetTransform(transform);

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, double>;
  const auto interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);
  filter->SetDefaultPixelValue(0);

  // const double outputSpacing[2] = { 1.0, 1.0 };
  // const double outputOrigin[2] = { 0.0, 0.0 };

  UShortImageType::DirectionType direction;
  direction.SetIdentity();
  filter->SetOutputDirection(direction);

  // ResampledImgType2D::SizeType outSize;

  auto movingSpacing = imgSpacing;
  auto movingOrigin = imgOrigin;
  auto movingSize = imgSize;

  switch (enPlane) {
  case enPLANE::PLANE_AXIAL:
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?

    movingSpacing[2] = 1.0;
    movingOrigin[2] = pos;
    movingSize[2] = 1;
    // Resample Here! make corresponding 2D image for Moving image
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);

    break;
  case enPLANE::PLANE_FRONTAL:
    width = imgSize[0];
    height = imgSize[2];
    iCntSlice = imgSize[1];
    iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[1] = 1.0;
    movingOrigin[1] = pos;
    movingSize[1] = 1;

    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[2]);
    break;
  case enPLANE::PLANE_SAGITTAL:
    width = imgSize[1];
    height = imgSize[2];
    iCntSlice = imgSize[0];
    iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
    it.SetFirstDirection(1);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[0] = 1.0;
    movingOrigin[0] = pos;
    movingSize[0] = 1;

    YKFixed.SetSpacing(imgSpacing[1], imgSpacing[2]);
    break;

  default:
    std::cout << "default should not passed by" << std::endl;
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  }

  filter->SetOutputSpacing(movingSpacing);
  filter->SetOutputOrigin(movingOrigin);
  filter->SetSize(movingSize);
  filter->Update();

  YKFixed.CreateImage(width, height, 0);
  // std::cout << "Before MovingImg Creation " << std::endl;

  YKMoving.CreateImage(width, height, 0); // exactly same dimension

  itk::ImageRegionConstIterator<UShortImageType> itMoving(
      filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
  size_t cnt = 0;

  // this simple code will cause flip of the image in frontal and sagittal image
  for (itMoving.GoToBegin(); !itMoving.IsAtEnd(); ++itMoving) {
    YKMoving.m_pData[cnt] = itMoving.Get();
    ++cnt;
  }
  if (enPlane != enPLANE::PLANE_AXIAL) {
    YKMoving.EditImage_Flip();
  }

  it.GoToBegin();

  auto iNumSlice = 0;

  if (iReqSlice < 0 || iReqSlice >= static_cast<int>(iCntSlice)) {
    return;
  }

  while (!it.IsAtEnd()) {
    if (iNumSlice == iReqSlice) {
      size_t iNumHeight = 0;

      while (!it.IsAtEndOfSlice()) {
        size_t iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          const auto fixedImgVal = it.Get();

          if (enPlane == enPLANE::PLANE_AXIAL) {
            YKFixed.m_pData[iNumWidth + width * iNumHeight] = fixedImgVal;
          } else {
            YKFixed.m_pData[iNumWidth + width * (height - iNumHeight - 1)] =
                fixedImgVal;
          }

          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      break;
    }
    it.NextSlice();
    iNumSlice++;
  }

  //    std::cout << "YK Images were filled" << std::endl;
}

void CbctRecon::Draw2DFrom3D(UShortImageType::Pointer &pImg,
                             const enPLANE direction, const double pos,
                             YK16GrayImage &Output2D) const {
  if (pImg == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      pImg, pImg->GetBufferedRegion());

  auto imgSize = pImg->GetBufferedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::SizeType imgSizeBuf = pImg->GetBufferedRegion().GetSize();
  // //1016x1016 x z UShortImageType::SizeType imgSizeLargest =
  // pImg->GetLargestPossibleRegion().GetSize(); //1016x1016 x z

  auto imgOrigin = pImg->GetOrigin();
  auto imgSpacing = pImg->GetSpacing();

  auto width = imgSize[0];
  auto height = imgSize[1];
  auto iCntSlice = imgSize[2];
  auto iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);

  switch (direction) {
  case enPLANE::PLANE_AXIAL:
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    Output2D.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  case enPLANE::PLANE_FRONTAL:
    width = imgSize[0];
    height = imgSize[2];
    iCntSlice = imgSize[1];
    iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(2); // y?
    Output2D.SetSpacing(imgSpacing[0], imgSpacing[2]);
    break;
  case enPLANE::PLANE_SAGITTAL:
    width = imgSize[1];
    height = imgSize[2];
    iCntSlice = imgSize[0];
    iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
    it.SetFirstDirection(1);  // x?
    it.SetSecondDirection(2); // y?
    Output2D.SetSpacing(imgSpacing[1], imgSpacing[2]);
    break;
  default:
    std::cout << "default should not be passed by" << std::endl;
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    Output2D.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  }

  Output2D.CreateImage(width, height, 0);
  it.GoToBegin();

  auto iNumSlice = 0;

  if (iReqSlice < 0 || iReqSlice >= static_cast<int>(iCntSlice)) {
    return;
  }

  while (!it.IsAtEnd()) {
    if (iNumSlice == iReqSlice) {
      auto iNumHeight = 0;
      while (!it.IsAtEndOfSlice()) {
        auto iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          // double tmpVal = it.Get()*multiplyFactor;
          const double tmpVal = it.Get();

          if (direction == enPLANE::PLANE_AXIAL) {
            Output2D.m_pData[iNumWidth + width * iNumHeight] = tmpVal;
          } else {
            Output2D.m_pData[iNumWidth + width * (height - iNumHeight - 1)] =
                tmpVal;
          }

          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      break;
    }
    it.NextSlice();
    iNumSlice++;
  }
}

void CbctRecon::RegisterImgDuplication(const enREGI_IMAGES src,
                                       const enREGI_IMAGES target) {

  UShortImageType::Pointer tmpSrc;

  switch (src) {
  case enREGI_IMAGES::REGISTER_REF_CT:
    tmpSrc = m_spRefCTImg;
    break;
  default:
    std::cerr << "You are using a non-valid target!" << std::endl;
    return;
  }

  if (tmpSrc == nullptr) {
    std::cout << "src image is empty" << std::endl;
    return;
  }

  // Duplication for registration. Starting point is manual Rigid CT image
  using DuplicatorType = itk::ImageDuplicator<UShortImageType>;
  auto duplicator = DuplicatorType::New();
  duplicator->SetInputImage(tmpSrc);
  duplicator->Update();

  switch (target) {
  case enREGI_IMAGES::REGISTER_MANUAL_RIGID:
    m_spManualRigidCT = duplicator->GetOutput();
    break;
  default:
    std::cerr << "You are using a non-valid target!" << std::endl;
  }
  // Duplication for : End
}

template <size_t N>
constexpr bool any_of_str_in_str(const std::array<std::string, N> &test_strings,
                                 const std::string &string_to_search) {
  for (auto &&test_str : test_strings) {
    if (string_to_search.find(test_str) != std::string::npos) {
      return true;
    }
  }
  return false;
}

void CbctRecon::FindAllRelevantPaths(
    const fs::path &pathProjHisDir) // called following SLT_SetHisDir
{
  // in case of eletka, img_UID
  // std::string aa;
  // std::cout<< "ddd " << aa.toLocal8Bit().constData() << std::endl;

  m_strDCMUID.clear();
  m_strPathPatientDir.clear();
  m_strPatientDirName.clear();
  m_strPathFRAME_DBF.clear();
  m_strPathIMAGE_DBF.clear();
  m_strPathGeomXML.clear();
  m_strPathPlanCTDir.clear();
  m_strPathRS.clear();
  m_strPathRS_CBCT.clear();
  m_strPathElektaINI.clear();
  m_strPathElektaINIXVI2.clear();
  m_strPathPlan.clear();
  m_strPathIMAGES.clear();

  auto curHisDir = fs::path(pathProjHisDir);
  auto movingDir = fs::path(pathProjHisDir);
  m_projFormat = enProjFormat::HIS_FORMAT;

  const auto cur_his_path_str = curHisDir.string();

  std::array<std::string, 4> dir_prefixes{{"img_", "fwd_", "sca_", "cor_"}};

  if (!any_of_str_in_str(dir_prefixes, cur_his_path_str)) {
    if (cur_his_path_str.find("Scan0") != std::string::npos) {
      std::cout << "XML set by guessing: Scan0/../ProjectionInfo.xml"
                << std::endl;
      m_strPathGeomXML =
          fs::absolute(curHisDir).parent_path() / "ProjectionInfo.xml";
      std::cout << "Patient DIR set to: Scan0/../../" << std::endl;
      m_projFormat = enProjFormat::HND_FORMAT;
      return;
    }
    if (fs::absolute(curHisDir).string().find("Acquisitions")) {
      std::cout << "XML set by guessing: Acquisitions/../Scan.xml" << std::endl;
      m_strPathGeomXML =
          fs::absolute(curHisDir).parent_path().parent_path() / "Scan.xml";
      std::cout << "Patient DIR set to: Acquisitions/../../../" << std::endl;
      m_projFormat = enProjFormat::XIM_FORMAT;
      return;
    }

    std::cout << "Projection folder should have format [img_UID]" << std::endl;
    std::cout << "XML file cannot be made" << std::endl;
    return;
  }

  const auto underscore_pos = cur_his_path_str.find("_");
  m_strDCMUID = cur_his_path_str.substr(underscore_pos);

  // m_strDCMUID = cur_his_path.right(cur_his_path.length() - 4);

  if (!movingDir.has_parent_path()) // projDir ==> IMAGES
  {
    std::cout << "no upper dir" << std::endl;
    return;
  }
  m_strPathIMAGES = fs::absolute(movingDir);

  if (!movingDir.parent_path()
           .has_parent_path()) // IMAGES ==> patient_402-02-78
  {
    std::cout << "no upper dir" << std::endl;
    return;
  }
  auto tmpDir_PatientFolder{
      fs::absolute(movingDir.parent_path().parent_path())};

  if (!movingDir.parent_path()
           .parent_path()
           .has_parent_path()) // patient_402-02-78 ==> Data folder where DBF
                               // files are.
  {
    std::cout << "no upper dir" << std::endl;
    return;
  }

  m_strPatientDirName = tmpDir_PatientFolder.parent_path().filename();
  m_strPathPatientDir = fs::absolute(tmpDir_PatientFolder.parent_path());

  if (tmpDir_PatientFolder.has_parent_path()) {
    m_strPathDirDefault = m_strPathPatientDir;
  }

  // option 1: already made rtk xml file
  const auto tmpPathRTKGeometry =
      m_strPathPatientDir / std::string("ElektaGeom_" + m_strDCMUID + ".xml");

  // option 2
  const auto pathXVIGeometryXML = fs::absolute(curHisDir) / "_Frames.xml";

  if (fs::exists(tmpPathRTKGeometry)) // The best option: rtk geometry file is
                                      // already existing
  {
    std::cout << "RTK XLM file is found" << std::endl;
    m_strPathGeomXML = tmpPathRTKGeometry;
  } else if (fs::exists(pathXVIGeometryXML)) // 2nd option:_Frames.xml already
                                             // exists in each projection folder
                                             // for > XVI5.0.2
  {
    std::cout << "XVI XLM file is found" << std::endl;
    // YKdebug: Later, it
    m_strPathGeomXML = pathXVIGeometryXML;

    // Later, it should be

    // XVIXMLReader.SetFile(pathXVIGeometryXML)
    // XVIXMLReader.GenerateData()

    // rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter
    // xmlWriter->SetFilename( tmpPathRTKGeometry ) //as option 1
    // xmlWriter->SetObject(XVIXMLReader->GetGeometry());
    // TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
    // m_strPathGeomXML = tmpPathRTKGeometry;
  } else {
    auto tmpStrPath1 = m_strPathPatientDir;
    auto tmpStrPath2 = m_strPathPatientDir;

    // 1st priority: DBF files saved in each "patient" folder --> just in case
    // data are collected separately
    auto fInfo_FrameDBF = tmpStrPath1 / "FRAME.DBF";
    auto fInfo_ImageDBF = tmpStrPath2 / "IMAGE.DBF";

    if (!fs::exists(fInfo_FrameDBF) || !fs::exists(fInfo_ImageDBF)) {
      std::cout << "No found in the patient folder. DBF files can be saved in "
                   "each individual patient as well. Continues to search them "
                   "again in root folder(standard)"
                << std::endl;

      fInfo_FrameDBF = m_strPathPatientDir / "FRAME.DBF";
      fInfo_ImageDBF = m_strPathPatientDir / "IMAGE.DBF";

      if (!fs::exists(fInfo_FrameDBF) || !fs::exists(fInfo_ImageDBF)) {
        std::cout << "DBF files were not found" << std::endl;
        std::cout << "XML file cannot be made" << std::endl;
        return;
      }
    } else {
      std::cout << "DBF files are found in the individual patient directory."
                << std::endl;
    }
    m_strPathFRAME_DBF = fs::absolute(fInfo_FrameDBF);
    m_strPathIMAGE_DBF = fs::absolute(fInfo_ImageDBF);

    m_strPathGeomXML.clear();
    m_strPathGeomXML = crl::MakeElektaXML(
        m_strPathIMAGE_DBF.string(), m_strPathFRAME_DBF.string(),
        m_strDCMUID); // if DBF files exist but UID is not
                      // found, it will crash

    if (m_strPathGeomXML.empty()) {
      std::cout << "No releated data in DBF file" << std::endl;
      return;
    }
  }

  // std::cout << "Root folder: " <<
  // tmpDir_RootFolder.absolutePath().toLocal8Bit().constData() << std::endl;
  // std::cout << "Root folder2: " <<
  // tmpDir_RootFolder.absolutePath().append("\\FRAME.DBF").toLocal8Bit().constData()
  // << std::endl;  std::cout <<
  // fInfo_FrameDBF.absoluteFilePath().toLocal8Bit().constData() << std::endl;

  // GenerateXMLFunc.

  // Search for the geometry XML file: naming convention:
  // ElektaGeom_DICOMUID.xml  after Generation of the XML from DBF files
  //  movingDir = tmpDir_IMAGES;
  // tmpDir_PatientFolder;

  auto enDirStructure_Type = 2;
  // 0: in patient DIR --> 3 folders(CT_SET, DICOM_PLAN, IMAGES)
  // 1: // Patient DIR ==> IMAGES --> CT_SET / DICOM_PLAN
  // 2: NO CT image

  auto tmpDIR_CTSET = fs::absolute(tmpDir_PatientFolder / "CT_SET");

  if (fs::exists(tmpDIR_CTSET)) {
    enDirStructure_Type = 0;
  } else {
    auto tmpStrPathCTSET = m_strPathIMAGES;
    tmpDIR_CTSET = tmpStrPathCTSET / "CT_SET";

    if (fs::exists(tmpDIR_CTSET)) {
      enDirStructure_Type = 1;
    }
  }

  // std::string strPathCTSet = m_strPathIMAGES.append("/CT_SET");

  // switch (enDirStructure_Type)
  // {
  // case 0:
  //  movingDir = tmpDir_IMAGES;
  // break;
  // case 1:
  // movingDir = tmpDir_IMAGES;
  // break;
  // case 2:
  // std::cout << "No CT DICOM folder exist. Proceeding w/o CT" << std::endl;
  // break;
  // }
  //

  if (enDirStructure_Type != 2) {
    for (auto &listDir : fs::directory_iterator(tmpDIR_CTSET)) {
      if (listDir.is_directory()) {
        m_strPathPlanCTDir = listDir.path();
      }
    }
  }

  // for (int i = 0 ; i < listDir.size() ; i++)
  // {
  ////std::cout << listDir.at(i).absolutePath().toLocal8Bit().constData() <<
  /// std::endl; //this returns Dir, not itself
  // std::string tmpPath = listDir.at(i).absoluteFilePath();
  // }

  for (auto &listFile : fs::directory_iterator(tmpDIR_CTSET)) {
    if (listFile.is_regular_file()) {
      auto ext = listFile.path().extension().string();
      if (listFile.path().stem().string().find("RS") != std::string::npos &&
          (ext.compare("DCM") == 0 || ext.compare("dcm") == 0)) {
        m_strPathRS = listFile.path();
      }
    }
  }

  auto tmpDIR_DCM_Plan = fs::absolute(tmpDir_PatientFolder) / "DICOM_PLAN";

  if (fs::exists(tmpDIR_DCM_Plan)) {
    enDirStructure_Type = 0;
  } else {
    auto tmpStrPathCTSET = m_strPathIMAGES;
    tmpDIR_DCM_Plan = tmpStrPathCTSET / "DICOM_PLAN";

    if (fs::exists(tmpDIR_DCM_Plan)) {
      enDirStructure_Type = 1;
    } else {
      enDirStructure_Type = 2;
    }
  }

  if (enDirStructure_Type != 2) {
    for (auto &listFileDCMPlan : fs::directory_iterator(tmpDIR_DCM_Plan)) {
      if (listFileDCMPlan.is_regular_file()) {
        if (listFileDCMPlan.path().extension().compare("DCM") == 0 ||
            listFileDCMPlan.path().extension().compare("dcm") == 0) {
          m_strPathPlan = fs::absolute(listFileDCMPlan.path());
          break; // get first one only
        }
      }
    }
  }

  fs::path movingDirCBCTRS;

  if (enDirStructure_Type == 0) {
    movingDirCBCTRS = std::move(tmpDir_PatientFolder);
  } else if (enDirStructure_Type == 1) {
    movingDirCBCTRS = m_strPathIMAGES;
  }

  movingDirCBCTRS /= "CBCT_RS";

  if (!fs::exists(movingDirCBCTRS)) {
    std::cout << "no CBCT_RS dir exists. Proceed with out CBCT RS image"
              << std::endl;
  } else {
    for (auto &listFile : fs::directory_iterator(movingDirCBCTRS)) {
      if (listFile.is_regular_file() &&
          (listFile.path().extension().compare("DCM") == 0 ||
           listFile.path().extension().compare("DCM") == 0)) {
        m_strPathRS_CBCT = fs::absolute(listFile.path());
        break;
      }
    }
  }

  const auto pathAcqParamDir = pathProjHisDir / "Reconstruction";

  if (fs::exists(pathAcqParamDir)) {

    auto iMinNameLength = 9999;

    auto iMaxNameLength = 0;
    auto iCnt_INIXVI = 0;

    fs::path strPathINIXVI_long;
    for (const auto &listFileAcqParam :
         fs::directory_iterator(pathAcqParamDir)) {
      // suffix:*.tar.gz ==> gz only
      if (listFileAcqParam.is_regular_file() &&
          (listFileAcqParam.path().extension().compare("INI") == 0 ||
           listFileAcqParam.path().extension().compare("ini") == 0)) {
        auto tmpPath = fs::absolute(listFileAcqParam.path());

        if (tmpPath.string().size() < iMinNameLength) {
          iMinNameLength = tmpPath.string().size();
          m_strPathElektaINI = tmpPath;
        }
      }

      auto StrSuffix = listFileAcqParam.path().stem().extension().string() +
                       listFileAcqParam.path().extension().string();

      if (StrSuffix.compare("INI.XVI") == 0 ||
          StrSuffix.compare("ini.xvi") == 0) {
        iCnt_INIXVI++;

        auto tmpPath2 = fs::absolute(listFileAcqParam.path());

        if (tmpPath2.string().size() > iMaxNameLength) {
          iMaxNameLength = tmpPath2.string().size();
          strPathINIXVI_long = tmpPath2;
        }
      }
    }

    if (iCnt_INIXVI == 2) {
      m_strPathElektaINIXVI2 = strPathINIXVI_long;
    }
  }

  std::cerr << "m_strDCMUID: " << m_strDCMUID << "\n"
            << "m_strPathPatientDir: " << m_strPathPatientDir << "\n"
            << "m_strPatientDirName: " << m_strPatientDirName << "\n"
            << "m_strPathFRAME_DBF: " << m_strPathFRAME_DBF << "\n"
            << "m_strPathIMAGE_DBF: " << m_strPathIMAGE_DBF << "\n"
            << "m_strPathGeomXML: " << m_strPathGeomXML << "\n"
            << "m_strPathPlanCTDir: " << m_strPathPlanCTDir << "\n"
            << "m_strPathRS: " << m_strPathRS << "\n"
            << "m_strPathRS_CBCT: " << m_strPathRS_CBCT << "\n"
            << "m_strPathPlan: " << m_strPathPlan << "\n"
            << "m_strPathElektaINI: " << m_strPathElektaINI << "\n"
            << "m_strPathElektaINIXVI2: " << m_strPathElektaINIXVI2 << "\n";
}

void CbctRecon::SaveProjImageAsHIS(FloatImageType::Pointer &spProj3D,
                                   std::vector<YK16GrayImage> arrYKImage,
                                   const fs::path &strSavingFolder,
                                   const double resampleF) const {
  std::cout << "Starting Saving files" << std::endl;

  FloatImageType::Pointer targetImg3D;
  const auto restoreResampleF = 1.0 / resampleF;

  if (fabs(resampleF - 1.0) > 0.001) {
    std::cout << "restore the  resampled image by applying a factor of "
              << restoreResampleF << std::endl;
    ResampleItkImage(spProj3D, targetImg3D, restoreResampleF);
  } else {
    targetImg3D = spProj3D;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it_FwdProj(
      targetImg3D, targetImg3D->GetBufferedRegion());

  it_FwdProj.SetFirstDirection(0);
  it_FwdProj.SetSecondDirection(1);
  it_FwdProj.GoToBegin();

  for (auto it = arrYKImage.begin();
       it != arrYKImage.end() && !it_FwdProj.IsAtEnd(); ++it) {

    auto crntFileName = it->m_strFilePath.filename();
    auto crntPath = strSavingFolder / crntFileName;

    {
      std::ofstream fd(crntPath, std::ios::binary);
      fd.write(it->m_pElektaHisHeader, sizeof(it->m_pElektaHisHeader));
      // fwrite(it->m_pElektaHisHeader, 100, 1, fd);
      // this buffer only include header info

      // Search matching slice using slice iterator for m_spProjCTImg
      while (!it_FwdProj.IsAtEndOfSlice()) {
        while (!it_FwdProj.IsAtEndOfLine()) {
          auto tmpVal = static_cast<unsigned short>(it_FwdProj.Get());
          tmpVal = 65535 - tmpVal; // inverse is done here

          fd << tmpVal;
          ++it_FwdProj;
        }
        it_FwdProj.NextLine();
      }
    } // fd.close
    it_FwdProj.NextSlice();
  }

  std::cout << "Saving completed" << std::endl;
}

// spProjRaw3D: raw intensity value (0-65535), spProjCT3D: raw intensity value
// (0-65535)
void CbctRecon::GenScatterMap_PriorCT(FloatImageType::Pointer &spProjRaw3D,
                                      FloatImageType::Pointer &spProjCT3D,
                                      FloatImageType::Pointer &spProjScat3D,
                                      double medianRadius,
                                      const double gaussianSigma,
                                      const bool bSave) const {
  // Scatter map: should be 2D to use 2D median, Gaussian filters
  if (m_iCntSelectedProj < 1) {
    std::cout << "error: no count of proj image" << std::endl;
    return;
  }

  if (spProjRaw3D == nullptr || spProjCT3D == nullptr) {
    std::cout << "error: proj image 3D is not ready" << std::endl;
    return;
  }

  using SizeType = FloatImageType::SizeType;
  auto size1 = spProjRaw3D->GetBufferedRegion().GetSize();
  auto size2 = spProjCT3D->GetBufferedRegion().GetSize();

  std::cout << "Raw3DProj Size= " << size1 << std::endl;
  std::cout << "spProjCT Size= " << size2 << std::endl;

  auto bHighResolMacro = false; // raw imag= 1024, scattermap = 512
  if (size1[0] != size2[0] || size1[1] != size2[1] || size1[2] != size2[2]) {
    std::cout << "Raw and CT projection dimension are not matching. under the "
                 "high resolution macro?"
              << std::endl;
    // Why 2.0 shouldn't it be 1/downresfactor ?
    if (size1[0] ==
            static_cast<SizeType::SizeValueType>(qRound(size2[0] * 2.0)) &&
        size1[1] ==
            static_cast<SizeType::SizeValueType>(qRound(size2[1] * 2.0)) &&
        size1[2] == size2[2]) {
      bHighResolMacro = true;
    } else {
      return;
    }
  }

  FloatImageType::Pointer spTmpProjRaw3D;

  if (bHighResolMacro) {
    std::cout << "bHighResolMacro is unexpectedly on" << std::endl;
    ResampleItkImage(spProjRaw3D, spTmpProjRaw3D, 0.5);
  } else {
    spTmpProjRaw3D = spProjRaw3D;
  }

  crl::AllocateByRef<FloatImageType, FloatImageType>(spTmpProjRaw3D,
                                                     spProjScat3D);

  auto imgSize = spTmpProjRaw3D->GetBufferedRegion().GetSize();

  // dimension of the spProjRaw3D
  const int iSizeZ = imgSize[2];
  using ImageType = FloatImage2DType;

  for (auto i = 0; i < iSizeZ; i++) {
    ImageType::Pointer spImg2DRaw;
    ImageType::Pointer spImg2DPrim;

    crl::Get2DFrom3D(spTmpProjRaw3D, spImg2DRaw, i,
                     enPLANE::PLANE_AXIAL); // simple conversion between ushort
                                            // 3D to float 2D (using casting,
                                            // not log): input/output: 0-65535
    crl::Get2DFrom3D(spProjCT3D, spImg2DPrim, i, enPLANE::PLANE_AXIAL);

    // The OpenCL version: ~49ms CPU ~76ms
#ifndef _WIN32
    ImageType::Pointer spImg2DScat = OpenCL_LogItoI_subtract_median_gaussian(
        spImg2DRaw, spImg2DPrim, medianRadius, gaussianSigma);
#else

    // CPU version if OpenCL is causing problems:
    using convert_filter_type =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     crl::LineInt2Intensity>;
    auto convert_filter = convert_filter_type::New();
    convert_filter->SetInput(spImg2DRaw);

    auto convert_filter_2 = convert_filter_type::New();
    convert_filter_2->SetInput(spImg2DPrim);

    auto subtract_filter =
        itk::SubtractImageFilter<ImageType, ImageType>::New();
    subtract_filter->SetInput1(convert_filter->GetOutput());
    subtract_filter->SetInput2(convert_filter_2->GetOutput());

    using MedianFilterType = itk::MedianImageFilter<ImageType, ImageType>;
    auto medianFilterY = MedianFilterType::New();
    MedianFilterType::InputSizeType radiusY{};
    radiusY[0] = 0;
    radiusY[1] = medianRadius;
    medianFilterY->SetRadius(radiusY);
    medianFilterY->SetInput(subtract_filter->GetOutput());

    auto medianFilterX = MedianFilterType::New();
    MedianFilterType::InputSizeType radiusX{};
    radiusX[0] = medianRadius;
    radiusX[1] = 0;
    medianFilterX->SetRadius(radiusX);
    medianFilterX->SetInput(medianFilterY->GetOutput());

    auto gaussian_filter =
        itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>::New();
    gaussian_filter->SetInput(medianFilterX->GetOutput());

    itk::SmoothingRecursiveGaussianImageFilter<
        ImageType, ImageType>::SigmaArrayType gauss_sigma;
    const auto sca_size = spImg2DRaw->GetBufferedRegion().GetSize();
    if (sca_size[0] == sca_size[1]) {
      gauss_sigma[0] = gaussianSigma;
      gauss_sigma[1] = gauss_sigma[0] * 0.75;
      gaussian_filter->SetSigmaArray(gauss_sigma);
    } else {
      gaussian_filter->SetSigma(gaussianSigma);
    }
    /* Converting back seems like a waste of compute power and may cause
    problems with log(<0) auto convert_back_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     Intensity2LineInt>::New();
    convert_back_filter->SetInput(gaussian_filter->GetOutput());
    convert_back_filter->Update();
        */
    gaussian_filter->Update();
    ImageType::Pointer spImg2DScat = gaussian_filter->GetOutput();
#endif

    // float to unsigned short
    crl::Set2DTo3D<FloatImageType>(
        spImg2DScat, spProjScat3D, i,
        enPLANE::PLANE_AXIAL); // input/Output: 0-65535 intensity valuesno
                               // mu_t to intensity converion is involved

    const auto unit = qRound(iSizeZ / 10.0);
    if (i % unit == 0) {
      std::cout << "Generating scatter map: "
                << i / static_cast<double>(unit) * 10.0 << " % is done"
                << std::endl;
    }
  } // end of for

  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << "Patient DIR Path: " << m_strPathPatientDir << std::endl;

    if (fs::is_empty(m_strPathPatientDir) &&
        m_projFormat == enProjFormat::HIS_FORMAT) {
      std::cout << "File save error!: No patient DIR name" << std::endl;
      return;
    }

    // Get current folder
    fs::path crntDir;
    if (m_projFormat == enProjFormat::HIS_FORMAT) {
      crntDir = m_strPathPatientDir / "IMAGES"; // current Proj folder
    } else {
      // stolen from registration class: m_strPathPlastimatch definition
      auto crntDir =
          fs::current_path(); // folder where current exe file exists.
      auto crntPathStr = fs::absolute(crntDir);
      auto dirName = crntPathStr / "plm_tmp";

      if (!fs::exists(dirName)) {
        if (!fs::create_directory(dirName)) {
          std::cerr << "Could not create " << dirName << std::endl;
          // Not enough reason to fail
        }
      }
      crntDir = dirName;
    }

    // Make a sub directory

    if (!fs::exists(crntDir)) {
      std::cout << "File save error: The specified folder does not exist."
                << std::endl;
      return;
    }

    const auto scatDirName = "sca_" + m_strDCMUID;

    if (!fs::create_directory(crntDir / scatDirName)) {
      std::cout << "Scatter map directory seems to exist already. Files will "
                   "be overwritten."
                << std::endl;
    }

    auto strSavingFolder = crntDir / scatDirName;
    if (m_projFormat == enProjFormat::HIS_FORMAT) {
      SaveProjImageAsHIS(spProjScat3D, m_arrYKBufProj, strSavingFolder,
                         m_fResampleF);
      std::cerr << "Scatter saved in Intensity values as HIS at: "
                << strSavingFolder << "\n";
    } else {
      using imagewritertype = itk::ImageFileWriter<FloatImageType>;
      auto imagewriter = imagewritertype::New();
      imagewriter->SetInput(spProjScat3D);
      imagewriter->SetFileName((strSavingFolder / "scatter.mha").string());
      imagewriter->Update();
      std::cerr << "Scatter saved in Intensity values as mha at: "
                << strSavingFolder << "\n";
    }
  }
}

void CbctRecon::ScatterCorr_PrioriCT(FloatImageType::Pointer &spProjRaw3D,
                                     FloatImageType::Pointer &spProjScat3D,
                                     FloatImageType::Pointer &m_spProjCorr3D,
                                     int postMedian, const bool bSave) const {
  // Scatter map: should be 2D to use 2D median, Gaussian filters
  if (m_iCntSelectedProj < 1) {
    std::cout << "error: no count of proj image" << std::endl;
    return;
  }

  if (spProjRaw3D == nullptr || spProjScat3D == nullptr) {
    std::cout << "Error: proj image 3D is not ready" << std::endl;
    return;
  }

  auto size1 = spProjRaw3D->GetBufferedRegion().GetSize();
  auto size2 = spProjScat3D->GetBufferedRegion().GetSize();

  std::cout << "Raw3DProj Size= " << size1 << std::endl;
  std::cout << "spProjScat3D Size= " << size2 << std::endl;

  auto bHighResolMacro = false;

  if (size1[0] != size2[0] || size1[1] != size2[1] || size1[2] != size2[2]) {
    std::cout << "Raw and scatter projection dimension are not matching. under "
                 "the high resolution macro?"
              << std::endl;

    if (static_cast<int>(size1[0]) == qRound(size2[0] * 2.0) &&
        static_cast<int>(size1[1]) == qRound(size2[1] * 2.0) &&
        size1[2] == size2[2]) {
      bHighResolMacro = true;
    } else {
      return;
    }
  }

  FloatImageType::Pointer spTmpProjScat3D;

  if (bHighResolMacro) {
    ResampleItkImage(spProjScat3D, spTmpProjScat3D, 2.0);
  } else {
    spTmpProjScat3D = spProjScat3D;
  }

  crl::AllocateByRef<FloatImageType, FloatImageType>(spProjRaw3D,
                                                     m_spProjCorr3D);

  auto imgSize = spProjRaw3D->GetBufferedRegion().GetSize();

  const int iSizeZ = imgSize[2];

  // std::cout << "resample factor " << resF2D << std::endl;.

  for (auto i = 0; i < iSizeZ; i++) {
    FloatImage2DType::Pointer spImg2DRaw;
    FloatImage2DType::Pointer spImg2DScat;

    crl::Get2DFrom3D(spProjRaw3D, spImg2DRaw, i, enPLANE::PLANE_AXIAL);
    crl::Get2DFrom3D(spTmpProjScat3D, spImg2DScat, i, enPLANE::PLANE_AXIAL);

    if (bHighResolMacro) {
      postMedian = postMedian * 2;
    }
    if (postMedian >= 2) {
      postMedian = qRound(postMedian / 2.0);
    }

#ifndef _WIN32
    auto spImg2DCorr = OpenCL_LogItoI_subtract_median_ItoLogI(
        spImg2DRaw, spImg2DScat, postMedian);
#else
    using ImageType = FloatImage2DType;
    auto convert_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     crl::LineInt2Intensity>::New();
    convert_filter->SetInput(spImg2DRaw);
    /* The conversion seems redundant, see comment in GenScatterMap
        auto convert_filter_2 =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     LineInt2Intensity>::New();
    convert_filter_2->SetInput(spImg2DScat);*/

    auto subtract_filter =
        itk::SubtractImageFilter<ImageType, ImageType>::New();
    subtract_filter->SetInput1(convert_filter->GetOutput());
    subtract_filter->SetInput2(spImg2DScat); // convert_filter_2->GetOutput());

    auto median_filter = itk::MedianImageFilter<ImageType, ImageType>::New();
    median_filter->SetInput(subtract_filter->GetOutput());
    median_filter->SetRadius(postMedian);

    auto convert_back_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     crl::Intensity2LineInt>::New();
    convert_back_filter->SetInput(median_filter->GetOutput());
    convert_back_filter->Update();
    ImageType::Pointer spImg2DCorr = convert_back_filter->GetOutput();
#endif
    /*
    auto subtract_filter =
        itk::SubtractImageFilter<FloatImage2DType, FloatImage2DType,
                                 FloatImage2DType>::New();
    subtract_filter->SetInput1(spImg2DRaw);
    subtract_filter->SetInput2(spImg2DScat);
    subtract_filter->Update();
    spImg2DCorr = subtract_filter->GetOutput();

    // Post Median filtering


    if (postMedian >= 2) // YK2015
    {
      using MedianFilterType =
          itk::MedianImageFilter<FloatImage2DType, FloatImage2DType>;
      auto medianFilter = MedianFilterType::New();
      MedianFilterType::InputSizeType radius{};

      radius[0] = qRound(postMedian / 2.0);
      radius[1] = radius[0];

      medianFilter->SetRadius(radius);
      medianFilter->SetInput(spImg2DCorr);
      medianFilter->Update();
      spImg2DCorr = medianFilter->GetOutput();
    }*/

    crl::Set2DTo3D<FloatImageType>(spImg2DCorr, m_spProjCorr3D, i,
                                   enPLANE::PLANE_AXIAL); // float2D to USHORT

    const auto unit = qRound(iSizeZ / 10.0);
    if (i % unit == 0) {
      std::cout << "Applying scatter correction: "
                << i / static_cast<double>(unit) * 10.0 << " % is done"
                << std::endl;
    }
  }

  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << "Patient DIR Path: " << m_strPathPatientDir << std::endl;

    if (fs::is_empty(m_strPathPatientDir)) {
      std::cout << "File save error!: No patient DIR name" << std::endl;
      return;
    }

    // Get current folder
    const auto crntDir = m_strPathPatientDir / "IMAGES"; // current Proj folder

    if (!fs::exists(crntDir)) {
      // Probably IMAGES subdir didn't exist, let's just use the selected dir
      // anyway
      if (!fs::exists(crntDir.parent_path())) {
        std::cerr << "File save error: The specified folder does not exist.\n";
        return;
      }
    }

    const auto scatDirName = "cor_" + m_strDCMUID;

    if (!fs::create_directory(crntDir / scatDirName)) {
      std::cout << "Corrected projection directory seems to exist already. "
                   "Files will be overwritten."
                << std::endl;
    }

    auto strSavingFolder = fs::absolute(crntDir) / scatDirName;

    if (m_projFormat == enProjFormat::HIS_FORMAT) {
      if (!bHighResolMacro) {
        SaveProjImageAsHIS(m_spProjCorr3D, m_arrYKBufProj, strSavingFolder,
                           m_fResampleF);
      } else {
        SaveProjImageAsHIS(m_spProjCorr3D, m_arrYKBufProj, strSavingFolder,
                           1.0);
      }
    } else {
      crl::saveImageAsMHA<FloatImageType>(
          m_spProjCorr3D, (strSavingFolder / "corr_proj.mha").string());
    }
  }
  // spProjScat3D->Initialize(); //memory release
}

// it works! new memory will be allocated for spTarImg
void CbctRecon::ResampleItkImage(FloatImageType::Pointer &spSrcImg,
                                 FloatImageType::Pointer &spTarImg,
                                 const double resFactor) const {
  if (spSrcImg == nullptr) {
    return;
  }

  //  std::cout << "original Origin: " << spSrcImg2D->GetOrigin() <<
  //  std::endl;
  using ResampleImageFilterType =
      itk::ResampleImageFilter<FloatImageType, FloatImageType, float>;
  auto resample = ResampleImageFilterType::New();

  resample->SetOutputDirection(spSrcImg->GetDirection());

  using TransformType = itk::AffineTransform<float, 3>;
  auto transform = TransformType::New();

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<FloatImageType, float>;
  const auto interpolator = InterpolatorType::New();
  resample->SetInterpolator(interpolator);
  if ((m_projFormat == enProjFormat::HIS_FORMAT &&
       DEFAULT_ELEKTA_PROJ_HEIGHT ==
           spSrcImg->GetBufferedRegion().GetSize()[1]) ||
      (m_projFormat != enProjFormat::HIS_FORMAT &&
       DEFAULT_VARIAN_PROJ_HEIGHT ==
           spSrcImg->GetBufferedRegion().GetSize()[1])) {
    resample->SetDefaultPixelValue(50);
  } else {
    resample->SetDefaultPixelValue(0);
  }

  auto inputSize = spSrcImg->GetLargestPossibleRegion().GetSize();
  FloatImageType::SizeType outputSize{};
  outputSize[0] = qRound(inputSize[0] * resFactor);
  outputSize[1] = qRound(inputSize[1] * resFactor);
  outputSize[2] = inputSize[2];
  resample->SetSize(outputSize);

  FloatImageType::SpacingType outputSpacing;
  outputSpacing[0] =
      spSrcImg->GetSpacing()[0] *
      (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] =
      spSrcImg->GetSpacing()[1] *
      (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
  outputSpacing[2] = spSrcImg->GetSpacing()[2];
  resample->SetOutputSpacing(outputSpacing);

  const auto outputOrigin = spSrcImg->GetOrigin(); // Float image
  resample->SetOutputOrigin(outputOrigin);

  resample->SetInput(spSrcImg);
  transform->SetIdentity();
  resample->SetTransform(transform);

  resample->Update();

  // resample->GetOutput()->SetOrigin(prevOrigin);
  spTarImg = resample->GetOutput(); // is it copied? or replaced?
}

void CbctRecon::ResampleItkImage(UShortImageType::Pointer &spSrcImg,
                                 UShortImageType::Pointer &spTarImg,
                                 const double resFactor) const {
  if (spSrcImg == nullptr) {
    return;
  }

  //  std::cout << "original Origin: " << spSrcImg2D->GetOrigin() <<
  //  std::endl;
  using ResampleImageFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType, float>;
  auto resample = ResampleImageFilterType::New();

  resample->SetOutputDirection(spSrcImg->GetDirection());

  using TransformType = itk::AffineTransform<float, 3>;
  auto transform = TransformType::New();

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, float>;
  const auto interpolator = InterpolatorType::New();
  resample->SetInterpolator(interpolator);

  resample->SetDefaultPixelValue(50);

  auto inputSize = spSrcImg->GetLargestPossibleRegion().GetSize();
  UShortImageType::SizeType outputSize{};
  outputSize[0] = qRound(inputSize[0] * resFactor);
  outputSize[1] = qRound(inputSize[1] * resFactor);
  outputSize[2] = inputSize[2];
  resample->SetSize(outputSize);

  UShortImageType::SpacingType outputSpacing;
  outputSpacing[0] =
      spSrcImg->GetSpacing()[0] *
      (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] =
      spSrcImg->GetSpacing()[1] *
      (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
  outputSpacing[2] = spSrcImg->GetSpacing()[2];
  resample->SetOutputSpacing(outputSpacing);

  const auto outputOrigin = spSrcImg->GetOrigin(); // Float
                                                   // image
  resample->SetOutputOrigin(outputOrigin);

  resample->SetInput(spSrcImg);
  transform->SetIdentity();
  resample->SetTransform(transform);

  resample->Update();

  // resample->GetOutput()->SetOrigin(prevOrigin);
  spTarImg = resample->GetOutput(); // is it copied? or replaced?
}

void CbctRecon::ResampleItkImage2D(FloatImage2DType::Pointer &spSrcImg2D,
                                   FloatImage2DType::Pointer &spTarImg2D,
                                   const double resFactor) const {
  if (spSrcImg2D == nullptr) {
    std::cout << "ERROR! SrcImage is empty" << std::endl;
    return;
  }

  //  std::cout << "original Origin: " << spSrcImg2D->GetOrigin() <<
  //  std::endl;
  using ResampleImageFilterType =
      itk::ResampleImageFilter<FloatImage2DType, FloatImage2DType, float>;
  auto resample = ResampleImageFilterType::New();

  resample->SetOutputDirection(spSrcImg2D->GetDirection());

  // outputSpacing[2] = spSrcImg2D->GetSpacing()[2];

  // std::cout << "Output spacing: " << outputSpacing << std::endl;

  // FloatImageType2D::Pointer input = spSrcImg2D;
  // FloatImageType2D::PointType prevOrigin = input->GetOrigin(); //-204.6 -
  // 204.6  0

  // input->SetOrigin(outputOrigin);

  /* std::cout << "outputSize " << outputSize << std::endl;
   std::cout << "OutputSpacing " << outputSpacing << std::endl;
   std::cout << "OutputOrigin " << outputOrigin << std::endl;*/

  // Resample the image
  // typedef itk::IdentityTransform<float, 2> TransformType;

  using TransformType = itk::AffineTransform<float, 2>;
  auto transform = TransformType::New();

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<FloatImage2DType, float>;
  const auto interpolator = InterpolatorType::New();
  resample->SetInterpolator(interpolator);

  resample->SetDefaultPixelValue(50);

  auto inputSize = spSrcImg2D->GetLargestPossibleRegion().GetSize();
  FloatImage2DType::SizeType outputSize{};
  outputSize[0] = qRound(inputSize[0] * resFactor);
  outputSize[1] = qRound(inputSize[1] * resFactor);
  resample->SetSize(outputSize);

  FloatImage2DType::SpacingType outputSpacing;
  outputSpacing[0] =
      spSrcImg2D->GetSpacing()[0] *
      (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] =
      spSrcImg2D->GetSpacing()[1] *
      (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
  resample->SetOutputSpacing(outputSpacing);

  const auto outputOrigin = spSrcImg2D->GetOrigin(); // Float image
  resample->SetOutputOrigin(outputOrigin);

  resample->SetInput(spSrcImg2D);
  transform->SetIdentity();
  resample->SetTransform(transform);

  resample->Update();

  // resample->GetOutput()->SetOrigin(prevOrigin);
  spTarImg2D = resample->GetOutput(); // is it copied? or replaced?

  // std::cout << "resampled Origin: " << spTarImg2D->GetOrigin() <<
  // std::endl;
}

void CbctRecon::AfterScatCorrectionMacro(const bool use_cuda,
                                         const bool use_opencl,
                                         const bool save_dicom,
                                         FDK_options &fdk_options) {
  // Original projection file can be replaced by the corrected one
  // Current projection map (float) used for the reconstruction is:
  // m_spProjImg3DFloat and this is resampled one
  m_spProjImg3DFloat =
      m_spProjImgCorr3D; // ConvertIntensity2LineInt(m_spProjImgCorr3D);

  // Do reconstruction

  // Regardeless of previous setting, The Truncation should not be applied!

  // Truncation is invalidated inside the function
  if (use_cuda) {
    DoReconstructionFDK<enDeviceType::CUDA_DEVT>(
        enREGI_IMAGES::REGISTER_COR_CBCT, fdk_options);
  } else if (use_opencl) {
    DoReconstructionFDK<enDeviceType::OPENCL_DEVT>(
        enREGI_IMAGES::REGISTER_COR_CBCT, fdk_options);
  } else {
    DoReconstructionFDK<enDeviceType::CPU_DEVT>(
        enREGI_IMAGES::REGISTER_COR_CBCT, fdk_options);
  }

  // Save Image as DICOM
  if (save_dicom) {
    // Get current folder
    const auto crntDir = m_strPathPatientDir / "IMAGES" /
                         ("cor_" + m_strDCMUID); // current Proj folder
    const auto SubDirName = "Reconstruction"s;
    auto savingFolder = crntDir / SubDirName;
    if (!fs::create_directory(crntDir / SubDirName)) {
      std::cout
          << "DICOM dir seems to exist already. Files will be overwritten."
          << std::endl;
    }
    auto updated_text_ct = "PriorCT_ScatterCorr"s;
    crl::SaveUSHORTAsSHORT_DICOM(m_spScatCorrReconImg, m_strDCMUID,
                                 updated_text_ct, savingFolder);
    // Export as DICOM (using plastimatch) folder?
  }
  std::cout << "Exiting AfterScatCorrectionMacro.";
}

int CbctRecon::CropSkinUsingThreshold(const int threshold,
                                      const int erode_radius,
                                      const int dilate_radius) {
  std::cout << "Overwriting of values below threshold to air ";
  if (m_spCrntReconImg == nullptr) {
    return -1;
  }

  using threshFilterType =
      itk::BinaryThresholdImageFilter<UShortImageType, UShortImageType>;
  auto threshFilter = threshFilterType::New();
  threshFilter->SetInput(m_spCrntReconImg);

  threshFilter->SetOutsideValue(0);
  threshFilter->SetInsideValue(1);
  threshFilter->SetLowerThreshold(threshold);
  threshFilter->Update();
  const UShortImageType::Pointer spCrntImgMask = threshFilter->GetOutput();
  using iteratorType = itk::ImageRegionIteratorWithIndex<UShortImageType>;
  iteratorType it(spCrntImgMask, spCrntImgMask->GetBufferedRegion());
  auto imgDims = spCrntImgMask->GetBufferedRegion().GetSize();

  it.GoToBegin();
  while (!it.IsAtEnd()) {
    const int z_idx = it.GetIndex()[2];
    if (z_idx == static_cast<int>(imgDims[2] - 10) || z_idx == 10) {
      it.Set(1.0f);
    }
    ++it;
  }

  using HoleFillingFilterType = itk::BinaryFillholeImageFilter<UShortImageType>;
  auto HoleFillingFilter = HoleFillingFilterType::New();
  HoleFillingFilter->SetForegroundValue(1);
  HoleFillingFilter->SetFullyConnected(false);
  std::cout << "Threshold filtering.. ";
  HoleFillingFilter->SetInput(spCrntImgMask);

  using StructElementType =
      itk::BinaryBallStructuringElement<USHORT_PixelType, 3>;
  using ErodeFilterType =
      itk::BinaryErodeImageFilter<UShortImageType, UShortImageType,
                                  StructElementType>;
  using DilateFilterType =
      itk::BinaryDilateImageFilter<UShortImageType, UShortImageType,
                                   StructElementType>;
  auto binaryErode = ErodeFilterType::New();
  binaryErode->SetErodeValue(1);
  StructElementType erodeStructElement;
  erodeStructElement.SetRadius(erode_radius);
  erodeStructElement.CreateStructuringElement();
  binaryErode->SetKernel(erodeStructElement);
  std::cout << "filling holes.. ";
  binaryErode->SetInput(HoleFillingFilter->GetOutput());

  auto binaryDilate = DilateFilterType::New();
  binaryDilate->SetDilateValue(1);
  StructElementType dilateStructElement;
  dilateStructElement.SetRadius(dilate_radius);
  dilateStructElement.CreateStructuringElement();
  binaryDilate->SetKernel(dilateStructElement);
  std::cout << "eroding dirt.. ";
  binaryDilate->SetInput(binaryErode->GetOutput());
  std::cout << "Skin mask is being created..." << std::endl;

  using MaskFilterType =
      itk::MaskImageFilter<UShortImageType, UShortImageType, UShortImageType>;
  auto MaskFilter = MaskFilterType::New();
  MaskFilter->SetMaskingValue(0);
  std::cout << "Dilating.. ";
  MaskFilter->SetMaskImage(binaryDilate->GetOutput());
  if (m_spCrntReconImg == m_spRawReconImg) {
    MaskFilter->SetInput(m_spRawReconImg);
    MaskFilter->Update();
    m_spRawReconImg = MaskFilter->GetOutput();
    return 1;
  }
  if (m_spCrntReconImg == m_spRefCTImg) {
    MaskFilter->SetInput(m_spRefCTImg);
    MaskFilter->Update();
    m_spRefCTImg = MaskFilter->GetOutput();
    return 2;
  }
  if (m_spCrntReconImg == m_spScatCorrReconImg) {
    MaskFilter->SetInput(m_spScatCorrReconImg);
    MaskFilter->Update();
    m_spScatCorrReconImg = MaskFilter->GetOutput();
    return 3;
  }
  return 0;
}

// Below version is optimized for many points and much faster
void CbctRecon::ExportAngularWEPL_byFile(fs::path &strPathOutput,
                                         const double fAngleStart,
                                         const double fAngleEnd,
                                         const double fAngleGap) {
  if (strPathOutput.empty()) {
    return;
  }

  if (m_vPOI_DCM.empty()) {
    std::cout << "No POI data is prepared. Load them first" << std::endl;
    return;
  }

  if (m_spRawReconImg == nullptr) {
    std::cout << "Error: no Raw Recon image is found" << std::endl;
    return;
  }

  if (m_spScatCorrReconImg == nullptr) {
    std::cout << "Warning: no ScatCorrReconImg is found" << std::endl;
  }

  if (m_spManualRigidCT == nullptr) {
    std::cout << "Warning: no ManualRigidCT is found" << std::endl;
  }
  if (m_spAutoRigidCT == nullptr) {
    std::cout << "Warning: no AutoRigidCT is found" << std::endl;
  }
  if (m_spDeformedCT_Final == nullptr) {
    std::cout << "Warning: no DeformedCT is found" << std::endl;
  }

  std::vector<WEPLData> vOutputWEPL_manual;
  std::vector<WEPLData> vOutputWEPL_auto_rigid;
  std::vector<WEPLData> vOutputWEPL_deform;
  std::vector<WEPLData> vOutputWEPL_rawCBCT;
  std::vector<WEPLData> vOutputWEPL_corCBCT;
#pragma omp parallel sections
  {
#pragma omp section
    {
      GetAngularWEPL_window(m_spRawReconImg, fAngleGap, fAngleStart, fAngleEnd,
                            vOutputWEPL_rawCBCT,
                            true); // mandatory
      std::cout << "Done: (RAW)";
    }
#pragma omp section
    {
      if (m_spScatCorrReconImg != nullptr) {
        try {
          GetAngularWEPL_window(m_spScatCorrReconImg, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_corCBCT, true);
          std::cout << " (COR)";
        } catch (std::exception &e) {
          std::cout << " (COR) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
#pragma omp section
    {
      if (m_spManualRigidCT != nullptr) {
        try {
          GetAngularWEPL_window(m_spManualRigidCT, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_manual, true);
          std::cout << " (MAN)";
        } catch (std::exception &e) {
          std::cout << " (MAN) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
#pragma omp section
    {
      if (m_spAutoRigidCT != nullptr) {
        try {
          GetAngularWEPL_window(m_spAutoRigidCT, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_auto_rigid, true);
          std::cout << " (AUT)";
        } catch (std::exception &e) {
          std::cout << " (AUT) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
#pragma omp section
    {
      if (m_spDeformedCT_Final != nullptr) {
        try {
          GetAngularWEPL_window(m_spDeformedCT_Final, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_deform, true);
          std::cout << " (DEF)";
        } catch (std::exception &e) {
          std::cout << " (DEF) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
  }
  std::cout << std::endl;

  std::cout << "Saving results...";

  std::ofstream fout;
  fout.open(strPathOutput);

  const auto cntWEPL = vOutputWEPL_rawCBCT.size();

  fout << "Point Index"
       << "\t"
       << "Gantry Angle"
       << "\t"
       << "Sample Number"
       << "\t"
       << "RawCBCT"
       << "\t";

  if (m_spScatCorrReconImg != nullptr &&
      vOutputWEPL_corCBCT.size() == cntWEPL) {
    fout << "CorrCBCT"
         << "\t";
  }
  if (m_spManualRigidCT != nullptr && vOutputWEPL_manual.size() == cntWEPL) {
    fout << "ManualRigidCT"
         << "\t";
  }
  if (m_spAutoRigidCT != nullptr && vOutputWEPL_auto_rigid.size() == cntWEPL) {
    fout << "AutoRigidCT"
         << "\t";
  }
  if (m_spDeformedCT_Final != nullptr && vOutputWEPL_deform.size() == cntWEPL) {
    fout << "DeformedCT"
         << "\t";
  }
  fout << std::endl;

  for (size_t i = 0; i < cntWEPL; i++) {
    const auto cur_rawpoint = vOutputWEPL_rawCBCT.at(i);
    fout << static_cast<int64_t>(cur_rawpoint.ptIndex) << "\t"
         << cur_rawpoint.fGanAngle << "\t" << static_cast<int64_t>(i) << "\t"
         << cur_rawpoint.fWEPL << "\t";

    if (m_spScatCorrReconImg != nullptr &&
        vOutputWEPL_corCBCT.size() == cntWEPL) {
      fout << vOutputWEPL_corCBCT.at(i).fWEPL << "\t";
    }
    if (m_spManualRigidCT != nullptr && vOutputWEPL_manual.size() == cntWEPL) {
      fout << vOutputWEPL_manual.at(i).fWEPL << "\t";
    }
    if (m_spAutoRigidCT != nullptr &&
        vOutputWEPL_auto_rigid.size() == cntWEPL) {
      fout << vOutputWEPL_auto_rigid.at(i).fWEPL << "\t";
    }
    if (m_spDeformedCT_Final != nullptr &&
        vOutputWEPL_deform.size() == cntWEPL) {
      fout << vOutputWEPL_deform.at(i).fWEPL << "\t";
    }

    fout << std::endl;
  }
  fout.close();
  std::cout << "done!" << std::endl;
}

void CbctRecon::GetAngularWEPL_window(UShortImageType::Pointer &spUshortImage,
                                      const float fAngleGap,
                                      const float fAngleStart,
                                      const float fAngleEnd,
                                      std::vector<WEPLData> &vOutputWEPLData,
                                      const bool bAppend) {
  if (spUshortImage == nullptr) {
    return;
  }

  if (fAngleGap <= 0) {
    return;
  }

  if (!bAppend) {
    vOutputWEPLData.clear();
  }

  const auto wepl_image = crl::wepl::ConvertUshort2WeplFloat(spUshortImage);

  const auto fullAngle = fAngleEnd - fAngleStart;
  const auto sizeAngles = qRound(fullAngle / fAngleGap);

  const std::array<double, 3> pixel_size = {{spUshortImage->GetSpacing()[0],
                                             spUshortImage->GetSpacing()[1],
                                             spUshortImage->GetSpacing()[2]}};

  const auto couch = 0.0;

  for (auto i = 0; i < sizeAngles; ++i) {
    const auto gantry = fAngleStart + i * fAngleGap;
    const auto basis = crl::wepl::get_basis_from_angles(gantry, couch);
    size_t loop_idx = 0;
    // int z_slice = -5000; // for progress and debug
    for (auto &poi_it : m_vPOI_DCM) {
      loop_idx++;
      UShortImageType::PointType cur_point;
      cur_point[0] = poi_it.x;
      cur_point[1] = poi_it.y;
      cur_point[2] = poi_it.z;

      UShortImageType::IndexType cur_idx{};
      if (!spUshortImage->TransformPhysicalPointToIndex(cur_point, cur_idx)) {
        // cur_point not in image
        continue;
      }
      /* // Turn on when debugging
      if (cur_idx[2] != z_slice) {
              z_slice = cur_idx[2];
              std::cout << z_slice << ", " << std::endl;
      }
      */
      const std::array<size_t, 3> point_id = {
          {static_cast<size_t>(cur_idx[0]), static_cast<size_t>(cur_idx[1]),
           static_cast<size_t>(cur_idx[2])}};

      WEPLData wepl_data{};
      wepl_data.fWEPL =
          crl::wepl::WEPL_from_point(point_id, basis, pixel_size, wepl_image);
      wepl_data.ptIndex = loop_idx;
      wepl_data.fGanAngle = gantry;

      vOutputWEPLData.push_back(wepl_data);
    }
  }
}

void CbctRecon::GetAngularWEPL_SinglePoint(
    UShortImageType::Pointer &spUshortImage, const float fAngleGap,
    const float fAngleStart, const float fAngleEnd, const VEC3D &calcPt,
    const size_t curPtIdx, std::vector<WEPLData> &vOutputWEPLData,
    const bool bAppend) const {
  if (spUshortImage == nullptr) {
    return;
  }

  if (fAngleGap <= 0) {
    return;
  }

  const auto wepl_image = crl::wepl::ConvertUshort2WeplFloat(spUshortImage);

  const auto fullAngle = fAngleEnd - fAngleStart;
  const auto sizeAngles = static_cast<size_t>(qRound(fullAngle / fAngleGap));

  UShortImageType::PointType calc_point;
  calc_point[0] = calcPt.x;
  calc_point[1] = calcPt.y;
  calc_point[2] = calcPt.z;
  UShortImageType::IndexType calc_idx{};

  if (!spUshortImage->TransformPhysicalPointToIndex(calc_point, calc_idx)) {
    std::cerr << "Point was outside image!" << std::endl;
  }

  // 1) Generate parms according to the angle e.g 360 parms
  const std::array<size_t, 3> isoTarget = {{static_cast<size_t>(calc_idx[0]),
                                            static_cast<size_t>(calc_idx[1]),
                                            static_cast<size_t>(calc_idx[2])}};

  std::cout << "Target: ( " << isoTarget[0] << ", " << isoTarget[1] << ", "
            << isoTarget[2] << " ), sizeAngles: " << sizeAngles << std::endl;

  const std::string stdout_file = "WEPL_stdout.txt";

  std::ofstream ofs(stdout_file); // Open stdout_file for writing
  if (ofs.is_open()) {
    std::cerr
        << "couldn't open file: " << stdout_file << " for writing!\n"
        << "Are you running this app from a folder without write permissions?"
        << std::endl;
    return;
  }
  const std::array<double, 3> pixel_size = {{wepl_image->GetSpacing()[0],
                                             wepl_image->GetSpacing()[1],
                                             wepl_image->GetSpacing()[2]}};

  auto region = wepl_image->GetLargestPossibleRegion();

  std::vector<WEPLData> stArrWEPL(sizeAngles);

  size_t i = 0;
  for (auto &it : stArrWEPL) {
    // YKTEMP Should be updated according to recent update of plastimatch
    const auto curAngle = fAngleStart + i * fAngleGap;

    const auto basis = crl::wepl::get_basis_from_angles(curAngle, 0.0);

    ofs << std::fixed << std::setprecision(3) << curAngle << ", [" << basis[0]
        << ", " << basis[1] << ", " << basis[2] << "]: ";

    it.fWEPL = crl::wepl::WEPL_from_point(isoTarget, basis, pixel_size,
                               wepl_image); // get_rgdepth
    it.fGanAngle = curAngle;
    it.ptIndex = curPtIdx;
    ofs << std::fixed << std::setprecision(5) << it.fWEPL << "\n";
    i++;
  }

  ofs.close();

  if (!bAppend) {
    vOutputWEPLData.clear();
  }
  vOutputWEPLData.insert(vOutputWEPLData.end(), stArrWEPL.begin(),
                         stArrWEPL.end());
  /*
for (int i = 0; i < sizeAngles; i++)
{
  vOutputWEPLData.push_back(stArrWEPL[i]);
}*/

  // delete[] stArrWEPL;
}

void CbctRecon::GeneratePOIData(const bool AnteriorToPosterior,
                                const double table_posY) // it fills m_vPOI_DCM
{
  if (!m_vPOI_DCM.empty()) {
    m_vPOI_DCM.clear();
  }

  auto imgSize = m_spCrntReconImg->GetLargestPossibleRegion().GetSize();
  const VEC3D imgDims = {static_cast<double>(imgSize[0]),
                         static_cast<double>(imgSize[1]),
                         static_cast<double>(imgSize[2])};
  const VEC3D imgSpac{m_spCrntReconImg->GetSpacing()[0],
                      m_spCrntReconImg->GetSpacing()[1],
                      m_spCrntReconImg->GetSpacing()[2]};

  if (AnteriorToPosterior) {
    for (size_t k = 2; k < imgSize[2] - 2; k++) {
      for (size_t i = 2; i < imgSize[0] - 2; i++) {
        VEC3D fPOI = {
            static_cast<double>(i) * imgSpac.x - (imgSpac.x * imgDims.x) / 2.,
            table_posY,
            static_cast<double>(k) * imgSpac.z - (imgSpac.z * imgDims.z) / 2.};
        m_vPOI_DCM.push_back(fPOI);
      }
    }
  } else {
    for (size_t k = 2; k < imgSize[2] - 2; k++) {
      for (size_t j = 2; j < imgSize[1] - 2; j++) {
        VEC3D fPOI = {
            2.,
            static_cast<double>(j) * imgSpac.y - (imgSpac.y * imgDims.y) / 2.,
            static_cast<double>(k) * imgSpac.z - (imgSpac.z * imgDims.z) / 2.};
        m_vPOI_DCM.push_back(fPOI);
      }
    }
  }
  const auto last_point = m_vPOI_DCM.back();
  std::cout << "POI data generated! last value: [" << last_point.x << ", "
            << last_point.y << ", " << last_point.z << "]" << std::endl;
}

void CbctRecon::LoadExternalFloatImage(fs::path &strPath,
                                       const bool bConversion) {
  using ReaderType = itk::ImageFileReader<FloatImageType>;
  auto reader = ReaderType::New();

  // std::string filePath = strPath;

  reader->SetFileName(strPath.string());
  reader->Update();

  FloatImageType::Pointer spCrntImg = reader->GetOutput();

  // Float image
  std::cout << "Float image has been loaded" << std::endl;

  if (bConversion) {
    crl::TransformationRTK2IEC(spCrntImg);
  }

  using AbsImageFilterType =
      itk::AbsImageFilter<FloatImageType, FloatImageType>;
  auto absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(spCrntImg); // 20140206 modified it was a buug

  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  auto multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  auto castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());
  castFilter->Update();
  m_spRawReconImg = castFilter->GetOutput();
}

// Only can be used for m_spRawRecon
void CbctRecon::MedianFilterByGUI(
    const UShortImageType::SizeType &indexRadius) {
  using FilterType = itk::MedianImageFilter<UShortImageType, UShortImageType>;
  auto medFilter = FilterType::New();

  // this is radius. 1 --> median window 3
  std::cout << "Post median(3D) filtering is under progress..Size(radius X Y "
               "Z) is = "
            << indexRadius << std::endl;

  medFilter->SetRadius(indexRadius);
  medFilter->SetInput(m_spCrntReconImg);
  medFilter->Update();

  m_spCrntReconImg = medFilter->GetOutput();
  std::cout << "median filtering has been done" << std::endl;
}

void CbctRecon::Export2DDoseMapAsMHA(fs::path &strPath) const {
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_spCrntReconImg == nullptr) {
    return;
  }

  if (strPath.empty()) {
    return;
  }

  const auto originLeft = static_cast<double>(m_spCrntReconImg->GetOrigin()[0]);
  const auto originTop =
      static_cast<double>(m_spCrntReconImg->GetOrigin()[1]); // not sure...

  const auto spacingX = static_cast<double>(m_spCrntReconImg->GetSpacing()[0]);
  const auto spacingY =
      static_cast<double>(m_spCrntReconImg->GetSpacing()[1]); // not sure...

  // Export float 2D image
  auto doseImg2D = FloatImage2DType::New();
  FloatImage2DType::SizeType doseSize{};
  doseSize[0] = m_dspYKReconImage->m_iWidth;
  doseSize[1] = m_dspYKReconImage->m_iHeight;

  FloatImage2DType::IndexType doseStart{};
  doseStart[0] = 0;
  doseStart[1] = 0;

  FloatImage2DType::RegionType doseRegion;
  doseRegion.SetSize(doseSize);
  doseRegion.SetIndex(doseStart);

  FloatImage2DType::SpacingType doseSpacing;
  doseSpacing[0] = spacingX;
  doseSpacing[1] = spacingY;

  FloatImage2DType::PointType doseOrigin;
  doseOrigin[0] = originLeft;
  doseOrigin[1] = originTop;

  doseImg2D->SetRegions(doseRegion);
  doseImg2D->SetSpacing(doseSpacing);
  doseImg2D->SetOrigin(doseOrigin);

  doseImg2D->Allocate();
  doseImg2D->FillBuffer(0);

  const auto factor_ushort2float = 0.01; // cGy --> Gy

  itk::ImageRegionIterator<FloatImage2DType> it(
      doseImg2D, doseImg2D->GetLargestPossibleRegion());

  auto pixel_val = 0.0f;
  size_t i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    pixel_val = static_cast<double>(m_dspYKReconImage->m_pData[i]) *
                factor_ushort2float;
    it.Set(pixel_val);
    ++i;
  }
  // YK201502
  using WriterType = itk::ImageFileWriter<FloatImage2DType>;
  auto writer = WriterType::New();
  writer->SetFileName(strPath.string());
  writer->SetUseCompression(true);
  writer->SetInput(doseImg2D);
  writer->Update();

  std::cout << "File was exported successfully" << std::endl;
}

void CbctRecon::ExportProjGeometryTXT(fs::path &strPath) const {
  // if (!m_spFullGeometry)
  //	return;

  if (m_spCustomGeometry ==
      nullptr) { // will be filled after Projection load button is pushed
    return;
  }

  if (strPath.empty()) {
    return;
  }

  const int cntAngle = m_spCustomGeometry->GetGantryAngles().size();
  const int cntShiftX = m_spCustomGeometry->GetProjectionOffsetsX().size();
  const int cntShiftY = m_spCustomGeometry->GetProjectionOffsetsY().size();

  if (cntAngle <= 0) {
    std::cout << "Error! no angle std::vector is found" << std::endl;
    return;
  }

  if (cntAngle != cntShiftX || cntAngle != cntShiftY) {
    std::cout << "Error! Angle number and shift number are not matching."
              << std::endl;
    return;
  }

  auto itShiftX = m_spCustomGeometry->GetProjectionOffsetsX().begin();
  auto itShiftY = m_spCustomGeometry->GetProjectionOffsetsY().begin();

  std::ofstream fout;
  fout.open(strPath);

  fout << "MV_Gantry_Angle"
       << "	"
       << "PanelShiftX(mm)"
       << "	"
       << "PanelShiftY(mm)" << std::endl;

  for (auto &itAng : m_spCustomGeometry->GetGantryAngles()) {
    fout << itAng << "	" << *itShiftX << "	" << *itShiftY << std::endl;

    ++itShiftX;
    ++itShiftY;
  }

  fout.close();
}

bool CbctRecon::LoadXVIGeometryFile(const fs::path &filePath) {

  m_spFullGeometry = GeometryType::New();

  itk::DOMNode::Pointer output_dom_object;
  itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
  reader->SetFileName(filePath.string());
  try {
    reader->Update();
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n"
              << "In XVI xml reader\n";
    return false;
  }
  output_dom_object = reader->GetOutput();

  m_vExcludeProjIdx.clear();
  auto iIdx = 0;

  auto dom = output_dom_object->Find(std::to_string(iIdx));
  while (dom != nullptr) {

    if (dom->GetName() == "Frame") {

      auto flxData = crl::XML_parseFrameForXVI5(dom);
      // m_vRetroFlexmap.push_back(flxData);

      if (flxData.fGanAngle < 0) {
        flxData.fGanAngle = flxData.fGanAngle + 360.0;
      }

      flxData.fPanelOffsetX = -flxData.fPanelOffsetX;
      flxData.fPanelOffsetY = -flxData.fPanelOffsetY;

      if (!flxData.bKV_On) {
        m_vExcludeProjIdx.push_back(iIdx);
      }

      /*               if (flxData.bKV_On)
                         m_vExcludeProjIdx.push_back(iIdx);*/

      ////Image qual test
      // flxData.fGanAngle = -flxData.fGanAngle;
      // if (flxData.fGanAngle < 0)
      //	flxData.fGanAngle = flxData.fGanAngle + 360.0;
      // flxData.fPanelOffsetX = -flxData.fPanelOffsetX;
      // flxData.fPanelOffsetY = -flxData.fPanelOffsetY;

      m_spFullGeometry->AddProjection(1000.0, 1536.0, flxData.fGanAngle,
                                      flxData.fPanelOffsetX,
                                      flxData.fPanelOffsetY, // Flexmap
                                      0.0, 0.0,  // In elekta, these are 0
                                      0.0, 0.0); // In elekta, these are 0
    }
    iIdx++;
    dom = output_dom_object->Find(std::to_string(iIdx));
  }
  return true;
}

void CbctRecon::GenerateCylinderMask(UShortImageType::Pointer &spImgCanvas,
                                     const float fDcmPosX, const float fDcmPosY,
                                     const float fRadius) const {
  if (spImgCanvas == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  auto origin = spImgCanvas->GetOrigin();
  auto spacing = spImgCanvas->GetSpacing();
  // UShortImageType::SizeType size =
  // spImgCanvas->GetBufferedRegion().GetSize();

  // itk::ImageSliceConstIteratorWithIndex<FloatImageType> it (m_spReconImg,
  // m_spReconImg->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      spImgCanvas, spImgCanvas->GetBufferedRegion());

  // ImageSliceConstIteratorWithIndex<ImageType> it( image,
  // image->GetRequestedRegion() );  UShortImageType::SizeType imgSize =
  // spImgCanvas->GetRequestedRegion().GetSize(); //1016x1016 x z

  // int width = imgSize[0];
  // int height = imgSize[1];

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  auto iNumSlice = 0;

  // int i = 0;//height
  // int j = 0; // width

  while (!it.IsAtEnd()) {
    auto iPosY = 0;
    while (!it.IsAtEndOfSlice()) {
      auto iPosX = 0;
      while (!it.IsAtEndOfLine()) {
        // Calculate physical position

        const auto crntPhysX = iPosX * static_cast<double>(spacing[0]) +
                               static_cast<double>(origin[0]);
        const auto crntPhysY = iPosY * static_cast<double>(spacing[1]) +
                               static_cast<double>(origin[1]);

        if (pow(crntPhysX - fDcmPosX, 2.0) + pow(crntPhysY - fDcmPosY, 2.0) >=
            pow(fRadius, 2.0)) {
          //(*it) = (unsigned short)0; //air value
          it.Set(0);
        } else {
          it.Set(1);
        }

        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }
}

float CbctRecon::GetMeanIntensity(UShortImageType::Pointer &spImg,
                                  const float sphereR,
                                  float *sdIntensity) const {
  if (spImg == nullptr) {
    return -1.0;
  }

  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  auto origin = spImg->GetOrigin();
  auto spacing = spImg->GetSpacing();
  // UShortImageType::SizeType size = spImg->GetBufferedRegion().GetSize();

  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      spImg, spImg->GetBufferedRegion());
  // UShortImageType::SizeType imgSize =
  // spImg->GetRequestedRegion().GetSize();
  // //1016x1016 x z

  // int width = imgSize[0];
  // int height = imgSize[1];

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  auto iNumSlice = 0;

  auto pixSum = 0.0;
  auto iCnt = 0;

  while (!it.IsAtEnd()) {
    auto iPosY = 0;
    while (!it.IsAtEndOfSlice()) {
      auto iPosX = 0;
      while (!it.IsAtEndOfLine()) {
        // Calculate physical position

        const auto crntPhysX = iPosX * static_cast<double>(spacing[0]) +
                               static_cast<double>(origin[0]);
        const auto crntPhysY = iPosY * static_cast<double>(spacing[1]) +
                               static_cast<double>(origin[1]);
        const auto crntPhysZ = iNumSlice * static_cast<double>(spacing[2]) +
                               static_cast<double>(origin[2]);

        if (pow(crntPhysX, 2.0) + pow(crntPhysY, 2.0) + pow(crntPhysZ, 2.0) <
            pow(sphereR, 2.0)) {
          pixSum = pixSum + static_cast<double>(it.Get());
          iCnt++;
        }
        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }

  float meanIntensity;
  if (iCnt > 0) {
    meanIntensity = pixSum / static_cast<double>(iCnt);
  } else {
    meanIntensity = -1.0;
  }

  if (sdIntensity == nullptr) {
    return meanIntensity;
  }

  auto devSum = 0.0;
  it.GoToBegin();

  iNumSlice = 0;

  while (!it.IsAtEnd()) {
    auto iPosY = 0;
    while (!it.IsAtEndOfSlice()) {
      auto iPosX = 0;
      while (!it.IsAtEndOfLine()) {
        // Calculate physical position

        const auto crntPhysX = iPosX * static_cast<double>(spacing[0]) +
                               static_cast<double>(origin[0]);
        const auto crntPhysY = iPosY * static_cast<double>(spacing[1]) +
                               static_cast<double>(origin[1]);
        const auto crntPhysZ = iNumSlice * static_cast<double>(spacing[2]) +
                               static_cast<double>(origin[2]);

        if (pow(crntPhysX, 2.0) + pow(crntPhysY, 2.0) + pow(crntPhysZ, 2.0) <
            pow(sphereR, 2.0)) {
          devSum =
              devSum + pow(static_cast<double>(it.Get()) - meanIntensity, 2.0);
        }
        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }

  if (iCnt > 0) {
    *sdIntensity = sqrt(devSum / static_cast<double>(iCnt));
  } else {
    *sdIntensity = -1.0;
  }

  return meanIntensity;
}

bool CbctRecon::ResortCBCTProjection(
    std::vector<int> &vIntPhaseBinSelected, fs::path &strPathForXML,
    fs::path &strPathProjRoot, std::string &strUID,
    std::vector<float> &vFloatPhaseFull, GeometryType::Pointer &spGeomFull,
    std::vector<std::string> &vProjPathsFull) const {
  if (vIntPhaseBinSelected.empty()) {
    return false;
  }

  const auto NumOfPhaseFull = vFloatPhaseFull.size();
  const auto NumOfGeomFull = spGeomFull->GetGantryAngles().size();
  const auto NumOfProjFileFull = vProjPathsFull.size();

  if (NumOfPhaseFull != NumOfGeomFull || NumOfGeomFull != NumOfProjFileFull) {
    std::cout << "Num of data is not matching:"
              << " NumOfPhaseFull= " << NumOfPhaseFull
              << " NumOfGeomFull= " << NumOfGeomFull
              << " NumOfProjFileFull= " << NumOfProjFileFull << std::endl;
    return false;
  }

  // Check Root dir is set

  if (strUID.length() < 1) {
    return false;
  }

  auto dirSaveXML = strPathForXML;
  auto dirSaveProj = strPathProjRoot;

  if (!fs::exists(dirSaveXML) || !fs::exists(dirSaveProj)) {
    std::cerr << "Error! Directories don't exist\n";
    ;
    return false;
  }

  auto strUID_Endfix = std::accumulate(
      vIntPhaseBinSelected.begin(), vIntPhaseBinSelected.end(), "P"s,
      [](std::string acc_str, auto phase) {
        std::array<char, 2> strNum{{'0', '0'}};
        // [ptr, ec] =
        std::to_chars(strNum.data(), strNum.data() + strNum.size(), phase, 10);
        // strNum = QString("%1").arg(vIntPhaseBinSelected.at(i), 2, 10, zero);
        return std::move(acc_str) + std::string(strNum.data(), strNum.size());
      });
  strUID_Endfix = strUID_Endfix + "P"; // UID...P00102030405060P
  const auto strNewUID = strUID + strUID_Endfix;

  // Create a subDir
  const auto strSubDirName = "img_" + strNewUID;

  const auto projDir = strPathProjRoot / strSubDirName;
  if (!fs::create_directory(projDir)) {
    std::cerr << "Could not make subdir" << std::endl;
    return false;
  }

  if (!fs::exists(strPathForXML)) {
    std::cout << "no XML Dir exists" << std::endl;
    return false;
  }

  // strPathProj
  // strPathForXML
  std::vector<size_t> vSelectedIdxTemp;
  std::vector<size_t> vSelectedIdxFin;

  for (auto& phase : vIntPhaseBinSelected) {
    AppendInPhaseIndex(phase, vFloatPhaseFull,
                       vSelectedIdxTemp);
  }
  // Remove redandancy

  sort(vSelectedIdxTemp.begin(), vSelectedIdxTemp.end()); // hopefully,
                                                          // ascending
  std::cout << "sorting check" << std::endl;
  std::cout << "0 " << vSelectedIdxTemp.at(0) << std::endl;
  std::cout << "1 " << vSelectedIdxTemp.at(1) << std::endl;

  auto prevVal = -1;
  for (auto &it : vSelectedIdxTemp) {
    if (static_cast<long long signed>(it) > prevVal) {
      vSelectedIdxFin.push_back(it);
    }

    prevVal = static_cast<int>(it);
  }

  auto spSubGeometry = GeometryType::New();

  for (auto &itIdx : vSelectedIdxFin) {
    std::cout << "cur Idx=" << itIdx << std::endl;
    // 9 parameters are required
    const auto curSID = spGeomFull->GetSourceToIsocenterDistances().at(itIdx);
    const auto curSDD = spGeomFull->GetSourceToDetectorDistances().at(itIdx);
    const auto curGantryAngle = spGeomFull->GetGantryAngles().at(itIdx);

    const auto curProjOffsetX = spGeomFull->GetProjectionOffsetsX().at(itIdx);
    const auto curProjOffsetY = spGeomFull->GetProjectionOffsetsY().at(itIdx);

    const auto curOutOfPlaneAngles =
        spGeomFull->GetOutOfPlaneAngles().at(itIdx);
    const auto curInPlaneAngles = spGeomFull->GetInPlaneAngles().at(itIdx);

    const auto curSrcOffsetX = spGeomFull->GetSourceOffsetsX().at(itIdx);
    const auto curSrcOffsetY = spGeomFull->GetSourceOffsetsY().at(itIdx);

    spSubGeometry->AddProjection(
        curSID, curSDD, curGantryAngle, curProjOffsetX,
        curProjOffsetY,                        // Flexmap
        curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
        curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0
  }
  // Export spSubGeometry
  auto xmlWriter = rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();

  const auto geomFileName = "ElektaGeom_" + strNewUID + ".xml";
  auto geomFilePath = strPathForXML / geomFileName;

  xmlWriter->SetFilename(geomFilePath.string());
  xmlWriter->SetObject(spSubGeometry);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
  // Copy selected his files to a different folder

  for (auto &itIdx : vSelectedIdxFin) {
    std::string strPathProjOriginal = vProjPathsFull.at(itIdx).c_str();
    // Copy this file to target dir

    fs::path fInfo(strPathProjOriginal);
    auto strPathProjNew = projDir / fInfo.filename();
    fs::copy(fs::absolute(fInfo), strPathProjNew);
  }

  //    std::vector<float>& vFloatPhaseFull, GeometryType::Pointer&
  //    spGeomFull, std::vector<std::string>& vProjPathsFull
  std::cout << vSelectedIdxFin.size() << " files were copied." << std::endl;

  return true;
}

void CbctRecon::AppendInPhaseIndex(const int iPhase,
                                   const std::vector<float> &vFloatPhaseFull,
                                   std::vector<size_t> &vOutputIndex,
                                   const int margin) const {

  int startPhase1;
  int endPhase1;

  int startPhase2;
  int endPhase2;

  size_t i = 0;
  for (const auto &phase : vFloatPhaseFull) {
    const auto iCurPhase = qRound(phase * 100.0);
    // determine wether it is within the range

    if (iPhase < margin) // if 5 --> 0 ~ 10%, IF 4--> 99 ~ 09
    {
      startPhase2 = iPhase + 100 - margin;
      endPhase2 = 100;

      startPhase1 = 0;
      endPhase1 = iPhase + margin;
    } else {
      startPhase1 = iPhase - margin;
      endPhase1 = iPhase + margin;

      startPhase2 = 1; // reverse
      endPhase2 = 0;
    }

    if ((iCurPhase >= startPhase1 && iCurPhase <= endPhase1) ||
        (iCurPhase >= startPhase2 && iCurPhase <= endPhase2)) {
      vOutputIndex.push_back(i);
    }
    ++i;
  }
}

void CbctRecon::LoadShort3DImage(fs::path &filePath,
                                 const enREGI_IMAGES enTarget) {
  if (!fs::exists(filePath)) {
    return;
  }

  UShortImageType::Pointer spImg;

  if (!crl::LoadShortImageToUshort(filePath, spImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  switch (enTarget) {
  case enREGI_IMAGES::REGISTER_RAW_CBCT:
    m_spRawReconImg = spImg;
    break;
  case enREGI_IMAGES::REGISTER_COR_CBCT:
    m_spScatCorrReconImg = spImg;
    break;
  case enREGI_IMAGES::REGISTER_MANUAL_RIGID:
    m_spManualRigidCT = spImg;
    break;
  case enREGI_IMAGES::REGISTER_AUTO_RIGID:
    m_spAutoRigidCT = spImg;
    break;
  case enREGI_IMAGES::REGISTER_DEFORM_FINAL:
    m_spDeformedCT_Final = spImg;
    break;
  default:
    m_spRawReconImg = spImg;
    break;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  auto imageCalculatorFilter2 = ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(spImg);
  imageCalculatorFilter2->Compute();

  const auto minVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMinimum());
  const auto maxVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMaximum());

  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  auto imgDim = spImg->GetBufferedRegion().GetSize();
  auto spacing = spImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  m_spCrntReconImg = spImg.GetPointer();

  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);
}

void CbctRecon::GetWEPLDataFromSingleFile(const fs::path &filePath,
                                          std::vector<VEC3D> &vPOI,
                                          std::vector<WEPLData> &vOutputWEPL,
                                          const double fAngleStart,
                                          const double fAngleEnd) const {

  const auto iCntPOI = vPOI.size();

  if (iCntPOI < 1) {
    return;
  }

  const auto fAngleGap = 1.0f;

  UShortImageType::Pointer spImg;

  auto strFilePath = filePath;
  if (!crl::LoadShortImageToUshort(strFilePath, spImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
    return;
  }

  for (size_t i = 0; i < iCntPOI; i++) {
    const auto curPOI = vPOI.at(i);
    // append mode
    GetAngularWEPL_SinglePoint(spImg, fAngleGap, fAngleStart, fAngleEnd, curPOI,
                               i, vOutputWEPL, true); // mandatory
  }
}

void CbctRecon::ScatterCorPerProjRef(const double scaMedian,
                                     const double scaGaussian,
                                     const int postScatMedianSize,
                                     const bool use_cuda, const bool use_opencl,
                                     const bool save_dicom,
                                     FDK_options &fdk_options) // load text file
{
  if (m_strListPerProjRefVol.empty()) {
    std::cout << "Error! Ref Vol list is not ready yet. Load it first"
              << std::endl;
    return;
  }

  // Find the mask file
  // std::string strPath_mskSkinCT_final;
  // std::string strPath_mskSkinCT_autoRegi_exp =
  // m_pDlgRegistration->m_strPathPlastimatch +
  // "/msk_skin_CT_autoRegi_exp.mha"; QFileInfo
  // maskInfoAuto(strPath_mskSkinCT_autoRegi_exp);

  // std::string strPath_mskSkinCT_manualRegi_exp =
  // m_pDlgRegistration->m_strPathPlastimatch +
  // "/msk_skin_CT_manRegi_exp.mha"; QFileInfo
  // maskInfoManual(strPath_mskSkinCT_manualRegi_exp);

  // if (maskInfoAuto.exists()) //if the mask file is not prepared, give up
  // the skin removal
  //{
  //    strPath_mskSkinCT_final = strPath_mskSkinCT_autoRegi_exp;
  //}
  // else
  //{
  //    std::cout << "Mask file of auto-registration is not prepared. Use
  //    manual regi-mask instead" << std::endl;

  //    if (maskInfoManual.exists())
  //    {
  //        strPath_mskSkinCT_final = strPath_mskSkinCT_manualRegi_exp;
  //    }
  //    else
  //    {
  //        std::cout << "Mask file of manual registration is not prepared.
  //        Skip skin removal!" << std::endl; return;
  //    }
  //}

  ////std::cout << "Plastimatch Path " <<
  /// m_strPathPlastimatch.toLocal8Bit().constData() << std::endl;

  // if (m_pDlgRegistration->m_strPathPlastimatch.length() < 1)
  //{
  //    std::cout << "NO plastimatch Dir was defined. CorrCBCT will not be
  //    saved automatically" << std::endl; return;
  //}
  // Forward proj

  // Make a canvas

  if (m_spProjImgRaw3D == nullptr) {
    std::cout << "ERRORRR! m_spProjImgRaw3D" << std::endl;
    return;
  }

  auto spProjImgCT3D = FloatImageType::New(); // later
  const auto projCT_size =
      m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize(); // 1024 1024 350
  const auto projCT_idxStart =
      m_spProjImgRaw3D->GetLargestPossibleRegion().GetIndex(); // 0 0 0
  const auto projCT_spacing = m_spProjImgRaw3D->GetSpacing();  // 0.4 0.4 1.0
  const auto projCT_origin =
      m_spProjImgRaw3D->GetOrigin(); //-204.6 -204.6 -174.5

  FloatImageType::RegionType projCT_region;
  projCT_region.SetSize(projCT_size);
  projCT_region.SetIndex(projCT_idxStart);

  spProjImgCT3D->SetRegions(projCT_region);
  spProjImgCT3D->SetSpacing(projCT_spacing);
  spProjImgCT3D->SetOrigin(projCT_origin);

  spProjImgCT3D->Allocate();
  spProjImgCT3D->FillBuffer(0);

  // YKTEMP
  const auto proj_size = spProjImgCT3D->GetBufferedRegion().GetSize();
  std::cout << "ProjImgCT Size = " << proj_size[0] << ", " << proj_size[1]
            << ", " << proj_size[2] << std::endl;
  std::cout << "ProjImgCT origin = " << spProjImgCT3D->GetOrigin()[0] << ", "
            << spProjImgCT3D->GetOrigin()[1] << ", "
            << spProjImgCT3D->GetOrigin()[2] << std::endl;
  std::cout << "ProjImgCT spacing = " << spProjImgCT3D->GetSpacing()[0] << ", "
            << spProjImgCT3D->GetSpacing()[1] << ", "
            << spProjImgCT3D->GetSpacing()[2] << std::endl;

  const auto iCntRefVol = m_strListPerProjRefVol.size();

  if (iCntRefVol < 1) {
    std::cout << "Error! no volume data for loading" << std::endl;
    return;
  }

  const auto flexCnt =
      m_spCustomGeometry->GetGantryAngles().size();
  if (flexCnt != iCntRefVol) {
    std::cout << "Error! flex count doesn't match" << std::endl;
    return;
  }

  for (auto i = 0ull; i < iCntRefVol; i++) {
    // Load volume: Short image
    auto spOutputShort_raw = ShortImageType::New();
    // ShortImageType::Pointer spOutputShort_threshold =
    // ShortImageType::New();
    auto spOutputUshort = UShortImageType::New();
    // UShortImageType::Pointer spOutputUshort_register =
    // UShortImageType::New();
    auto spUshortRotated = UShortImageType::New();
    auto spAttFloat = FloatImageType::New();

    fs::path strDirPath{m_strListPerProjRefVol.at(i)};

    if (!crl::LoadShortImageDirOrFile(strDirPath, spOutputShort_raw)) {
      std::cout << "Error! in " << i
                << " th image. File couldn't be found. Path= "
                << strDirPath << std::endl;
      return;
    }

    crl::ConvertShort2Ushort(spOutputShort_raw, spOutputUshort);

    crl::RotateImgBeforeFwd(spOutputUshort,
                       spUshortRotated); // IEC to RTK w/ kVGantry
    crl::ConvertUshort2AttFloat(spUshortRotated, spAttFloat);

    const auto curMVAngle = m_spCustomGeometry->GetGantryAngles().at(i);
    const auto curPanelOffsetX =
        m_spCustomGeometry->GetProjectionOffsetsX().at(i);
    const auto curPanelOffsetY =
        m_spCustomGeometry->GetProjectionOffsetsY().at(i);

#if USE_CUDA
    if (use_cuda) {
      SingleForwardProjection<CUDAFloatImageType>(
          spAttFloat, curMVAngle, curPanelOffsetX, curPanelOffsetY,
          spProjImgCT3D, i);
    } else
#endif
    {
      SingleForwardProjection<FloatImageType>(spAttFloat, curMVAngle,
                                              curPanelOffsetX, curPanelOffsetY,
                                              spProjImgCT3D, i);
    }
    std::cout << "Proj: " << i << "/" << iCntRefVol << std::endl;
  }

  std::cout << "Generating scatter map is ongoing..." << std::endl;

  std::cout << "To account for the mAs values, the intensity scale factor of "
            << crl::GetRawIntensityScaleFactor(m_strRef_mAs, m_strCur_mAs)
            << "will be multiplied during scatter correction to avoid negative "
               "scatter"
            << std::endl;
  auto cast_ushort_to_float =
      itk::CastImageFilter<UShortImageType, FloatImageType>::New();
  cast_ushort_to_float->SetInput(m_spProjImgRaw3D);
  cast_ushort_to_float->Update();
  FloatImageType::Pointer spProjImg3DFloat = cast_ushort_to_float->GetOutput();

  GenScatterMap_PriorCT(spProjImg3DFloat, spProjImgCT3D, m_spProjImgScat3D,
                        scaMedian, scaGaussian,
                        false); // void GenScatterMap2D_PriorCT()
  spProjImgCT3D->Initialize();  // memory saving

  std::cout << "Scatter correction is in progress..." << std::endl;

  ScatterCorr_PrioriCT(spProjImg3DFloat, m_spProjImgScat3D, m_spProjImgCorr3D,
                       postScatMedianSize, true);
  m_spProjImgScat3D->Initialize(); // memory saving

  auto cast_float_to_ushort =
      itk::CastImageFilter<FloatImageType, UShortImageType>::New();
  cast_float_to_ushort->SetInput(spProjImg3DFloat);
  cast_float_to_ushort->Update();
  m_spProjImgRaw3D = cast_float_to_ushort->GetOutput();

  std::cout << "AfterCorrectionMacro is ongoing..." << std::endl;
  AfterScatCorrectionMacro(use_cuda, use_opencl, save_dicom, fdk_options);
  std::cout << "FINISHED!Scatter correction: CBCT DICOM files are saved"
            << std::endl;
}

// double CbctRecon::CropSkinUsingRS(UShortImageType::Pointer& spImgUshort,
// std::string& strPathRS, double cropMargin )
//{
//   //if (m_pParent->m_strPathRS.isEmpty())
//	//return;
//  /* End of [1]Segment air region*/
//
//  //plastimatch convert --input E:\PlastimatchData\DicomEg\OLD\RS.dcm
//  --output-ss-img E:\PlastimatchData\DicomEg\OLD\ssimg_all2.mha
//  --output-ss-list E:\PlastimatchData\DicomEg\OLD\sslist_all.txt
//  --referenced-ct E:\PlastimatchData\DicomEg\OLD\CT
//
//  //plm_clp_parse (&parms, &parse_fn, &usage_fn, argc, argv, 1)
//  //plastimatch segment --input E:\PlastimatchData\DicomEg\OLD\CT
//  --output-img E:\PlastimatchData\DicomEg\OLD\msk_bubbles_oldCT.mha
//  --lower-threshold -600
//
//  //do_command_warp(argc, argv);
//
//
//
//
//}
