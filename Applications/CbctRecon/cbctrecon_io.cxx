/*All IO functions used with cbctrecon*/
#include "cbctrecon.h"

// std
#include <iostream> // for operator<<, basic_o...
#include <string>   // for string
#include <vector>   // for vector

// ITK
#include "gdcmUIDGenerator.h"
#include "itkGDCMImageIO.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkThresholdImageFilter.h"

// RTK
#include "rtkElektaSynergyGeometryReader.h"
#include "rtkMacro.h" // for TRY_AND_EXIT_ON_ITK...
#include "rtkProjectionsReader.h"
#include "rtkThreeDCircularProjectionGeometryXMLFileReader.h" // for ThreeDCircularProje...
#include "rtkThreeDCircularProjectionGeometryXMLFileWriter.h" // for ThreeDCircularProje...

// PLM
#include "dcmtk_rt_study.h"
#include "plm_image.h"
#include <rt_study_metadata.h>

// local
#include "StructureSet.h"
#include "cbctrecon_io.h"

QString MakeElektaXML(const QString &filePath_ImageDBF,
                      const QString &filePath_FrameDBF,
                      const QString &DICOM_UID) {
  std::cout << "Elekta geometry XML file is being generated." << std::endl;
  // Define FRAME.DBF path
  rtk::ElektaSynergyGeometryReader::Pointer reader =
      rtk::ElektaSynergyGeometryReader::New();
  // string strDicomUID = DICOM_UID.toLocal8Bit().constData();
  // //DICOM_UID.toStdString()  string strDicomUID = DICOM_UID.toStdString();
  // string strDbfImg = filePath_ImageDBF.toStdString();
  // string strDbfFrame = filePath_FrameDBF.toStdString();

  QFileInfo info = QFileInfo(filePath_ImageDBF);
  const QString dirPath = info.absolutePath();

  const QString fileName = "ElektaGeom_" + DICOM_UID + ".xml";

  QString str_output = dirPath + "/" + fileName;

  reader->SetDicomUID(DICOM_UID.toLocal8Bit().constData());
  reader->SetImageDbfFileName(filePath_ImageDBF.toLocal8Bit().constData());
  reader->SetFrameDbfFileName(filePath_FrameDBF.toLocal8Bit().constData());

  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->UpdateOutputData());

  // Write
  rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter =
      rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
  xmlWriter->SetFilename(str_output.toLocal8Bit().constData());
  xmlWriter->SetObject(reader->GetGeometry());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile())

  std::cout << "Reading succeed" << std::endl;

  return str_output;
}

// Get the projection geometry
void CbctRecon::LoadRTKGeometryFile(const char *filePath) {
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader =
      rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(filePath);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())
  std::cout << "Geometry reading succeed" << std::endl;

  m_spFullGeometry = geometryReader->GetOutputObject();

  // fullGeometry->GetGantryAngles();
  const int geoDataSize =
      m_spFullGeometry->GetGantryAngles().size(); // This is MV gantry angle!!!
  std::cout << "Geometry data size(projection gantry angles): " << geoDataSize
            << std::endl;
  if (geoDataSize < 1) {
    return;
  }

  // CW: continuously ascending except 360 - 0 interface, no negative value
  // CCW: continuously descending except 0 - 360 interface, no negative value

  // m_bScanDirectionCW

  /*for (int i = 0 ; i<geoDataSize ; i++)
  {
  std::cout << "Projection gantry angle " << i << ": " <<
  m_spFullGeometry->GetGantryAngles().at(i) << std::endl;
  }
  std::cout << "Coordination is following IEC coordinate (e.g. 190 deg = RPO,
  170 LPO, 30 = LAO supposing HFS) no and kV source position (not MV gantry)" <<
  std::endl;*/

  std::vector<double> vTempConvAngles;

  const auto itBegin = (m_spFullGeometry->GetGantryAngles()).begin();
  const auto itEnd = (m_spFullGeometry->GetGantryAngles()).end();

  for (auto it = itBegin; it != itEnd; ++it) {
    double tmpAngle = (*it);

    if (tmpAngle > 180.0) {
      tmpAngle = tmpAngle - 360.0;
    }

    vTempConvAngles.push_back(tmpAngle);
  }

  // compare 2 points in the middle of the angle list
  const auto iLowerIdx = static_cast<int>(geoDataSize * 1.0 / 3.0);
  const auto iUpperIdx = static_cast<int>(geoDataSize * 2.0 / 3.0);

  if (vTempConvAngles.at(iLowerIdx) <
      vTempConvAngles.at(iUpperIdx)) // ascending
  {
    m_bScanDirectionCW = true;
    std::cout << "The scan direction is CW" << std::endl;
  }

  else {
    m_bScanDirectionCW = false;
    std::cout << "The scan direction is CCW" << std::endl;
  }

  std::cout << "AngularGaps Size: "
            << m_spFullGeometry
                   ->GetAngularGaps(m_spFullGeometry->GetSourceAngles())
                   .size()
            << std::endl;
}

bool LoadShortImageToUshort(QString &strPath,
                            UShortImageType::Pointer &pUshortImage) {
  using ReaderType = itk::ImageFileReader<ShortImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  // QString fileName = QFileDialog::getOpenFileName(this, "Open Image","",
  // "Plan CT file (*.mha)",0,0);

  if (strPath.length() < 1) {
    return false;
  }

  reader->SetFileName(strPath.toLocal8Bit().constData());
  reader->Update();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(reader->GetOutput());
  imageCalculatorFilter->Compute();

  const auto minVal0 = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal0 = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Original Min and Max Values are	" << minVal0 << "	"
            << maxVal0 << std::endl;

  bool bNKI = false;
  if (minVal0 > -600) // impossible for normal Short image. IN NKI, always -512.
                      // don't know why
  {
    bNKI = true;
  }

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilter =
      ThresholdImageFilterType::New();

  if (!bNKI) {
    thresholdFilter->SetInput(reader->GetOutput());
    thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
    thresholdFilter->SetOutsideValue(-1024);
    thresholdFilter->Update();
  } else {
    thresholdFilter->SetInput(reader->GetOutput());
    thresholdFilter->ThresholdOutside(0, 4095); //--> 0 ~ 4095
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();
  }

  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();

  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Current Min and Max Values are	" << minVal << "	"
            << maxVal << std::endl;

  USHORT_PixelType outputMinVal, outputMaxVal;
  if (!bNKI) {
    outputMinVal = static_cast<USHORT_PixelType>(minVal + 1024);
    outputMaxVal = static_cast<USHORT_PixelType>(maxVal + 1024);
  } else {
    outputMinVal = static_cast<USHORT_PixelType>(minVal);
    outputMaxVal = static_cast<USHORT_PixelType>(maxVal);
  }

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();
  pUshortImage = spRescaleFilter->GetOutput();

  return true;
}

void ExportReconSHORT_HU(UShortImageType::Pointer &spUsImage,
                         QString &outputFilePath) {
  if (spUsImage == nullptr) {
    std::cout << " no image to export" << std::endl;
    return;
  }

  using DuplicatorType = itk::ImageDuplicator<UShortImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(spUsImage);
  duplicator->Update();
  const UShortImageType::Pointer clonedReconImage = duplicator->GetOutput();
  ShortImageType::Pointer clonedReconImageSHORT;

  using ThresholdImageFilterType = itk::ThresholdImageFilter<UShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilterAbove =
      ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(clonedReconImage);
  thresholdFilterAbove->ThresholdAbove(4095);
  thresholdFilterAbove->SetOutsideValue(4095);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
      ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(0);
  thresholdFilterBelow->SetOutsideValue(0);
  thresholdFilterBelow->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<UShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilterBelow->GetOutput());
  imageCalculatorFilter->Compute();
  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  const auto outputMinVal = static_cast<SHORT_PixelType>(minVal - 1024);
  const auto outputMaxVal = static_cast<SHORT_PixelType>(maxVal - 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilterBelow->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);

  spRescaleFilter->Update();
  // clonedReconImageSHORT = spRescaleFilter->GetOutput();

  // waterHU = 1024;
  /*
  using AddImageFilterType =
      itk::AddImageFilter<ShortImageType, ShortImageType, ShortImageType>;
  AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput1(spRescaleFilter->GetOutput());

  const int addingVal = 0; // 1024-680
  addImageFilter->SetConstant2(addingVal);
  addImageFilter->Update();
  */
  using WriterType = itk::ImageFileWriter<ShortImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFilePath.toLocal8Bit().constData());
  // writer->SetUseCompression(true);
  writer->SetUseCompression(false); // for plastimatch
  writer->SetInput(spRescaleFilter->GetOutput());

  std::cout << "Writing is under progress...: "
            << outputFilePath.toLocal8Bit().constData() << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;
}

bool CbctRecon::ReadDicomDir(QString &dirPath) {
  Dcmtk_rt_study drs(dirPath.toLocal8Bit().constData());
  drs.load_directory(); // parse_directory();

  Plm_image plmImg;
  auto tmp_img = drs.get_image();

  if (!tmp_img->have_image()) {
    return false;
  }
  std::cout << "PLM_imagetype: " << tmp_img->m_type << std::endl;
  plmImg.set(tmp_img);
  // plmImg.load_native(dirPath.toLocal8Bit().constData());

  auto planCT_ss = drs.get_rtss(); // dies at end of scope...
  if (!planCT_ss) {
    // ... so I copy to my own modern-C++ implementation
    m_structures->set_planCT_ss(planCT_ss.get());
  }

  const ShortImageType::Pointer spShortImg = plmImg.itk_short();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilter =
      ThresholdImageFilterType::New();

  thresholdFilter->SetInput(spShortImg);
  thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
  thresholdFilter->SetOutsideValue(-1024);
  thresholdFilter->Update();

  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();

  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Current Min and Max Values are	" << minVal << "	"
            << maxVal << std::endl;

  const auto outputMinVal = static_cast<USHORT_PixelType>(minVal + 1024);
  const auto outputMaxVal = static_cast<USHORT_PixelType>(maxVal + 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  // m_spRawReconImg = spRescaleFilter->GetOutput();
  m_spRefCTImg = spRescaleFilter->GetOutput();
  return true;
}

// From DlgRegistration

void ConvertUshort2Short(UShortImageType::Pointer &spImgUshort,
                         ShortImageType::Pointer &spImgShort) {
  using ThresholdImageFilterType = itk::ThresholdImageFilter<UShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilter =
      ThresholdImageFilterType::New();
  thresholdFilter->SetInput(spImgUshort);
  thresholdFilter->ThresholdOutside(0, 4096); //--> 0 ~ 4095
  thresholdFilter->SetOutsideValue(0);
  thresholdFilter->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<UShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();
  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  const auto outputMinVal = static_cast<SHORT_PixelType>(minVal - 1024);
  const auto outputMaxVal = static_cast<SHORT_PixelType>(maxVal - 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  spImgShort = spRescaleFilter->GetOutput();
}


void ConvertShort2Ushort(ShortImageType::Pointer &spInputImgShort,
                         UShortImageType::Pointer &spOutputImgUshort) {
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilterAbove =
      ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(spInputImgShort);
  thresholdFilterAbove->ThresholdAbove(3071);
  thresholdFilterAbove->SetOutsideValue(3071);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
      ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(-1024);
  thresholdFilterBelow->SetOutsideValue(-1024);
  thresholdFilterBelow->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilterBelow->GetOutput());
  imageCalculatorFilter->Compute();
  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  const auto outputMinVal =
      static_cast<UShortImageType::PixelType>(minVal + 1024);
  const auto outputMaxVal =
      static_cast<UShortImageType::PixelType>(maxVal + 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilterBelow->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  spOutputImgUshort = spRescaleFilter->GetOutput();
}

QString SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
                                QString &strPatientID, QString &strPatientName,
                                QString &strPathTargetDir) {
  if (spImg == nullptr) {
    return QString("");
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  Plm_image plm_img(spShortImg);

  QString endFix = "_DCM";

  QString newDirPath = strPathTargetDir + "/" + strPatientID + "_DCM";
  // QString newDirPath =
  //  strPathTargetDir + "/" + strPatientID + strPatientName + "_DCM";

  QDir dirNew(newDirPath);
  if (!dirNew.exists()) {
    if (!dirNew.mkdir(".")) {
      std::cerr << "Could not create dir for dicom" << std::endl;
      return QString("");
    }
  }

  Rt_study_metadata rsm;
  rsm.set_patient_id(strPatientID.toLocal8Bit().constData());
  rsm.set_patient_name(strPatientName.toLocal8Bit().constData());

  plm_img.save_short_dicom(newDirPath.toLocal8Bit().constData(), &rsm);

  return newDirPath.toLocal8Bit().constData();
}

QString SaveUSHORTAsSHORT_DICOM_gdcmITK(UShortImageType::Pointer &spImg,
                                        QString &strPatientID,
                                        QString &strPatientName,
                                        QString &strPathTargetDir) {
  if (spImg == nullptr) {
    return "";
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  QString newDirPath =
      strPathTargetDir + "/" + strPatientID + strPatientName + "_DCM";

  QDir dirNew(newDirPath);
  if (!dirNew.exists()) {
    if (!dirNew.mkdir(".")) {
      std::cerr << "Could not create dir for gdcm dicom" << std::endl;
      return QString("");
    }
  } else {
    if (dirNew.removeRecursively()) {
      QDir dirReNew(newDirPath);
      if (!dirNew.mkdir(".")) {
        std::cerr << "Could not create dir for gdcm dicom (2)" << std::endl;
        return QString("");
      }
    }
  }
  using OutputImageType =
      itk::Image<USHORT_PixelType,
                 2>; // because dicom is one 2d image for each slice-file
  using ImageIOType = itk::GDCMImageIO;
  using NamesGeneratorType = itk::NumericSeriesFileNames;

  UShortImageType::RegionType region = spShortImg->GetLargestPossibleRegion();
  UShortImageType::IndexType start = region.GetIndex();
  UShortImageType::SizeType size = region.GetSize();

  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  itk::MetaDataDictionary &dict = gdcmIO->GetMetaDataDictionary();
  std::string value = "CT";
  itk::EncapsulateMetaData<std::string>(dict, "0008|0060", value); // Modality
  value = "DERIVED\\SECONDARY\\AXIAL"; // This is virtually always correct when
                                       // using ITK to write an image
  itk::EncapsulateMetaData<std::string>(dict, "0008|0008", value); // Image Type
  value = "SI";
  itk::EncapsulateMetaData<std::string>(dict, "0008|0064",
                                        value); // Conversion Type
  const double value_double = spShortImg->GetSpacing()[2];
  std::ostringstream strs;
  strs << value_double;
  value = strs.str();
  std::cout << "slice spacing: " + value << std::endl;
  itk::EncapsulateMetaData<std::string>(dict, "0018|0050",
                                        value); // SliceThickness
  itk::EncapsulateMetaData<std::string>(dict, "0018|0088",
                                        '-' + value); // SpacingBetweenSlices

  gdcm::UIDGenerator stduid;
  const std::string studyUID = stduid.Generate();
  std::cout << studyUID << std::endl;
  itk::EncapsulateMetaData<std::string>(dict, "0020|000d", studyUID);

  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
  namesGenerator->SetStartIndex(static_cast<itk::SizeValueType>(start[2]));
  namesGenerator->SetEndIndex(static_cast<itk::SizeValueType>(start[2]) +
                              size[2] - 1);
  namesGenerator->SetIncrementIndex(1);
  namesGenerator->SetSeriesFormat(newDirPath.toStdString() + "/CT." + studyUID +
                                  ".%d.dcm");

  using SeriesWriterType =
      itk::ImageSeriesWriter<ShortImageType, OutputImageType>;
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput(spShortImg);
  seriesWriter->SetImageIO(gdcmIO);
  seriesWriter->SetFileNames(namesGenerator->GetFileNames());

  try {
    seriesWriter->Update();
  } catch (itk::ExceptionObject &excp) {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return ""; // EXIT_FAILURE;
  }
  std::cerr << "Alledgedly writing the series was successful to dir: "
            << newDirPath.toStdString() << std::endl;
  return newDirPath.toLocal8Bit().constData();
}

QString get_output_options(const UShortImageType::Pointer &m_spFixed) {

  const QString str_fixed_origin =
      QString("%1,%2,%3") // done per image because CT might be different from
                          // reconstructed CBCT
          .arg(m_spFixed->GetOrigin()[0])
          .arg(m_spFixed->GetOrigin()[1])
          .arg(m_spFixed->GetOrigin()[2]);
  const QString str_fixed_dimension =
      QString("%1,%2,%3")
          .arg(m_spFixed->GetBufferedRegion().GetSize()[0])
          .arg(m_spFixed->GetBufferedRegion().GetSize()[1])
          .arg(m_spFixed->GetBufferedRegion().GetSize()[2]);
  const QString str_fixed_spacing = QString("%1,%2,%3")
                                        .arg(m_spFixed->GetSpacing()[0])
                                        .arg(m_spFixed->GetSpacing()[1])
                                        .arg(m_spFixed->GetSpacing()[2]);
  const QString str_fixed_direction = QString("%1,%2,%3,%4,%5,%6,%7,%8,%9")
                                          .arg(m_spFixed->GetDirection()[0][0])
                                          .arg(m_spFixed->GetDirection()[0][1])
                                          .arg(m_spFixed->GetDirection()[0][2])
                                          .arg(m_spFixed->GetDirection()[1][0])
                                          .arg(m_spFixed->GetDirection()[1][1])
                                          .arg(m_spFixed->GetDirection()[1][2])
                                          .arg(m_spFixed->GetDirection()[2][0])
                                          .arg(m_spFixed->GetDirection()[2][1])
                                          .arg(m_spFixed->GetDirection()[2][2]);

  return QString(" --origin %1 --spacing %2 --dimension %3 --direction %4")
      .arg(str_fixed_origin)
      .arg(str_fixed_spacing)
      .arg(str_fixed_dimension)
      .arg(str_fixed_direction);
}

bool GetCouchShiftFromINIXVI(QString &strPathINIXVI, VEC3D *pTrans,
                             VEC3D *pRot) {
  QFileInfo fInfo(strPathINIXVI);
  if (!fInfo.exists()) {
    return false;
  }

  std::ifstream fin;
  fin.open(strPathINIXVI.toLocal8Bit().constData());

  if (fin.fail()) {
    return false;
  }

  char str[MAX_LINE_LENGTH];

  float couch_Lat_cm = 0.0;
  float couch_Long_cm = 0.0;
  float couch_Vert_cm = 0.0;

  float couch_Pitch = 0.0;
  float couch_Yaw = 0.0;
  float couch_Roll = 0.0;

  bool bFound = false;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString tmpStr = QString(&str[0]);
    QStringList strListParam = tmpStr.split("=");

    QString tagName, strVal;

    if (strListParam.count() == 2) {
      tagName = strListParam.at(0);
      strVal = strListParam.at(1);
      tagName = tagName.trimmed();
      strVal = strVal.trimmed();

      if (tagName == "CouchShiftLat") {
        couch_Lat_cm = strVal.toFloat();
        bFound = true;
      } else if (tagName == "CouchShiftLong") {
        couch_Long_cm = strVal.toFloat();
      } else if (tagName == "CouchShiftHeight") {
        couch_Vert_cm = strVal.toFloat();
      } else if (tagName == "CouchPitch") {
        couch_Pitch = strVal.toFloat();
      } else if (tagName == "CouchRoll") {
        couch_Yaw = strVal.toFloat();
      } else if (tagName == "CouchYaw") {
        couch_Roll = strVal.toFloat();
      }
    }
  }
  fin.close();

  if (!bFound) {
    return false;
  }

  // Warning!! dicom convention!
  pTrans->x = couch_Lat_cm * 10.0; // sign should be checked
  // pTrans->y = couch_Vert_cm*10.0; //sign should be checked // IEC-->DICOM is
  // already accounted for..but sign!
  pTrans->y = couch_Vert_cm * (-10.0); // consistent with Tracking software
  pTrans->z = couch_Long_cm * 10.0;    // sign should be checked

  pRot->x = couch_Pitch;
  pRot->y = couch_Yaw;
  pRot->z = couch_Roll;
  // x,y,z: dicom
  return true;
}

bool GetXrayParamFromINI(QString &strPathINI, float &kVp, float &mA,
                         float &ms) {
  QFileInfo info = QFileInfo(strPathINI);

  kVp = 0.0;
  mA = 0.0;
  ms = 0.0;

  if (!info.exists()) {
    return false;
  }

  // TubeMA=64.0000
  // TubeKV = 120.0000
  // TubeKVLength = 40.0000
  std::ifstream fin;
  fin.open(strPathINI.toLocal8Bit().constData());

  if (fin.fail()) {
    return false;
  }

  char str[MAX_LINE_LENGTH];

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString tmpStr = QString(&str[0]);
    QStringList strListParam = tmpStr.split("=");

    QString tagName;
    QString strVal;

    if (strListParam.count() == 2) {
      tagName = strListParam.at(0);
      strVal = strListParam.at(1);
      tagName = tagName.trimmed();
      strVal = strVal.trimmed();

      if (tagName == "TubeMA") {
        mA = strVal.toFloat();
      }
      if (tagName == "TubeKVLength") {
        ms = strVal.toFloat();
      }
      if (tagName == "TubeKV") {
        kVp = strVal.toFloat();
      }
    }
  }
  fin.close();

  return !(kVp == 0 || mA == 0 || ms == 0);
}

// Projection image Median filtering using CUDA

// Dir or File
bool LoadShortImageDirOrFile(QString &strPathDir,
                             ShortImageType::Pointer &spOutputShortImg) {
  QFileInfo fInfo(strPathDir);
  if (!fInfo.exists()) {
    return false;
  }

  Plm_image plmImg;
  plmImg.load_native(strPathDir.toLocal8Bit().constData());
  const ShortImageType::Pointer spShortImg = plmImg.itk_short();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(spShortImg);
  imageCalculatorFilter->Compute();

  /* double minVal0 = (double)(imageCalculatorFilter->GetMinimum());
  double maxVal0 = (double)(imageCalculatorFilter->GetMaximum());*/

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;

  ThresholdImageFilterType::Pointer thresholdFilterAbove =
      ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(spShortImg);
  thresholdFilterAbove->ThresholdAbove(3072);
  thresholdFilterAbove->SetOutsideValue(3072);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
      ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(-1024);
  thresholdFilterBelow->SetOutsideValue(-1024);
  thresholdFilterBelow->Update();

  spOutputShortImg = thresholdFilterBelow->GetOutput();
  std::cout << "Image file was loaded" << std::endl;

  return true;
}
