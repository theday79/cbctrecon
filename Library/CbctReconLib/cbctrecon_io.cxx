/*All IO functions used with cbctrecon*/
#include "cbctrecon.h"

// std
#include <iostream> // for operator<<, basic_o...
#include <string>   // for string
#include <vector>   // for vector

// Qt
#include <qxmlstream.h>

// ITK
#ifdef OF
#undef OF
#endif // OF

#include "gdcmReader.h"
#include "gdcmAttribute.h"
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
#undef TIMEOUT
#undef CUDA_FOUND
#include "dcmtk_rt_study.h"
#include "plm_image.h"
#include <rt_study_metadata.h>

// local
#include "StructureSet.h"
#include "cbctrecon_io.h"
#include "cbctrecon_types.h"

QString MakeElektaXML(const QString &filePath_ImageDBF,
                      const QString &filePath_FrameDBF,
                      const QString &DICOM_UID) {
  std::cout << "Elekta geometry XML file is being generated." << std::endl;
  // Define FRAME.DBF path
  auto reader = rtk::ElektaSynergyGeometryReader::New();
  // string strDicomUID = DICOM_UID.toLocal8Bit().constData();
  // //DICOM_UID.toStdString()  string strDicomUID = DICOM_UID.toStdString();
  // string strDbfImg = filePath_ImageDBF.toStdString();
  // string strDbfFrame = filePath_FrameDBF.toStdString();

  auto info = QFileInfo(filePath_ImageDBF);
  const auto dirPath = info.absolutePath();

  const auto fileName = "ElektaGeom_" + DICOM_UID + ".xml";

  auto str_output = dirPath + "/" + fileName;

  reader->SetDicomUID(DICOM_UID.toLocal8Bit().constData());
  reader->SetImageDbfFileName(filePath_ImageDBF.toLocal8Bit().constData());
  reader->SetFrameDbfFileName(filePath_FrameDBF.toLocal8Bit().constData());

  reader->UpdateOutputData();

  // Write
  auto xmlWriter = rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
  xmlWriter->SetFilename(str_output.toLocal8Bit().constData());
  xmlWriter->SetObject(reader->GetGeometry());
  xmlWriter->WriteFile();

  std::cout << "Reading succeed" << std::endl;

  return str_output;
}

FLEXDATA XML_parseFrameForXVI5(QXmlStreamReader &xml) {

  auto tmpResult =
      FLEXDATA{/*fGanAngle =*/0.0, /*fPanelOffsetX =*/0.0,
               /*fPanelOffsetY =*/0.0, /*bKV_On =*/true, /*bMV_On =*/false};

  /* Let's check that we're really getting a person. */
  if (xml.tokenType() != QXmlStreamReader::StartElement &&
      xml.name() == "Frame") {
    return tmpResult;
  }
  /* Let's get the attributes for person */
  // QXmlStreamAttributes attributes = xml.attributes();
  /* Let's check that person has id attribute. */
  // if (attributes.hasAttribute("id")) {
  //	/* We'll add it to the map. */
  //	person["id"] = attributes.value("id").toString();
  //}
  /* Next element... */
  xml.readNext();
  /*
   * We're going to loop over the things because the order might change.
   * We'll continue the loop until we hit an EndElement named person.
   */
  while (!(xml.tokenType() == QXmlStreamReader::EndElement &&
           xml.name() == "Frame")) {
    auto tmpXmlName = xml.name();
    auto strTmpXMLName = QString(tmpXmlName.toLocal8Bit().constData());
    // int tmpType = (int)(xml.tokenType());

    QString tmpStr;
    if (xml.tokenType() == QXmlStreamReader::StartElement) {
      /* We've found first name. */
      if (xml.name() == "Seq") {
        tmpStr = XML_GetSingleItemString(xml);
      }
      /* We've found surname. */
      else if (xml.name() == "DeltaMS") {
        tmpStr = XML_GetSingleItemString(xml);
      }
      /* We've found email. */
      else if (xml.name() == "HasPixelFactor") {
        tmpStr = XML_GetSingleItemString(xml);
      }
      /* We've found website. */
      else if (xml.name() == "PixelFactor") {
        tmpStr = XML_GetSingleItemString(xml);
      } else if (xml.name() == "GantryAngle") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fGanAngle = tmpStr.toFloat();
      } else if (xml.name() == "Exposed") {
        tmpStr = XML_GetSingleItemString(xml);
        if (tmpStr == "True") {
          tmpResult.bKV_On = true;
        } else {
          tmpResult.bKV_On = false;
        }
      }
      /*else if (xml.name() == "Exposed") {
          tmpStr = XML_GetSingleItemString(xml);
          if (tmpStr == "True")
              tmpResult.bKV_On = true;
          else
              tmpResult.bKV_On = false;
      }*/
      else if (xml.name() == "MVOn") {
        tmpStr = XML_GetSingleItemString(xml);
        if (tmpStr == "True") {
          tmpResult.bMV_On = true;
        } else {
          tmpResult.bMV_On = false;
        }
      } else if (xml.name() == "UCentre") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fPanelOffsetX = tmpStr.toFloat();
      } else if (xml.name() == "VCentre") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fPanelOffsetY = tmpStr.toFloat();
      } else if (xml.name() == "Inactive") {
        tmpStr = XML_GetSingleItemString(xml);
      }
    }
    xml.readNext();
  }
  return tmpResult;
}

// Get the projection geometry
void CbctRecon::LoadRTKGeometryFile(const char *filePath) {
  auto geometryReader =
      rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(filePath);
  geometryReader->GenerateOutputInformation();
  std::cout << "Geometry reading succeed" << std::endl;

  m_spFullGeometry = geometryReader->GetOutputObject();

  // fullGeometry->GetGantryAngles();
  const auto geoDataSize =
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

  const auto itBegin = m_spFullGeometry->GetGantryAngles().begin();
  const auto itEnd = m_spFullGeometry->GetGantryAngles().end();

  for (auto it = itBegin; it != itEnd; ++it) {
    auto tmpAngle = *it;

    if (tmpAngle > 180.0) {
      tmpAngle = tmpAngle - 360.0;
    }

    vTempConvAngles.push_back(tmpAngle);
  }

  // compare 2 points in the middle of the angle list
  const auto iLowerIdx = static_cast<size_t>(geoDataSize * 1.0 / 3.0);
  const auto iUpperIdx = static_cast<size_t>(geoDataSize * 2.0 / 3.0);

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
  auto reader = ReaderType::New();

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
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(reader->GetOutput());
  imageCalculatorFilter->Compute();

  const auto minVal0 = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal0 = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Original Min and Max Values are	" << minVal0 << "	"
            << maxVal0 << std::endl;

  auto bNKI = false;
  if (minVal0 > -600) // impossible for normal Short image. IN NKI, always -512.
                      // don't know why
  {
    bNKI = true;
  }

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  auto thresholdFilter = ThresholdImageFilterType::New();

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
  auto spRescaleFilter = RescaleFilterType::New();
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
  auto duplicator = DuplicatorType::New();
  duplicator->SetInputImage(spUsImage);
  duplicator->Update();
  const UShortImageType::Pointer clonedReconImage = duplicator->GetOutput();
  ShortImageType::Pointer clonedReconImageSHORT;

  using ThresholdImageFilterType = itk::ThresholdImageFilter<UShortImageType>;
  auto thresholdFilterAbove = ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(clonedReconImage);
  thresholdFilterAbove->ThresholdAbove(4095);
  thresholdFilterAbove->SetOutsideValue(4095);

  auto thresholdFilterBelow = ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(0);
  thresholdFilterBelow->SetOutsideValue(0);
  thresholdFilterBelow->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<UShortImageType>;
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilterBelow->GetOutput());
  imageCalculatorFilter->Compute();
  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  const auto outputMinVal = static_cast<SHORT_PixelType>(minVal - 1024);
  const auto outputMaxVal = static_cast<SHORT_PixelType>(maxVal - 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType>;
  auto spRescaleFilter = RescaleFilterType::New();
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
  auto writer = WriterType::New();
  writer->SetFileName(outputFilePath.toLocal8Bit().constData());
  // writer->SetUseCompression(true);
  writer->SetUseCompression(false); // for plastimatch
  writer->SetInput(spRescaleFilter->GetOutput());

  std::cout << "Writing is under progress...: "
            << outputFilePath.toLocal8Bit().constData() << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;
}

DCM_MODALITY get_dcm_modality(const QString& filename){
  gdcm::Reader reader;
  reader.SetFileName(filename.toLocal8Bit().constData());
  if (!reader.Read())
  {
    std::cerr << "Reading dicom: " << filename.toStdString() << " failed!\n";
    return RTUNKNOWN;
  }
  gdcm::File &file = reader.GetFile();
  gdcm::DataSet &ds = file.GetDataSet();
  gdcm::Attribute<0x0008, 0x0060> at_modality;
  at_modality.SetFromDataElement(ds.GetDataElement(at_modality.GetTag()));
  const auto modality = at_modality.GetValue();
  if (modality.compare("RTIMAGE") == 0 ||
      modality.compare("CT") == 0){
    return RTIMAGE;
  }
  if (modality.compare("RTDOSE") == 0){
    return RTDOSE;
  }
  if (modality.compare("RTSTRUCT") == 0){
    return RTSTRUCT;
  }
  if (modality.compare("RTPLAN") == 0){
    return RTPLAN;
  }
  if (modality.compare("RTRECORD") == 0){
    return RTRECORD;
  }
  std::cerr << "Modality was: " << modality << "\n";
  return RTUNKNOWN;
}

std::unique_ptr<Rtss_modern> load_rtstruct(const QString& filename){

  auto reader = gdcm::Reader();
  reader.SetFileName(filename.toLocal8Bit().constData());
  if (!reader.Read())
  {
    std::cerr << "Reading dicom rtstruct: " << filename.toStdString() << " failed!\n";
    return nullptr;
  }

  gdcm::File &file = reader.GetFile();
  gdcm::DataSet &ds = file.GetDataSet();

  gdcm::Attribute<0x0008, 0x0060> at_modality;
  at_modality.SetFromDataElement(ds.GetDataElement(at_modality.GetTag()));
  const auto modality = at_modality.GetValue();
  if (modality != "RTSTRUCT"){
      std::cerr << "Modality was not RTSTRUCT, it was: " << modality << "\n";
      return nullptr;
  }

  auto rt_struct = std::make_unique<Rtss_modern>();
  rt_struct->num_structures = 0;

  const gdcm::DataElement &roi_seq_tag = ds.GetDataElement(gdcm::Tag(0x3006, 0x0020));
  auto roi_seq = roi_seq_tag.GetValueAsSQ();
  for (auto it_roi = roi_seq->Begin(); it_roi != roi_seq->End(); ++it_roi){
    auto rt_roi = std::make_unique<Rtss_roi_modern>();
    auto at_roi_number = gdcm_attribute_from<0x3006, 0x0022>(it_roi);
    rt_roi->id = static_cast<size_t>(at_roi_number.GetValue());

    auto at_roi_name = gdcm_attribute_from<0x3006, 0x0026>(it_roi);
    rt_roi->name = at_roi_name.GetValue();

    rt_struct->slist.emplace_back(std::move(rt_roi));
    rt_struct->num_structures++;

  }

  const gdcm::DataElement &roi_contour_seq_tag = ds.GetDataElement(gdcm::Tag(0x3006, 0x0039));
  auto roi_contour_seq = roi_contour_seq_tag.GetValueAsSQ();
  auto i = 0U;
  for (auto it_roi_contour = roi_contour_seq->Begin(); it_roi_contour != roi_contour_seq->End(); ++it_roi_contour){
    auto at_roi_contour_number = gdcm_attribute_from<0x3006, 0x0084>(it_roi_contour);
    if (static_cast<int>(rt_struct->slist.at(i).id) != at_roi_contour_number.GetValue()){
      std::cerr << "ID mismatch: " << rt_struct->slist.at(i).id << " vs " << at_roi_contour_number.GetValue() << "\n"
                << "There might be something wrong with " << rt_struct->slist.at(i).name << "\n"
                << "Caution! As we continue anyway...\n";
    }
    auto at_roi_contour_colour = gdcm_attribute_from<0x3006, 0x002A>(it_roi_contour);
    const auto color = at_roi_contour_colour.GetValues();
    auto s_color = std::to_string(color[0]) + " "
                 + std::to_string(color[1]) + " "
                 + std::to_string(color[2]);
    rt_struct->slist.at(i).color = s_color.c_str();

    const auto& contour_seq_tag = it_roi_contour->GetDataElement(gdcm::Tag(0x3006, 0x0040));
    auto contour_seq = contour_seq_tag.GetValueAsSQ();
    auto j = 0U;
    rt_struct->slist.at(i).num_contours = contour_seq->GetLength();
    rt_struct->slist.at(i).pslist.resize(contour_seq->GetLength());
    for (auto it_contour = contour_seq->Begin(); it_contour != contour_seq->End(); ++it_contour){
      auto rt_contour = std::make_unique<Rtss_contour_modern>();

      auto at_contour_number_of_points = gdcm_attribute_from<0x3006, 0x0046>(it_contour);
      rt_contour->num_vertices = static_cast<unsigned long>(at_contour_number_of_points.GetValue());

      auto at_contour_points = gdcm_attribute_from<0x3006, 0x0050>(it_contour);
      const auto points = at_contour_points.GetValues();

      rt_contour->coordinates.resize(rt_contour->num_vertices);

      std::generate(std::begin(rt_contour->coordinates), std::end(rt_contour->coordinates),
                    [&points, k = 0]() mutable {
         auto vec = FloatVector{
                 static_cast<float>(points[k + 0]),
                 static_cast<float>(points[k + 1]),
                 static_cast<float>(points[k + 2])};
         k += 3;
         return vec;
      });


      rt_struct->slist.at(i).pslist.at(j++) = std::move(rt_contour);
    }

    i++;
  }


  return rt_struct;
}


bool CbctRecon::ReadDicomDir(QString &dirPath) {
  auto filenamelist = std::vector<std::string>();
  auto dir = QDir(dirPath);
  for (auto&& filename : dir.entryList(QStringList() << "*.dcm" << "*.DCM", QDir::Files)){
    const auto fullfilename = dir.absolutePath() + "/" + filename;
    auto modality = get_dcm_modality(fullfilename);
    switch (modality) {
    case RTIMAGE:
    case RTDOSE:
      filenamelist.push_back(fullfilename.toStdString());
      break;
    case RTSTRUCT:
      m_structures->set_planCT_ss(load_rtstruct(fullfilename));
      break;
    case RTPLAN:
      break; // Maybe some pre-loading for gPMC could be useful?
    case RTRECORD:
      break; // I haven't ever seen one IRL
    case RTUNKNOWN:
      std::cerr << "File: " << fullfilename.toStdString() << " was not of a recognizeable modality type!\n";
      break;
    }
  }

  ShortImageType::Pointer spShortImg;

  if (filenamelist.size() != 0){
    using dcm_reader_type = itk::ImageSeriesReader<ShortImageType>;
    auto dcm_reader = dcm_reader_type::New();
    dcm_reader->SetFileNames(filenamelist);
    dcm_reader->Update();
    spShortImg = dcm_reader->GetOutput();
  }
  else{
    Dcmtk_rt_study drs(dirPath.toLocal8Bit().constData());
    drs.load_directory(); // parse_directory();

    Plm_image plmImg;
    auto tmp_img = drs.get_image();

    if (!tmp_img){
        std::cerr << "Plastimach couldn't read image data!\n";
        return false;
    }
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

    spShortImg = plmImg.itk_short();
  }



  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  auto thresholdFilter = ThresholdImageFilterType::New();

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
  auto spRescaleFilter = RescaleFilterType::New();
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
  auto thresholdFilter = ThresholdImageFilterType::New();
  thresholdFilter->SetInput(spImgUshort);
  thresholdFilter->ThresholdOutside(0, 4096); //--> 0 ~ 4095
  thresholdFilter->SetOutsideValue(0);
  thresholdFilter->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<UShortImageType>;
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();
  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  const auto outputMinVal = static_cast<SHORT_PixelType>(minVal - 1024);
  const auto outputMaxVal = static_cast<SHORT_PixelType>(maxVal - 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType>;
  auto spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  spImgShort = spRescaleFilter->GetOutput();
}

void ConvertShort2Ushort(ShortImageType::Pointer &spInputImgShort,
                         UShortImageType::Pointer &spOutputImgUshort) {
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  auto thresholdFilterAbove = ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(spInputImgShort);
  thresholdFilterAbove->ThresholdAbove(3071);
  thresholdFilterAbove->SetOutsideValue(3071);

  auto thresholdFilterBelow = ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(-1024);
  thresholdFilterBelow->SetOutsideValue(-1024);
  thresholdFilterBelow->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();
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
  auto spRescaleFilter = RescaleFilterType::New();
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

  auto newDirPath = strPathTargetDir + "/" + strPatientID + "_DCM";
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

  auto newDirPath =
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

  auto region = spShortImg->GetLargestPossibleRegion();
  auto start = region.GetIndex();
  auto size = region.GetSize();

  auto gdcmIO = ImageIOType::New();
  auto &dict = gdcmIO->GetMetaDataDictionary();
  std::string value = "CT";
  itk::EncapsulateMetaData<std::string>(dict, "0008|0060", value); // Modality
  value = "DERIVED\\SECONDARY\\AXIAL"; // This is virtually always correct when
                                       // using ITK to write an image
  itk::EncapsulateMetaData<std::string>(dict, "0008|0008", value); // Image Type
  value = "SI";
  itk::EncapsulateMetaData<std::string>(dict, "0008|0064",
                                        value); // Conversion Type
  const auto value_double = spShortImg->GetSpacing()[2];
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

  auto namesGenerator = NamesGeneratorType::New();
  namesGenerator->SetStartIndex(static_cast<itk::SizeValueType>(start[2]));
  namesGenerator->SetEndIndex(static_cast<itk::SizeValueType>(start[2]) +
                              size[2] - 1);
  namesGenerator->SetIncrementIndex(1);
  namesGenerator->SetSeriesFormat(newDirPath.toStdString() + "/CT." + studyUID +
                                  ".%d.dcm");

  using SeriesWriterType =
      itk::ImageSeriesWriter<ShortImageType, OutputImageType>;
  auto seriesWriter = SeriesWriterType::New();
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

  const auto str_fixed_origin =
      QString("%1,%2,%3") // done per image because CT might be different from
                          // reconstructed CBCT
          .arg(m_spFixed->GetOrigin()[0])
          .arg(m_spFixed->GetOrigin()[1])
          .arg(m_spFixed->GetOrigin()[2]);
  const auto str_fixed_dimension =
      QString("%1,%2,%3")
          .arg(m_spFixed->GetBufferedRegion().GetSize()[0])
          .arg(m_spFixed->GetBufferedRegion().GetSize()[1])
          .arg(m_spFixed->GetBufferedRegion().GetSize()[2]);
  const auto str_fixed_spacing = QString("%1,%2,%3")
                                     .arg(m_spFixed->GetSpacing()[0])
                                     .arg(m_spFixed->GetSpacing()[1])
                                     .arg(m_spFixed->GetSpacing()[2]);
  const auto str_fixed_direction = QString("%1,%2,%3,%4,%5,%6,%7,%8,%9")
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

  auto bFound = false;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    auto tmpStr = QString(&str[0]);
    auto strListParam = tmpStr.split("=");

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
  pTrans->x = static_cast<double>(couch_Lat_cm) * 10.0; // sign should be checked
  // pTrans->y = couch_Vert_cm*10.0; //sign should be checked // IEC-->DICOM is
  // already accounted for..but sign!
  pTrans->y = static_cast<double>(couch_Vert_cm) * -10.0; // consistent with Tracking software
  pTrans->z = static_cast<double>(couch_Long_cm) * 10.0;  // sign should be checked

  pRot->x = static_cast<double>(couch_Pitch);
  pRot->y = static_cast<double>(couch_Yaw);
  pRot->z = static_cast<double>(couch_Roll);
  // x,y,z: dicom
  return true;
}

bool GetXrayParamFromINI(QString &strPathINI, float &kVp, float &mA,
                         float &ms) {
  auto info = QFileInfo(strPathINI);

  kVp = 0.0f;
  mA = 0.0f;
  ms = 0.0f;

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
    auto tmpStr = QString(&str[0]);
    auto strListParam = tmpStr.split("=");

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

  return !(kVp == 0.0f || mA == 0.0f || ms == 0.0f);
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
  const auto spShortImg = plmImg.itk_short();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(spShortImg);
  imageCalculatorFilter->Compute();

  /* double minVal0 = (double)(imageCalculatorFilter->GetMinimum());
  double maxVal0 = (double)(imageCalculatorFilter->GetMaximum());*/

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;

  auto thresholdFilterAbove = ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(spShortImg);
  thresholdFilterAbove->ThresholdAbove(3072);
  thresholdFilterAbove->SetOutsideValue(3072);

  auto thresholdFilterBelow = ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(-1024);
  thresholdFilterBelow->SetOutsideValue(-1024);
  thresholdFilterBelow->Update();

  spOutputShortImg = thresholdFilterBelow->GetOutput();
  std::cout << "Image file was loaded" << std::endl;

  return true;
}
