// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

/*All IO functions used with cbctrecon*/
#include "cbctrecon.h"

// std
#include <filesystem>
#include <iostream> // for operator<<, basic_o...
#include <string>   // for string
#include <vector>   // for vector

// Qt
#include <qxmlstream.h>

// ITK
#ifdef OF
#undef OF
#endif // OF

#include "gdcmAttribute.h"
#include "gdcmReader.h"
#include "gdcmUIDGenerator.h"
#include "gdcmWriter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkThresholdImageFilter.h"

// RTK
#include "rtkElektaSynergyGeometryReader.h"
#include "rtkProjectionsReader.h"
#include "rtkThreeDCircularProjectionGeometryXMLFileReader.h" // for ThreeDCircularProje...
#include "rtkThreeDCircularProjectionGeometryXMLFileWriter.h" // for ThreeDCircularProje...

// PLM
#undef TIMEOUT
#undef CUDA_FOUND
#include "dcmtk_rt_study.h"
#include "dicom_util.h"
#include "plm_image.h"
#include <rt_study_metadata.h>

// DCMTK
#undef NEUTRAL
#include "dcmtk/dcmrt/drtstrct.h"

// local
#include "StructureSet.h"
#include "cbctrecon_io.h"
#include "cbctrecon_types.h"

std::filesystem::path
MakeElektaXML(const std::filesystem::path &filePath_ImageDBF,
              const std::filesystem::path &filePath_FrameDBF,
              const std::string &DICOM_UID) {
  std::cout << "Elekta geometry XML file is being generated." << std::endl;
  // Define FRAME.DBF path
  auto reader = rtk::ElektaSynergyGeometryReader::New();
  // string strDicomUID = DICOM_UID.toLocal8Bit().constData();
  // //DICOM_UID.toStdString()  string strDicomUID = DICOM_UID.toStdString();
  // string strDbfImg = filePath_ImageDBF.toStdString();
  // string strDbfFrame = filePath_FrameDBF.toStdString();

  const auto dirPath =
      std::filesystem::absolute(filePath_ImageDBF.parent_path());

  const auto fileName = "ElektaGeom_" + DICOM_UID + ".xml";

  auto str_output = dirPath / fileName;

  reader->SetDicomUID(DICOM_UID);
  reader->SetImageDbfFileName(filePath_ImageDBF.string());
  reader->SetFrameDbfFileName(filePath_FrameDBF.string());

  reader->UpdateOutputData();

  // Write
  auto xmlWriter = rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
  xmlWriter->SetFilename(str_output.string());
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
    // auto strTmpXMLName = std::string(tmpXmlName.toLocal8Bit().constData());
    // int tmpType = (int)(xml.tokenType());

    std::string tmpStr;
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
        tmpResult.fGanAngle = from_string<float>(tmpStr).value_or(0.0f);
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
        tmpResult.fPanelOffsetX = from_string<float>(tmpStr).value_or(0.0f);
      } else if (xml.name() == "VCentre") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fPanelOffsetY = from_string<float>(tmpStr).value_or(0.0f);
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
  const auto angles = m_spFullGeometry->GetGantryAngles();
  const auto geoDataSize = angles.size(); // This is MV gantry angle!!!
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

  const auto itBegin = angles.begin();
  const auto itEnd = angles.end();

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

bool LoadShortImageToUshort(std::string &strPath,
                            UShortImageType::Pointer &pUshortImage) {
  using ReaderType = itk::ImageFileReader<ShortImageType>;
  auto reader = ReaderType::New();

  // std::string fileName = QFileDialog::getOpenFileName(this, "Open Image","",
  // "Plan CT file (*.mha)",0,0);

  if (strPath.length() < 1) {
    return false;
  }

  reader->SetFileName(strPath);
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

  USHORT_PixelType outputMinVal;
  USHORT_PixelType outputMaxVal;
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
                         std::string &outputFilePath) {
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
  writer->SetFileName(outputFilePath);
  // writer->SetUseCompression(true);
  writer->SetUseCompression(false); // for plastimatch
  writer->SetInput(spRescaleFilter->GetOutput());

  std::cout << "Writing is under progress...: " << outputFilePath << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;
}

DCM_MODALITY get_dcm_modality(const std::string &filename) {
  gdcm::Reader reader;
  reader.SetFileName(filename.c_str());
  if (!reader.Read()) {
    std::cerr << "Reading dicom: " << filename << " failed!\n";
    return DCM_MODALITY::RTUNKNOWN;
  }
  auto &file = reader.GetFile();
  auto &ds = file.GetDataSet();
  gdcm::Attribute<0x0008, 0x0060> at_modality;
  at_modality.SetFromDataElement(
      ds.GetDataElement(gdcm::Attribute<0x8, 0x60>::GetTag()));
  const auto modality = at_modality.GetValue();
  if (modality == "RTIMAGE" || modality == "CT" || modality == "MR") {
    return DCM_MODALITY::RTIMAGE;
  }
  if (modality == "RTDOSE") {
    return DCM_MODALITY::RTDOSE;
  }
  if (modality == "RTSTRUCT") {
    return DCM_MODALITY::RTSTRUCT;
  }
  if (modality == "RTPLAN") {
    return DCM_MODALITY::RTPLAN;
  }
  if (modality == "RTRECORD") {
    return DCM_MODALITY::RTRECORD;
  }
  std::cerr << "Modality was: " << modality << "\n";
  return DCM_MODALITY::RTUNKNOWN;
}

bool check_rtss_dicom_integrity(const gdcm::DataSet &ds) {
  const auto trefframerefsq = gdcm::Tag(0x3006, 0x0010);
  if (!ds.FindDataElement(trefframerefsq)) {
    return false;
  }
  const auto &refframerefsq = ds.GetDataElement(trefframerefsq);
  auto sqi0 = refframerefsq.GetValueAsSQ();
  if (!sqi0 || !sqi0->GetNumberOfItems()) {
    return false;
  }
  for (unsigned int pd = 0; pd < sqi0->GetNumberOfItems(); ++pd) {
    const auto &item0 = sqi0->GetItem(pd + 1); // Item start at #1
    const auto &nestedds0 = item0.GetNestedDataSet();
    // RTReferencedStudySequence

    gdcm::Tag trtrefstudysq(0x3006, 0x0012);
    if (!nestedds0.FindDataElement(trtrefstudysq)) {
      return false;
    }
    const auto &rtrefstudysq = nestedds0.GetDataElement(trtrefstudysq);
    auto sqi00 = rtrefstudysq.GetValueAsSQ();
    if (!sqi00 || !sqi00->GetNumberOfItems()) {
      return false;
    }
    assert(sqi00->GetNumberOfItems() == 1);
    for (unsigned int pd0 = 0; pd0 < sqi00->GetNumberOfItems(); ++pd0) {
      const auto &item00 = sqi00->GetItem(pd0 + 1); // Item start at #1
      const auto &nestedds00 = item00.GetNestedDataSet();

      // RTReferencedSeriesSequence
      auto trtrefseriessq = gdcm::Tag(0x3006, 0x0014);
      if (!nestedds00.FindDataElement(trtrefseriessq)) {
        return false;
      }
      const auto &rtrefseriessq = nestedds00.GetDataElement(trtrefseriessq);

      auto sqi000 = rtrefseriessq.GetValueAsSQ();
      if (!sqi000 || !sqi000->GetNumberOfItems()) {
        return false;
      }
      assert(sqi000->GetNumberOfItems() == 1);
      for (unsigned int pd00 = 0; pd00 < sqi000->GetNumberOfItems(); ++pd00) {
        const auto &item000 = sqi000->GetItem(pd00 + 1); // Item start at #1
        const auto &nestedds000 = item000.GetNestedDataSet();

        // ContourImageSequence
        gdcm::Tag tcontourimageseq(0x3006, 0x0016);
        if (!nestedds000.FindDataElement(tcontourimageseq)) {
          return false;
        }
        const auto &contourimageseq =
            nestedds000.GetDataElement(tcontourimageseq);
        const auto sqi0000 = contourimageseq.GetValueAsSQ();
        if (!sqi0000 || !sqi0000->GetNumberOfItems()) {
          return false;
        }
      }
    }
  }
  return true;
}

// From GDCM vtkGDCMPolyDataReader, with modification:
std::unique_ptr<Rtss_modern>
RequestData_RTStructureSetStorage(gdcm::Reader const &reader) {

  auto rt_struct = std::make_unique<Rtss_modern>();

  const auto &ds = reader.GetFile().GetDataSet();
  // (3006,0010) SQ (Sequence with undefined length #=1)     # u/l, 1
  // ReferencedFrameOfReferenceSequence (3006,0020) SQ (Sequence with explicit
  // length #=4)      # 370, 1 StructureSetROISequence (3006,0039) SQ (Sequence
  // with explicit length #=4)      # 24216, 1 ROIContourSequence
  const auto troicsq = gdcm::Tag(0x3006, 0x0039);
  if (!ds.FindDataElement(troicsq)) {
    return nullptr;
  }
  const auto tssroisq = gdcm::Tag(0x3006, 0x0020);
  if (!ds.FindDataElement(tssroisq)) {
    return nullptr;
  }

  const auto input_good = check_rtss_dicom_integrity(ds);
  if (!input_good) {
    std::cerr << "Rtss dicom input does not contain required fields\n";
    return nullptr;
  }

  const auto &roicsq = ds.GetDataElement(troicsq);
  // std::cout << roicsq << std::endl;
  // const gdcm::SequenceOfItems *sqi_debug = roicsq.GetSequenceOfItems();
  auto sqi = roicsq.GetValueAsSQ();
  if (!sqi || !sqi->GetNumberOfItems()) {
    return nullptr;
  }
  const auto &ssroisq = ds.GetDataElement(tssroisq);
  // const gdcm::SequenceOfItems *ssqi = ssroisq.GetSequenceOfItems();
  auto ssqi = ssroisq.GetValueAsSQ();
  if (!ssqi || !ssqi->GetNumberOfItems()) {
    return nullptr;
  }

  rt_struct->num_structures = sqi->GetNumberOfItems();
  rt_struct->slist.resize(rt_struct->num_structures);

  // For each Item in the DataSet create a vtkPolyData
  for (unsigned int pd = 0; pd < sqi->GetNumberOfItems(); ++pd) {
    // StructureSetROI structuresetroi;

    const auto &item = sqi->GetItem(pd + 1); // Item start at #1
    // std::cout << item << std::endl;
    const auto &sitem = ssqi->GetItem(pd + 1); // Item start at #1
    const auto &snestedds = sitem.GetNestedDataSet();
    // (3006,0026) ?? (LO) [date]                                    # 4,1 ROI
    // Name
    gdcm::Tag stcsq(0x3006, 0x0026);
    if (!snestedds.FindDataElement(stcsq)) {
      continue;
    }
    const auto &sde = snestedds.GetDataElement(stcsq);
    std::string s(sde.GetByteValue()->GetPointer(),
                  sde.GetByteValue()->GetLength());
    // structuresetroi.ROIName = s;
    auto roinumber = gdcm::Attribute<0x3006, 0x0022>();
    roinumber.SetFromDataSet(snestedds);
    gdcm::Attribute<0x3006, 0x0026> roiname;
    roiname.SetFromDataSet(snestedds);
    assert(s == roiname.GetValue());

    auto &rt_roi = rt_struct->slist.at(pd);
    rt_roi.name = roiname.GetValue();
    rt_roi.id = roinumber.GetValue();

    const auto &nestedds = item.GetNestedDataSet();
    // std::cout << nestedds << std::endl;
    //(3006,002a) IS [255\192\96]                              # 10,3 ROI
    // Display Color
    gdcm::Tag troidc(0x3006, 0x002a);
    auto color = gdcm::Attribute<0x3006, 0x002a>();
    auto hasColor =
        false; // so that color[0] isn't referenced if the color isn't present.
    if (nestedds.FindDataElement(troidc)) {
      const auto &decolor = nestedds.GetDataElement(troidc);
      color.SetFromDataElement(decolor);
      hasColor = true;
      // std::cout << "color: " << roinumber.GetValue() << " -> " << color[0] <<
      // "," << color[1] << "," << color[2] << std::endl;
    }
    if (hasColor) {
      rt_struct->slist.at(pd).color = std::to_string(color[0]) + " " +
                                      std::to_string(color[1]) + " " +
                                      std::to_string(color[2]);
    } else {
      rt_roi.color = "255 0 0";
    }
    //(3006,0040) SQ (Sequence with explicit length #=8)      # 4326, 1
    // ContourSequence
    gdcm::Tag tcsq(0x3006, 0x0040);
    if (!nestedds.FindDataElement(tcsq)) {
      // FIXME: What if a contour sequence is empty but the color is set to
      // -say- 0/255/0 Since we are skipping entirely the contour sequence (no
      // vtkCellArray) we will not save the color.  which means it will be
      // reported as 0/0/0 in the output DICOM file.
      continue;
    }
    const auto &csq = nestedds.GetDataElement(tcsq);
    // std::cout << csq << std::endl;

    // const gdcm::SequenceOfItems *sqi2 = csq.GetSequenceOfItems();
    auto sqi2 = csq.GetValueAsSQ();
    if (!sqi2) //|| !sqi2->GetNumberOfItems() )
    {
      continue;
    }
    const auto nitems = sqi2->GetNumberOfItems();
    // std::cout << nitems << std::endl;
    // this->SetNumberOfOutputPorts(nitems);

    if (nitems == 0) {
      continue;
    }
    rt_roi.pslist.resize(nitems);
    rt_roi.num_contours = nitems;

    for (unsigned int ii = 0; ii < nitems; ++ii) {
      const auto &item2 = sqi2->GetItem(ii + 1); // Item start at #1

      const auto &nestedds2 = item2.GetNestedDataSet();
      // std::cout << nestedds2 << std::endl;
      // (3006,0050) DS
      // [43.57636\65.52504\-10.0\46.043102\62.564945\-10.0\49.126537\60.714...
      // # 398,48 ContourData
      gdcm::Tag tcontourdata(0x3006, 0x0050);
      const auto &contourdata = nestedds2.GetDataElement(tcontourdata);
      // std::cout << contourdata << std::endl;

      // const gdcm::ByteValue *bv = contourdata.GetByteValue();
      gdcm::Attribute<0x3006, 0x0042> contgeotype;
      contgeotype.SetFromDataSet(nestedds2);
      assert(contgeotype.GetValue() == "CLOSED_PLANAR " ||
             contgeotype.GetValue() == "POINT " ||
             contgeotype.GetValue() == "OPEN_NONPLANAR");

      auto numcontpoints = gdcm::Attribute<0x3006, 0x0046>();
      numcontpoints.SetFromDataSet(nestedds2);

      if (contgeotype.GetValue() == "POINT ") {
        assert(numcontpoints.GetValue() == 1);
      }

      gdcm::Attribute<0x3006, 0x0050> at;
      at.SetFromDataElement(contourdata);

      if (contgeotype.GetValue() == "CLOSED_PLANAR " ||
          contgeotype.GetValue() == "OPEN_NONPLANAR") {
        // http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.8.6.html
        if (nestedds2.FindDataElement(gdcm::Tag(0x3006, 0x0016))) {
          const auto &contourimagesequence =
              nestedds2.GetDataElement(gdcm::Tag(0x3006, 0x0016));
          const auto contourimagesequence_sqi =
              contourimagesequence.GetValueAsSQ();
          assert(contourimagesequence_sqi &&
                 contourimagesequence_sqi->GetNumberOfItems() == 1);

          /*gdcm::Attribute<0x0008, 0x1150> classat;
          classat.SetFromDataSet(thenestedds);
          gdcm::Attribute<0x0008, 0x1155> instat;
          instat.SetFromDataSet(thenestedds);*/
        }
      }

      // newPts->SetNumberOfPoints( at.GetNumberOfValues() / 3 );
      // assert( at.GetNumberOfValues() % 3 == 0); // FIXME
      const auto pts = at.GetValues();
      const auto npts = at.GetNumberOfValues() / 3;
      assert(npts == static_cast<unsigned int>(numcontpoints.GetValue()));
      assert(npts * 3 == at.GetNumberOfValues());
      auto &rt_contour = rt_roi.pslist.at(ii);
      rt_contour.num_vertices = npts;
      rt_contour.coordinates.reserve(npts);

      for (size_t i = 0; i < npts * 3; i += 3) {
        auto vertix = FloatVector{static_cast<float>(pts[i + 0]),
                                  static_cast<float>(pts[i + 1]),
                                  static_cast<float>(pts[i + 2])};

        rt_contour.coordinates.push_back(vertix);
      }
      // Each Contour Data is in fact a Cell:
    }
  }
  rt_struct->ready = true;

  return rt_struct;
}

bool AlterData_RTStructureSetStorage(const QFile &input_file,
                                     const Rtss_modern *input_rt_struct,
                                     const QFile &output_file) {

  DRTStructureSetIOD rtstruct;
  DcmFileFormat fileformat;
  auto status =
      fileformat.loadFile(input_file.fileName().toLocal8Bit().constData());
  if (!status.good()) {
    std::cerr << "Could not open RT struct dcm file: "
              << input_file.fileName().toStdString() << "\n";
    return false;
  }
  status = rtstruct.read(*fileformat.getDataset());
  if (!status.good()) {
    std::cerr << "Could not read RT struct dcm file: "
              << input_file.fileName().toStdString() << "\n";
    return false;
  }
  fileformat.clear();
  //  Change data in rtstruct
  //  ROI contour seq: 3006, 0039
  auto &ss_seq = rtstruct.getStructureSetROISequence();
  ss_seq.gotoFirstItem();
  // ROI Structure Set Seq: 3006, 0020
  auto &roi_seq = rtstruct.getROIContourSequence();
  roi_seq.gotoFirstItem();
  for (auto &rt_roi : input_rt_struct->slist) {
    auto &item = roi_seq.getCurrentItem();
    auto &ss_item = ss_seq.getCurrentItem();
    OFString roi_name;
    // ROI name: 3006, 0026
    ss_item.getROIName(roi_name);
    const auto trimmed_roi_name = trim_string(rt_roi.name);
    const auto trimmed_roi_name_dcm =
        trim_string(std::string(roi_name.c_str()));
    if (trimmed_roi_name.find("WEPL") != std::string::npos) {
      roi_seq.gotoNextItem();
      ss_seq.gotoNextItem();
      continue;
    }
    std::cerr << "Writing " << rt_roi.name << " to dicom file!\n";
    ss_item.setROIName(rt_roi.name.c_str());
    // Contour Seq: 3006, 0040
    auto &contour_seq = item.getContourSequence();
    contour_seq.gotoFirstItem();
    for (auto &rt_contour : rt_roi.pslist) {
      auto &contour = contour_seq.getCurrentItem();

      auto data_str = std::string("");
      for (auto &coord : rt_contour.coordinates) {
        data_str += std::to_string(coord.x) + "\\" + std::to_string(coord.y) +
                    "\\" + std::to_string(coord.z) + "\\";
      }
      data_str.pop_back();
      // Contour data: 3006, 0050
      status = contour.setContourData(OFString(data_str.c_str()), true);
      if (!status.good()) {
        std::cerr << "Could not set contour data: " << status.text() << "\n";
      }
      contour_seq.gotoNextItem();
    }
    roi_seq.gotoNextItem();
    ss_seq.gotoNextItem();
  }
  //
  status = rtstruct.write(*fileformat.getDataset());
  if (!status.good()) {
    std::cerr << "Could not write RT struct dcm file: "
              << input_file.fileName().toStdString() << "\n";
    return false;
  }
  status =
      fileformat.saveFile(output_file.fileName().toLocal8Bit().constData());
  if (!status.good()) {
    std::cerr << "Could not save RT struct dcm file: "
              << output_file.fileName().toStdString() << "\n";
    return false;
  }
  return true;
}

std::unique_ptr<Rtss_modern> load_rtstruct(const std::string &filename) {

  auto reader = gdcm::Reader();
  reader.SetFileName(filename.c_str());
  if (!reader.Read()) {
    std::cerr << "Reading dicom rtstruct: " << filename << " failed!\n";
    return nullptr;
  }

  auto rt_struct = RequestData_RTStructureSetStorage(reader);
  if (rt_struct == nullptr) {
    std::cerr << "Could not read RT structures!!\n";
  }

  return rt_struct;
}

double get_dcm_offset_z(const std::string &filename) {
  gdcm::Reader reader;
  reader.SetFileName(filename.c_str());
  if (!reader.Read()) {
    std::cerr << "Reading dicom: " << filename << " failed!\n";
    return std::numeric_limits<double>::min();
  }
  auto &file = reader.GetFile();
  auto &ds = file.GetDataSet();
  gdcm::Attribute<0x0020, 0x0032> at_position;
  at_position.SetFromDataElement(
      ds.GetDataElement(gdcm::Attribute<0x20, 0x32>::GetTag()));
  const auto position = at_position.GetValues();
  return position[2];
}

std::vector<std::string> get_dcm_image_files(std::filesystem::path &dir) {

  using NamesGeneratorType = itk::GDCMSeriesFileNames;
  auto nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetGlobalWarningDisplay(false);
  nameGenerator->SetDirectory(std::filesystem::absolute(dir).string());

  const auto &seriesUID = nameGenerator->GetSeriesUIDs();
  auto seriesItr = seriesUID.begin();
  const auto seriesEnd = seriesUID.end();

  if (seriesItr == seriesEnd) {
    std::cerr << "No DICOMs in: " << std::filesystem::absolute(dir).string()
              << "\n";
    return {};
  }

  std::vector<std::string> fileNames;
  seriesItr = seriesUID.begin();
  while (seriesItr != seriesUID.end()) {
    const auto seriesIdentifier = seriesItr->c_str();
    ++seriesItr;
    const auto names = nameGenerator->GetFileNames(seriesIdentifier);
    // For most dicom series, it would be enough to return names here in the
    // first iteration
    std::copy(names.begin(), names.end(), std::back_inserter(fileNames));
  }

  std::vector<size_t> idxtopop;
  for (size_t i = 0; i < fileNames.size(); ++i) {
    const auto modality =
        get_dcm_modality(std::string(fileNames.at(i).c_str()));
    switch (modality) {
    case DCM_MODALITY::RTIMAGE:
      break;
    default:
      idxtopop.push_back(i);
      break;
    }
  }

  std::reverse(idxtopop.begin(), idxtopop.end());
  for (const auto &id : idxtopop) {
    fileNames.erase(fileNames.begin() + id);
  }
  std::sort(fileNames.begin(), fileNames.end(),
            [](std::string &file_a, std::string &file_b) {
              const auto pos_z_a = get_dcm_offset_z(file_a);
              const auto pos_z_b = get_dcm_offset_z(file_b);
              return pos_z_a < pos_z_b;
            });

  return fileNames;
}

bool CbctRecon::ReadDicomDir(std::string &dirPath) {

  auto dir = std::filesystem::path(dirPath);
  const auto filenamelist = get_dcm_image_files(dir);

  for (auto &&filename : std::filesystem::directory_iterator(dir)) {
    if (!filename.is_regular_file ||
        !filename.is_symlink) { // I guess symlinks should be allowed?
      continue;
    }
    if (filename.path().string().find("-hash-stamp") != std::string::npos) {
      continue; // Just so the test data is less annoying.
    }
    const auto fullfilename = filename.path();
    const auto modality = get_dcm_modality(fullfilename.string());
    switch (modality) {
    case DCM_MODALITY::RTIMAGE:
      break;
    case DCM_MODALITY::RTDOSE:
      // filenamelist.push_back(fullfilename.toStdString());
      break;
    case DCM_MODALITY::RTSTRUCT:
      m_structures->set_planCT_ss(load_rtstruct(fullfilename.string()));
      m_strPathRS = fullfilename;
      break;
    case DCM_MODALITY::RTPLAN:
      break; // Maybe some pre-loading for gPMC could be useful?
    case DCM_MODALITY::RTRECORD:
      break; // I haven't ever seen one IRL
    case DCM_MODALITY::RTUNKNOWN:
      std::cerr << "File: " << fullfilename.string()
                << " was not of a recognizeable modality type!\n";
      break;
    }
  }

  ShortImageType::Pointer spShortImg;

  if (!filenamelist.empty()) {
    using dcm_reader_type = itk::ImageSeriesReader<ShortImageType>;
    auto dcm_reader = dcm_reader_type::New();
    const auto dicom_io = itk::GDCMImageIO::New();
    dcm_reader->SetImageIO(dicom_io);
    dcm_reader->SetFileNames(filenamelist);
    dcm_reader->Update();
    spShortImg = dcm_reader->GetOutput();
    m_dcm_dir = dir;
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

void ConvertUshort2Short(const UShortImageType::Pointer &spImgUshort,
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

std::filesystem::path SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
                                              std::string &strPatientID,
                                              std::string &strPatientName,
                                              std::string &strPathTargetDir) {
  if (spImg == nullptr) {
    return std::string();
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  auto plm_img = Plm_image::New(spShortImg);

  auto newDirPath =
      std::filesystem::path(strPathTargetDir + "/" + strPatientID + "_DCM");
  // std::string newDirPath =
  //  strPathTargetDir + "/" + strPatientID + strPatientName + "_DCM";

  if (!std::filesystem::exists(newDirPath)) {
    if (!std::filesystem::create_directory(newDirPath)) {
      std::cerr << "Could not create dir for dicom" << std::endl;
      return {};
    }
  }

  auto rsm = Rt_study_metadata::New();
  rsm->set_patient_id(strPatientID);
  rsm->set_patient_name(strPatientName);

  dicom_save_short(newDirPath.string(), plm_img, rsm);

  return newDirPath;
}

std::filesystem::path SaveUSHORTAsSHORT_DICOM_gdcmITK(
    UShortImageType::Pointer &spImg, std::string &strPatientID,
    std::string &strPatientName, std::string &strPathTargetDir) {
  if (spImg == nullptr) {
    return "";
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  auto newDirPath = std::filesystem::path(
      strPathTargetDir + "/" + strPatientID + strPatientName + "_DCM");

  if (!std::filesystem::exists(newDirPath)) {
    if (!std::filesystem::create_directory(newDirPath)) {
      std::cerr << "Could not create dir for gdcm dicom" << std::endl;
      return {};
    }
  } else {
    if (std::filesystem::remove_all(newDirPath) != 0) {
      if (!std::filesystem::create_directory(newDirPath)) {
        std::cerr << "Could not create dir for gdcm dicom (2)" << std::endl;
        return {};
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
  namesGenerator->SetSeriesFormat(
      (newDirPath / std::string("CT." + studyUID + ".%d.dcm")).string());

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
            << newDirPath << std::endl;
  return newDirPath;
}

template <typename T> std::string stringify(T arg) {
  return std::to_string(arg);
}

template<> std::string stringify(const std::string& arg) { return arg; }

template <char SEP, typename... Args> std::string make_sep_str(Args&&... arg) {
  return ((stringify(std::forward<Args>(arg)) + SEP) + ...);
}

std::string get_output_options(const UShortImageType::Pointer &m_spFixed) {

  // done per image because CT might be different
  // from reconstructed CBCT
  const auto str_fixed_origin =
      make_sep_str<','>(m_spFixed->GetOrigin()[0], m_spFixed->GetOrigin()[1],
                         m_spFixed->GetOrigin()[2]);

  const auto out_size = m_spFixed->GetBufferedRegion().GetSize();
  const auto str_fixed_dimension = make_sep_str<','>(
                                       out_size[0],
                                       out_size[1],
                                       out_size[2]);
  const auto out_spacing = m_spFixed->GetSpacing();
  const auto str_fixed_spacing = make_sep_str<','>(
                                     out_spacing[0],
                                     out_spacing[1],
                                     out_spacing[2]);
  const auto out_direction = m_spFixed->GetDirection();
  const auto str_fixed_direction = make_sep_str<','>(
                                       out_direction[0][0],
                                       out_direction[0][1],
                                       out_direction[0][2],
                                       out_direction[1][0],
                                       out_direction[1][1],
                                       out_direction[1][2],
                                       out_direction[2][0],
                                       out_direction[2][1],
                                       out_direction[2][2]);

  return make_sep_str<' '>(
      " --origin", str_fixed_origin, "--spacing", str_fixed_spacing,
      "--dimension", str_fixed_dimension, "--direction", str_fixed_direction);
}

bool GetCouchShiftFromINIXVI(std::string &strPathINIXVI, VEC3D *pTrans,
                             VEC3D *pRot) {
  std::filesystem::path fInfo(strPathINIXVI);
  if (!std::filesystem::exists(fInfo)) {
    return false;
  }

  std::ifstream fin;
  fin.open(strPathINIXVI);

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
    auto tmpStr = std::string(&str[0]);
    auto strListParam = split_string<'='>(tmpStr);

    std::string tagName, strVal;

    if (strListParam.size() == 2) {
      tagName = strListParam.at(0);
      strVal = strListParam.at(1);
      tagName = trim_string(tagName);
      strVal = trim_string(strVal);

      if (tagName == "CouchShiftLat") {
        couch_Lat_cm = from_string<float>(strVal).value_or(0.0f);
        bFound = true;
      } else if (tagName == "CouchShiftLong") {
        couch_Long_cm = from_string<float>(strVal).value_or(0.0f);
      } else if (tagName == "CouchShiftHeight") {
        couch_Vert_cm = from_string<float>(strVal).value_or(0.0f);
      } else if (tagName == "CouchPitch") {
        couch_Pitch = from_string<float>(strVal).value_or(0.0f);
      } else if (tagName == "CouchRoll") {
        couch_Yaw = from_string<float>(strVal).value_or(0.0f);
      } else if (tagName == "CouchYaw") {
        couch_Roll = from_string<float>(strVal).value_or(0.0f);
      }
    }
  }
  fin.close();

  if (!bFound) {
    return false;
  }

  // Warning!! dicom convention!
  pTrans->x =
      static_cast<double>(couch_Lat_cm) * 10.0; // sign should be checked
  // pTrans->y = couch_Vert_cm*10.0; //sign should be checked // IEC-->DICOM is
  // already accounted for..but sign!
  pTrans->y = static_cast<double>(couch_Vert_cm) *
              -10.0; // consistent with Tracking software
  pTrans->z =
      static_cast<double>(couch_Long_cm) * 10.0; // sign should be checked

  pRot->x = static_cast<double>(couch_Pitch);
  pRot->y = static_cast<double>(couch_Yaw);
  pRot->z = static_cast<double>(couch_Roll);
  // x,y,z: dicom
  return true;
}

bool GetXrayParamFromINI(std::string &strPathINI, float &kVp, float &mA,
                         float &ms) {
  auto info = std::filesystem::path(strPathINI);

  kVp = 0.0f;
  mA = 0.0f;
  ms = 0.0f;

  if (!std::filesystem::exists(info)) {
    return false;
  }

  // TubeMA=64.0000
  // TubeKV = 120.0000
  // TubeKVLength = 40.0000
  std::ifstream fin;
  fin.open(strPathINI);

  if (fin.fail()) {
    return false;
  }

  char str[MAX_LINE_LENGTH];

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    auto tmpStr = std::string(&str[0]);
    auto strListParam = split_string<'='>(tmpStr);

    std::string tagName;
    std::string strVal;

    if (strListParam.size() == 2) {
      tagName = strListParam.at(0);
      strVal = strListParam.at(1);
      tagName = trim_string(tagName);
      strVal = trim_string(strVal);

      if (tagName == "TubeMA") {
        mA = from_string<float>(strVal).value_or(0.0f);
      }
      if (tagName == "TubeKVLength") {
        ms = from_string<float>(strVal).value_or(0.0f);
      }
      if (tagName == "TubeKV") {
        kVp = from_string<float>(strVal).value_or(0.0f);
      }
    }
  }
  fin.close();

  return !(fabs(kVp * mA * ms) < 0.0001f);
}

// Dir or File
bool LoadShortImageDirOrFile(std::filesystem::path &strPathDir,
                             ShortImageType::Pointer &spOutputShortImg) {
  if (!std::filesystem::exists(strPathDir)) {
    return false;
  }

  Plm_image plmImg;
  plmImg.load_native(strPathDir.string());
  const auto spShortImg = plmImg.itk_short();

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
