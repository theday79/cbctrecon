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
#include <itkDOMNodeXMLReader.h>

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
#include "free_functions.h"

namespace fs = std::filesystem;
using namespace std::literals;

namespace crl {

fs::path
MakeElektaXML(const fs::path &filePath_ImageDBF,
              const fs::path &filePath_FrameDBF,
              const std::string &DICOM_UID) {
  std::cout << "Elekta geometry XML file is being generated." << std::endl;
  // Define FRAME.DBF path
  auto reader = rtk::ElektaSynergyGeometryReader::New();
  // string strDicomUID = DICOM_UID.toLocal8Bit().constData();
  // //DICOM_UID.toStdString()  string strDicomUID = DICOM_UID.toStdString();
  // string strDbfImg = filePath_ImageDBF.toStdString();
  // string strDbfFrame = filePath_FrameDBF.toStdString();

  const auto dirPath =
      fs::absolute(filePath_ImageDBF.parent_path());

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

FLEXDATA XML_parseFrameForXVI5(itk::DOMNode::Pointer dom) {

  auto tmpResult =
      FLEXDATA{/*fGanAngle =*/0.0, /*fPanelOffsetX =*/0.0,
               /*fPanelOffsetY =*/0.0, /*bKV_On =*/true, /*bMV_On =*/false};

  /* Let's check that we're really getting a person. */
  if (dom->GetName() == "Frame") {
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
  auto idx = from_string<int>(dom->GetID()).value_or(0) + 1;
  auto next_dom = dom->Find(
      std::to_string(idx));
  /*
   * We're going to loop over the things because the order might change.
   * We'll continue the loop until we hit an EndElement named person.
   */
  while (!(next_dom == nullptr && next_dom->GetName() == "Frame")) {
    auto tmpXmlName = next_dom->GetName();
    // auto strTmpXMLName = std::string(tmpXmlName.toLocal8Bit().constData());
    // int tmpType = (int)(xml.tokenType());

    std::string tmpStr;
    /* We've found first name. */
    if (next_dom->GetName() == "Seq") {
      tmpStr = next_dom->GetTextChild()->GetText();
    }
    /* We've found surname. */
    else if (next_dom->GetName() == "DeltaMS") {
      tmpStr = next_dom->GetTextChild()->GetText();
    }
    /* We've found email. */
    else if (next_dom->GetName() == "HasPixelFactor") {
      tmpStr = next_dom->GetTextChild()->GetText();
    }
    /* We've found website. */
    else if (next_dom->GetName() == "PixelFactor") {
      tmpStr = next_dom->GetTextChild()->GetText();
    } else if (next_dom->GetName() == "GantryAngle") {
      tmpStr = next_dom->GetTextChild()->GetText();
      tmpResult.fGanAngle = from_string<float>(tmpStr).value_or(0.0f);
    } else if (next_dom->GetName() == "Exposed") {
      tmpStr = next_dom->GetTextChild()->GetText();
      if (tmpStr == "True") {
        tmpResult.bKV_On = true;
      } else {
        tmpResult.bKV_On = false;
      }
    }
    /*else if (next_dom->GetName() == "Exposed") {
        tmpStr = next_dom->GetTextChild()->GetText();
        if (tmpStr == "True")
            tmpResult.bKV_On = true;
        else
            tmpResult.bKV_On = false;
    }*/
    else if (next_dom->GetName() == "MVOn") {
      tmpStr = next_dom->GetTextChild()->GetText();
      if (tmpStr == "True") {
        tmpResult.bMV_On = true;
      } else {
        tmpResult.bMV_On = false;
      }
    } else if (next_dom->GetName() == "UCentre") {
      tmpStr = next_dom->GetTextChild()->GetText();
      tmpResult.fPanelOffsetX = from_string<float>(tmpStr).value_or(0.0f);
    } else if (next_dom->GetName() == "VCentre") {
      tmpStr = next_dom->GetTextChild()->GetText();
      tmpResult.fPanelOffsetY = from_string<float>(tmpStr).value_or(0.0f);
    } else if (next_dom->GetName() == "Inactive") {
      tmpStr = next_dom->GetTextChild()->GetText();
    }
    ++idx;
    next_dom = dom->Find(std::to_string(idx));
  }
  return tmpResult;
}

// Get the projection geometry
[[nodiscard]] rtk::ThreeDCircularProjectionGeometry::Pointer
LoadRTKGeometryFile(const fs::path &filePath) {
  auto geometryReader =
      rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(filePath.string());
  geometryReader->GenerateOutputInformation();
  std::cout << "Geometry reading succeed" << std::endl;

  auto spFullGeometry = geometryReader->GetOutputObject();

  // fullGeometry->GetGantryAngles();
  const auto angles = spFullGeometry->GetGantryAngles();
  const auto geoDataSize = angles.size(); // This is MV gantry angle!!!
  std::cout << "Geometry data size(projection gantry angles): " << geoDataSize
            << std::endl;
  if (geoDataSize < 1) {
    return nullptr;
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

  auto bScanDirectionCW = false;
  if (vTempConvAngles.at(iLowerIdx) <
      vTempConvAngles.at(iUpperIdx)) // ascending
  {
    bScanDirectionCW = true;
    std::cout << "The scan direction is CW" << std::endl;
  } else {
    std::cout << "The scan direction is CCW" << std::endl;
  }

  std::cout << "AngularGaps Size: "
            << spFullGeometry->GetAngularGaps(spFullGeometry->GetSourceAngles())
                   .size()
            << std::endl;
  return spFullGeometry;
}

bool LoadShortImageToUshort(fs::path &strPath,
                            UShortImageType::Pointer &pUshortImage) {
  using ReaderType = itk::ImageFileReader<ShortImageType>;
  auto reader = ReaderType::New();

  // std::string fileName = QFileDialog::getOpenFileName(this, "Open Image","",
  // "Plan CT file (*.mha)",0,0);

  if (strPath.empty()) {
    return false;
  }

  reader->SetFileName(strPath.string());
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

bool AlterData_RTStructureSetStorage(const fs::path &input_file,
                                     const Rtss_modern *input_rt_struct,
                                     const fs::path &output_file) {

  DRTStructureSetIOD rtstruct;
  DcmFileFormat fileformat;
  auto status =
      fileformat.loadFile(input_file.string().c_str());
  if (!status.good()) {
    std::cerr << "Could not open RT struct dcm file: "
              << input_file << "\n";
    return false;
  }
  status = rtstruct.read(*fileformat.getDataset());
  if (!status.good()) {
    std::cerr << "Could not read RT struct dcm file: "
              << input_file << "\n";
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
              << input_file << "\n";
    return false;
  }
  status =
      fileformat.saveFile(output_file.string().c_str());
  if (!status.good()) {
    std::cerr << "Could not save RT struct dcm file: "
              << output_file << "\n";
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

std::vector<std::string> get_dcm_image_files(fs::path &dir) {

  using NamesGeneratorType = itk::GDCMSeriesFileNames;
  auto nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetGlobalWarningDisplay(false);
  nameGenerator->SetDirectory(fs::absolute(dir).string());

  const auto &seriesUID = nameGenerator->GetSeriesUIDs();
  auto seriesItr = seriesUID.begin();
  const auto seriesEnd = seriesUID.end();

  if (seriesItr == seriesEnd) {
    std::cerr << "No DICOMs in: " << fs::absolute(dir).string()
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

bool ReadDicomDir(CbctRecon *p_cr, fs::path &dir) {

  const auto filenamelist = get_dcm_image_files(dir);

  for (auto &&filename : fs::directory_iterator(dir)) {
    if (!filename.is_regular_file() ||
        !filename.is_symlink()) { // I guess symlinks should be allowed?
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
      p_cr->m_structures->set_planCT_ss(load_rtstruct(fullfilename.string()));
      p_cr->m_strPathRS = fullfilename;
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
    p_cr->m_dcm_dir = dir;
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
  p_cr->m_spRefCTImg = spRescaleFilter->GetOutput();
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

fs::path
SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
                        std::string &strPatientID, std::string &strPatientName,
                        fs::path &strPathTargetDir) {
  if (spImg == nullptr) {
    return {};
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  auto plm_img = Plm_image::New(spShortImg);

  auto newDirPath = strPathTargetDir / (strPatientID + "_DCM");
  // std::string newDirPath =
  //  strPathTargetDir + "/" + strPatientID + strPatientName + "_DCM";

  if (!fs::exists(newDirPath)) {
    if (!fs::create_directory(newDirPath)) {
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

fs::path SaveUSHORTAsSHORT_DICOM_gdcmITK(
    UShortImageType::Pointer &spImg, std::string &strPatientID,
    std::string &strPatientName, fs::path &strPathTargetDir) {
  if (spImg == nullptr) {
    return {};
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  auto newDirPath = fs::path(
      strPathTargetDir / (strPatientID + strPatientName + "_DCM"));

  if (!fs::exists(newDirPath)) {
    if (!fs::create_directory(newDirPath)) {
      std::cerr << "Could not create dir for gdcm dicom" << std::endl;
      return {};
    }
  } else {
    if (fs::remove_all(newDirPath) != 0) {
      if (!fs::create_directory(newDirPath)) {
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

std::string get_output_options(const UShortImageType::Pointer &m_spFixed) {

  // done per image because CT might be different
  // from reconstructed CBCT
  const auto str_fixed_origin =
      make_sep_str<','>(m_spFixed->GetOrigin()[0], m_spFixed->GetOrigin()[1],
                        m_spFixed->GetOrigin()[2]);

  const auto out_size = m_spFixed->GetBufferedRegion().GetSize();
  const auto str_fixed_dimension =
      make_sep_str<','>(out_size[0], out_size[1], out_size[2]);
  const auto out_spacing = m_spFixed->GetSpacing();
  const auto str_fixed_spacing =
      make_sep_str<','>(out_spacing[0], out_spacing[1], out_spacing[2]);
  const auto out_direction = m_spFixed->GetDirection();
  const auto str_fixed_direction = make_sep_str<','>(
      out_direction[0][0], out_direction[0][1], out_direction[0][2],
      out_direction[1][0], out_direction[1][1], out_direction[1][2],
      out_direction[2][0], out_direction[2][1], out_direction[2][2]);

  return make_sep_str<' '>(
      " --origin"sv, str_fixed_origin, "--spacing"sv, str_fixed_spacing,
      "--dimension"sv, str_fixed_dimension, "--direction"sv, str_fixed_direction);
}

bool GetCouchShiftFromINIXVI(std::string &strPathINIXVI, VEC3D *pTrans,
                             VEC3D *pRot) {
  fs::path fInfo(strPathINIXVI);
  if (!fs::exists(fInfo)) {
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
    std::string_view tmpStr{&str[0], MAX_LINE_LENGTH};
    auto strListParam = split_string(tmpStr, "=");

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
  auto info = fs::path(strPathINI);

  kVp = 0.0f;
  mA = 0.0f;
  ms = 0.0f;

  if (!fs::exists(info)) {
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
    auto strListParam = split_string(tmpStr, "=");

    if (strListParam.size() == 2) {
      auto tagName = strListParam.at(0);
      auto strVal = strListParam.at(1);
      tagName = trim_string(tagName);
      strVal = trim_string(strVal);

      if (tagName == "TubeMA") {
        mA = from_sv(strVal, mA);
      }
      if (tagName == "TubeKVLength") {
        ms = from_sv(strVal, ms);
      }
      if (tagName == "TubeKV") {
        kVp = from_sv(strVal, kVp);
      }
    }
  }
  fin.close();

  return !(fabs(kVp * mA * ms) < 0.0001f);
}

// Dir or File
bool LoadShortImageDirOrFile(fs::path &strPathDir,
                             ShortImageType::Pointer &spOutputShortImg) {
  if (!fs::exists(strPathDir)) {
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

bool SaveDoseGrayImage(
    const fs::path &filePath, const int width, const int height,
    const double spacingX, const double spacingY, const double originLeft_mm,
    const double originTop_mm,
    unsigned short *pData) // export dose array to a specified file (16bit TIF)
{
  // Global variables
  const long m_iSubFileType = 0;
  const short m_iWidth = width;
  const short m_iHeight = height;
  const short m_iBitsPerSample = 16;
  const short m_iCompression = 1;
  const short m_iPhotometric = 0;
  const long m_iStripOffset = 1024;
  const short m_iSamplePerPixel = 1;
  const long m_iRowsPerStrip = height;
  const long m_iStripByteCnts = qRound(width * height * 2.0);

  const short m_iResolUnit = 2;
  const short m_iPgNum = 0; // or 1?
  const unsigned short m_iMinSampleVal = 0;
  const unsigned short m_iMaxSampleVal = 65535U; // old: 255
  const auto ten_mill = 10000000;
  RATIONAL m_rXResol{static_cast<long>(qRound(1 / spacingX * 25.4 *
                                              ten_mill)), // spacingX in dpi
                     ten_mill};
  RATIONAL m_rYResol{static_cast<long>(qRound(1 / spacingY * 25.4 * ten_mill)),
                     ten_mill}; // spacingY

  // double fLeftPosMM = -dataPt.x()*spacingX;
  // double fTopPosMM = dataPt.y()*spacingY;
  const auto fLeftPosMM = originLeft_mm;
  const auto fTopPosMM = -originTop_mm;

  RATIONAL m_rXPos{static_cast<long>(qRound(fLeftPosMM / 25.4 * ten_mill)),
                   ten_mill};
  RATIONAL m_rYPos{static_cast<long>(qRound(fTopPosMM / 25.4 * ten_mill)),
                   ten_mill};

  auto m_iNextOffset = 0;

  if (pData == nullptr) {
    return false;
  }

  // Set Center
  QPoint dataPt;
  dataPt.setX(qRound(m_iWidth / 2.0));
  dataPt.setY(qRound(m_iHeight / 2.0));

  {
    // FILE *fd = nullptr;
    std::ofstream fd(filePath, std::ios::binary);

    long MarkerUpper;
    long MarkerLower;

    MarkerUpper = 0x002A4949;
    MarkerLower = 0x00000008;

    fd << MarkerUpper << MarkerLower;
    // fwrite(&MarkerUpper, sizeof(long), 1, fd); // 4
    // fwrite(&MarkerLower, sizeof(long), 1, fd); // 8

    constexpr auto IFDSize = 18;

    fd << IFDSize;
    // fwrite(&IFDSize, sizeof(unsigned short), 1, fd); // 10

    std::array<TIFIFD, IFDSize> IFDarr;
    // IFDarr.reserve(IFDSize);

    int offsetX;
    auto offsetY = 0;

    const unsigned short data_type = 3;
    const auto data_cnt = 1;

    size_t i = 0;
    if (m_iSubFileType >= 0) {
      const auto tififd_tmp = TIFIFD{254, data_type, data_cnt, m_iSubFileType};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }

    if (m_iWidth >= 0) {
      const auto tififd_tmp = TIFIFD{256, data_type, data_cnt, m_iWidth};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iHeight >= 0) {
      const auto tififd_tmp = TIFIFD{257, data_type, data_cnt, m_iHeight};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iBitsPerSample >= 0) {
      const auto tififd_tmp =
          TIFIFD{258, data_type, data_cnt, m_iBitsPerSample};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iCompression >= 0) {
      const auto tififd_tmp = TIFIFD{259, data_type, data_cnt, m_iCompression};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iPhotometric >= 0) {
      const auto tififd_tmp = TIFIFD{262, data_type, data_cnt, m_iPhotometric};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iStripOffset >= 0) {
      const auto tififd_tmp =
          TIFIFD{273, 4, data_cnt, static_cast<int>(m_iStripOffset)};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iSamplePerPixel >= 0) {
      const auto tififd_tmp =
          TIFIFD{277, data_type, data_cnt, m_iSamplePerPixel};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iRowsPerStrip >= 0) {
      const auto tififd_tmp =
          TIFIFD{278, data_type, data_cnt, static_cast<int>(m_iRowsPerStrip)};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iStripByteCnts >= 0) {
      const auto tififd_tmp =
          TIFIFD{279, 4, data_cnt, static_cast<int>(m_iStripByteCnts)};
      IFDarr.at(i) = tififd_tmp;
      ++i;
      /*if (m_iSamplePerPixel == 1)
      dataVal = m_iStripByteCnts;
      else if (m_iSamplePerPixel == 3)
      dataVal = (int)(m_iStripByteCnts/3.0);
      */
    }
    if (m_rXResol.a != 0) {
      offsetX = 8 + 2 + 12 * IFDSize + 4;

      const auto tififd_tmp = TIFIFD{282, 5, data_cnt, offsetX}; // maximum
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_rYResol.a != 0) {
      offsetY = 8 + 2 + 12 * IFDSize + 4 + 8;

      const auto tififd_tmp = TIFIFD{283, 5, data_cnt, offsetY};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }

    if (m_rXPos.a != 0) {
      offsetX = 8 + 2 + 12 * IFDSize + 4 + 8 + 8;

      const auto tififd_tmp = TIFIFD{286, 5, data_cnt, offsetX}; // maximum
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_rYPos.a != 0) {
      offsetY = 8 + 2 + 12 * IFDSize + 4 + 8 + 8 + 8;

      const auto tififd_tmp = TIFIFD{287, 5, data_cnt, offsetY};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }

    ////// Do not insert if dataVal is initial value
    if (m_iMinSampleVal >= 0) {
      const auto tififd_tmp = TIFIFD{280, data_type, data_cnt, m_iMinSampleVal};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iMaxSampleVal >= 0) {
      const auto tififd_tmp = TIFIFD{281, data_type, data_cnt, m_iMaxSampleVal};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iResolUnit >= 0) {
      const auto tififd_tmp = TIFIFD{296, data_type, data_cnt, m_iResolUnit};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    if (m_iPgNum >= 0) {
      const auto tififd_tmp = TIFIFD{297, data_type, 2, m_iPgNum};
      IFDarr.at(i) = tififd_tmp;
      ++i;
    }
    for (auto &&it : IFDarr) {
      fd.write(reinterpret_cast<char *>(&it), sizeof(TIFIFD));
      // fwrite(&it, sizeof(TIFIFD), 1, fd);
    }
    fd << m_iNextOffset;
    // fwrite(&m_iNextOffset, 4, 1, fd);

    fd.write(reinterpret_cast<char *>(&m_rXResol), sizeof(RATIONAL));
    fd.write(reinterpret_cast<char *>(&m_rYResol), sizeof(RATIONAL));
    // fwrite(&m_rXResol, 8, 1, fd);
    // fwrite(&m_rYResol, 8, 1, fd);

    fd.write(reinterpret_cast<char *>(&m_rYPos), sizeof(RATIONAL));
    fd.write(reinterpret_cast<char *>(&m_rXPos), sizeof(RATIONAL));
    // fwrite(&m_rXPos, 8, 1, fd); // Used to be 10 instead of 1, but that
    // must've been a mistake fwrite(&m_rYPos, 8, 1, fd);

    const auto iDummySize = static_cast<size_t>(1024ull - (offsetY + 8));

    auto tmpDummy = std::valarray<char>(static_cast<char>(0), iDummySize);
    fd.write(&tmpDummy[0], sizeof(char) * iDummySize);
    // fwrite(&tmpDummy[0], sizeof(char), iDummySize, fd); // "Padded with zeros
    // until'"

    const auto imgSize = m_iWidth * m_iHeight;

    for (auto i = 0; i < imgSize; i++) {
      // fwrite(&pData[i], 2, 1, fd);
      fd.write(reinterpret_cast<char *>(&pData[i]), sizeof(pData[i]));
    }

    // fclose(fd);
  }

  return true;
}

} // namespace crl