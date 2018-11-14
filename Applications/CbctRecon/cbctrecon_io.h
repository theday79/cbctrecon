#ifndef CBCTRECON_IO_H
#define CBCTRECON_IO_H
/*All IO functions used with cbctrecon*/
#include "gdcmReader.h"
#include "gdcmAttribute.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"
#include "StructureSet.h"

enum DCM_MODALITY {
    RTIMAGE,
    RTDOSE,
    RTSTRUCT,
    RTPLAN,
    RTRECORD,
    RTUNKNOWN,
};

class QXmlStreamReader;

CBCTRECON_API QString MakeElektaXML(const QString &filePath_ImageDBF,
                                    const QString &filePath_FrameDBF,
                                    const QString &DICOM_UID);

CBCTRECON_API FLEXDATA XML_parseFrameForXVI5(QXmlStreamReader &xml);

CBCTRECON_API bool
LoadShortImageToUshort(QString &strPath,
                       UShortImageType::Pointer &pUshortImage);

CBCTRECON_API void ExportReconSHORT_HU(UShortImageType::Pointer &spUsImage,
                                       QString &outputFilePath);
CBCTRECON_API bool
LoadShortImageDirOrFile(QString &strPathDir,
                        ShortImageType::Pointer &spOutputShortImg);
bool CBCTRECON_API SaveDoseGrayImage(const char *filePath, int width,
                                     int height, double spacingX,
                                     double spacingY, double originLeft_mm,
                                     double originTop_mm,
                                     unsigned short *pData);

template <typename ImageType>
void saveImageAsMHA(typename ImageType::Pointer &image,
                    std::string filename = "Projections.mha") {
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(image);
  writer->SetFileName(filename);
  writer->SetUseCompression(true); // not exist in original code (rtkfdk)
  writer->Update();
}

DCM_MODALITY get_dcm_modality(QString& filename){
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
  if (modality == "RTIMAGE"){
      return RTIMAGE;
  }
  if (modality == "RTDOSE"){
      return RTDOSE;
  }
  if (modality == "RTSTRUCT"){
      return RTSTRUCT;
  }
  if (modality == "RTPLAN"){
      return RTPLAN;
  }
  if (modality == "RTRECORD"){
      return RTRECORD;
  }
  else {
      return RTUNKNOWN;
  }
}

std::unique_ptr<Rtss_modern> load_rtstruct(QString& filename){
  if (get_dcm_modality(filename) != RTSTRUCT){
     std::cerr << "This is not an RTSTRUCT: " << filename.toStdString() << "\n";
     return nullptr;
  }

  gdcm::Reader reader;
  reader.SetFileName(filename.toLocal8Bit().constData());
  if (!reader.Read())
  {
    std::cerr << "Reading dicom rtstruct: " << filename.toStdString() << " failed!\n";
    return nullptr;
  }

  gdcm::File &file = reader.GetFile();
  gdcm::DataSet &ds = file.GetDataSet();
  auto rt_struct = std::make_unique<Rtss_modern>();
  rt_struct->num_structures = 0;

  const gdcm::DataElement &roi_seq_tag = ds.GetDataElement(gdcm::Tag(0x3006, 0x0020));
  auto roi_seq = roi_seq_tag.GetValueAsSQ();
  for (auto it_roi = roi_seq->Begin(); it_roi != roi_seq->End(); ++it_roi){
    auto rt_roi = std::make_unique<Rtss_roi_modern>();
    gdcm::Attribute<0x3006, 0x0022> at_roi_number;
    at_roi_number.SetFromDataElement(it_roi->GetDataElement(at_roi_number.GetTag()));
    rt_roi->id = at_roi_number.GetValue();

    gdcm::Attribute<0x3006, 0x0026> at_roi_name;
    at_roi_name.SetFromDataElement(it_roi->GetDataElement(at_roi_name.GetTag()));
    rt_roi->name = at_roi_name.GetValue();

    rt_struct->slist.push_back(rt_roi.release());
    rt_struct->num_structures++;

  }

  const gdcm::DataElement &roi_contour_seq_tag = ds.GetDataElement(gdcm::Tag(0x3006, 0x0039));
  auto roi_contour_seq = roi_contour_seq_tag.GetValueAsSQ();
  auto i = 0U;
  for (auto it_roi_contour = roi_contour_seq->Begin(); it_roi_contour != roi_contour_seq->End(); ++it_roi_contour){
    gdcm::Attribute<0x3006, 0x0084> at_roi_contour_number;
    at_roi_contour_number.SetFromDataElement(it_roi_contour->GetDataElement(at_roi_contour_number.GetTag()));
    if (rt_struct->slist.at(i).id != at_roi_contour_number.GetValue()){
      std::cerr << "ID mismatch: " << rt_struct->slist.at(i).id << " vs " << at_roi_contour_number.GetValue() << "\n"
                << "There might be something wrong with " << rt_struct->slist.at(i).name << "\n"
                << "Caution! As we continue anyway...\n";
    }
    gdcm::Attribute<0x3006, 0x002A> at_roi_contour_colour;
    at_roi_contour_colour.SetFromDataElement(it_roi_contour->GetDataElement(at_roi_contour_colour.GetTag()));
    const auto color = at_roi_contour_colour.GetValues();
    auto s_color = std::to_string(color[0]) + " "
                 + std::to_string(color[1]) + " "
                 + std::to_string(color[2]);
    rt_struct->slist.at(i).color = s_color.c_str();

    const auto& contour_seq_tag = it_roi_contour->GetDataElement(gdcm::Tag(0x3006, 0x0040));
    auto contour_seq = contour_seq_tag.GetValueAsSQ();
    for (auto it_contour = contour_seq->Begin(); it_contour != contour_seq->End(); ++it_contour){
      auto rt_contour = std::make_unique<Rtss_contour_modern>();

      gdcm::Attribute<0x3006, 0x0048> at_contour_number;
      at_contour_number.SetFromDataElement(it_contour->GetDataElement(at_contour_number.GetTag()));

      rt_struct->slist.at(i).pslist.push_back(rt_contour.release());
    }

    i++;
  }


  return rt_struct;
}


void CBCTRECON_API ConvertUshort2Short(UShortImageType::Pointer &spImgUshort,
                                       ShortImageType::Pointer &spImgShort);

void CBCTRECON_API
ConvertShort2Ushort(ShortImageType::Pointer &spInputImgShort,
                    UShortImageType::Pointer &spOutputImgUshort);

QString CBCTRECON_API SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
                                              QString &strPatientID,
                                              QString &strPatientName,
                                              QString &strPathTargetDir);
QString CBCTRECON_API SaveUSHORTAsSHORT_DICOM_gdcmITK(
    UShortImageType::Pointer &spImg, QString &strPatientID,
    QString &strPatientName, QString &strPathTargetDir);

QString CBCTRECON_API
get_output_options(const UShortImageType::Pointer &m_spFixed);

#endif // CBCTRECON_IO_H
