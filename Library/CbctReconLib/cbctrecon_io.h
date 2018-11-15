#ifndef CBCTRECON_IO_H
#define CBCTRECON_IO_H
/*All IO functions used with cbctrecon*/
#ifdef OF
#undef OF
#endif //OF
#include "gdcmAttribute.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"
#include "StructureSet.h"


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

template<int group, int element, typename T>
auto gdcm_attribute_from(T parent) {
  gdcm::Attribute<group, element> attribute;
  attribute.SetFromDataElement(parent->GetDataElement(attribute.GetTag()));
  return attribute;
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
