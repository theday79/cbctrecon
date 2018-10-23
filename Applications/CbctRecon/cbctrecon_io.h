#ifndef CBCTRECON_IO_H
#define CBCTRECON_IO_H
/*All IO functions used with cbctrecon*/
#include "cbctrecon.h"
#include "cbctrecon_config.h"

CBCTRECON_API QString MakeElektaXML(const QString &filePath_ImageDBF,
                                    const QString &filePath_FrameDBF,
                                    const QString &DICOM_UID);

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

template <
    typename ImageType> // image by value, because we call this from std::thread
void saveImageAsMHA(typename ImageType::Pointer &image) {
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(image);
  writer->SetFileName("Projections.mha");
  writer->Update();
}

void CBCTRECON_API ConvertUshort2Short(UShortImageType::Pointer &spImgUshort,
                                       ShortImageType::Pointer &spImgShort);

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
