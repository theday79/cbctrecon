#ifndef CBCTRECON_IO_H
#define CBCTRECON_IO_H
/*All IO functions used with cbctrecon*/
#include "cbctrecon_config.h"
#include "cbctrecon.h"

bool CBCTRECON_API SaveDoseGrayImage(const char *filePath, int width,
                                     int height,
                       double spacingX, double spacingY, double originLeft_mm,
                       double originTop_mm, unsigned short *pData);


template <typename ImageType> // image by value, because we call this from std::thread
void saveImageAsMHA(typename ImageType::Pointer &image) {
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(image);
  writer->SetFileName("Projections.mha");
  writer->Update();
}

void ConvertUshort2Short(UShortImageType::Pointer &spImgUshort,
  ShortImageType::Pointer &spImgShort);

QString SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
  QString &strPatientID, QString &strPatientName,
  QString &strPathTargetDir);
QString SaveUSHORTAsSHORT_DICOM_gdcmITK(UShortImageType::Pointer &spImg,
  QString &strPatientID,
  QString &strPatientName,
  QString &strPathTargetDir);

QString get_output_options(const UShortImageType::Pointer &m_spFixed);

#endif // CBCTRECON_IO_H
